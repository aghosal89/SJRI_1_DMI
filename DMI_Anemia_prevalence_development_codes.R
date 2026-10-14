
#--------------------------------------------------------------------------------
# This script contains the codes that develops DMI from the dataset triangulated 
# from the NFHS-4, NSSO (68th round), Satellite-PM2.5 exposure datasets
#--------------------------------------------------------------------------------

# clear the memory first
rm(list = ls())
gc()

# Set the working directory
setwd("~/Documents/Aritra_sjri_files")

# read the coding libraries necessary for our analysis
library(dplyr)
library(ggplot2)
library(gridExtra)
library(writexl)
library(glmnet)
library(patchwork)
library(mgcv)
library(ggeffects)
library(corrplot)
library(car)

# data reading and cleaning
d4c1<- readRDS("nfhs4_nss68 (2).rds") # NFHS-4 data 
d4c2 <- readRDS("~/Documents/Aritra_sjri_files/NSSO_NFHS_nutrients_pc.rds")

d4c1<- d4c1  %>%
  # select columns
  select(caseid, psu, StateNFHS,
    PM2.5_Child, sex_child, WealthIndex, TypeOfLocality,
    Cooking_fuel, hh_member_smoke_frequency,
    Education, hb_child, age_child)

# select the per-capita nutrients, no missing observations were found in the data
d4c2<- d4c2 %>% 
  select(NFHS_ID, vitaminA = pc_VitaminA.RAE, 
         vitaminB12 = pc_VitaminB12, vitaminC = pc_TotalAscorbicAcid.VITC, 
         vitaminD= pc_VitaminD.D2.D3, vitaminE= pc_alphaTocopherolEquivalent.VITE, 
         selenium = pc_SeleniumSe.SE, zinc = pc_ZincZn.ZN, iron = pc_IronFe.FE)

d4c2$NFHS_ID <- as.character(d4c2$NFHS_ID)

# merge the datasets based on the column 'caseid'
d4c<- merge(d4c1, d4c2, by.x ="caseid", by.y= "NFHS_ID", all.x = T)

# select children in the age range 6-59 months
d4c<- d4c %>% filter(age_child >= 6) #221858, 21

# gather the names of the categorical variables in the model
cat_vars<- c("sex_child", "WealthIndex", "TypeOfLocality", "Cooking_fuel",
              "hh_member_smoke_frequency", "Education")

# gather the nutrients in the dataset
nutrients <- c("vitaminA", "vitaminB12", "vitaminC", "vitaminD", "vitaminE" ,
               "selenium", "zinc")

# gather the names of the numerical variables in the model
num_vars<- c(nutrients, "iron", "PM2.5_Child")

# to check if there are any missing data
anyNA(d4c) #TRUE

#--------------------------------------------------------------------------------
# Before further analysis, percentile trimming is performed removing the outliers 
# below 0.5% and above 99.5% for every variable in the model. Each variable from
# 'num_vars' is considered here.
#--------------------------------------------------------------------------------

d4c_cleaned<- d4c
keep<- rep(TRUE, nrow(d4c_cleaned))
for(v in c(num_vars)){
   q <- quantile(d4c_cleaned[[v]], c(0.005, 0.995), na.rm = TRUE)
  keep <- keep & d4c_cleaned[[v]] >= q[1] & d4c_cleaned[[v]] <= q[2]}
d4c_cleaned <- d4c_cleaned[keep, ]

# dimension of the trimmed data
dim(d4c_cleaned) #206574, 20

# counts of missingness -------------------------
all_vars<- c(num_vars, cat_vars) 

missing_summary <- data.frame(Variable = all_vars,
  Missing = sapply(d4c_cleaned[all_vars], \(x) sum(is.na(x))),
  Non_missing = sapply(d4c_cleaned[all_vars], \(x) sum(!is.na(x)))
)

missing_summary$Missing_pct<- round(100 * missing_summary$Missing 
                                    / nrow(d4c_cleaned), 2); missing_summary

# Since data missingness is problematic, we remove those rows
d4c_cleaned <- na.omit(d4c_cleaned)
dim(d4c_cleaned)  #177463, 20

# summary of the data
me<- apply(d4c_cleaned[, nutrients], 2, median)
iqr<- apply(d4c_cleaned[, nutrients], 2, IQR)
mn<- apply(d4c_cleaned[, nutrients], 2, mean)
std<- apply(d4c_cleaned[, nutrients], 2, sd)

as.data.frame(bind_rows(mn, std, me, iqr))

#-------------------------------------------------------------------------------
# codes to create the contents of the table 2 in the document
#-------------------------------------------------------------------------------

tab2 <- lapply(cat_vars, function(v) {
  tbl <- table(d4c_cleaned[[v]])
  pct <- prop.table(tbl) * 100
  data.frame(
    Frequency  = as.numeric(tbl),
    Percentage = round(as.numeric(pct), 1),
    row.names  = names(tbl)
  )
})
names(tab2) <- cat_vars

# add variable names
catvar_labels <- c(
  sex_child      = "Sex of Child",
  WealthIndex    = "Wealth Index",
  TypeOfLocality = "Type of Locality",
  Cooking_fuel   = "Cooking Fuel",
  hh_member_smoke_frequency  = "Smoking per household",
  Education      = "Education"
)

tab2_combined <- do.call(rbind, lapply(names(tab2), function(var) {
  df <- tab2[[var]]
  data.frame(Variable = c(catvar_labels[var], rep("", nrow(df) - 1)),
             Category   = rownames(df), Frequency  = df$Frequency,
             Percentage = df$Percentage)}
)
)

# save the file in excel file 
tab2_combined$`N(%)`<- paste(tab2_combined$Frequency, " (", 
                              tab2_combined$Percentage, ")", sep ="")
tab2_combined
write.csv(tab2_combined, "Table2_nfhs4_ChU5.csv")

# We save the categorical variables as factors for further analysis
d4c_cleaned <- d4c_cleaned %>%
  mutate(across(all_of(cat_vars), as.factor)) %>% 
  mutate(across(all_of(nutrients), ~as.numeric(scale(.x)), .names= "scaled_{.col}"))

#--------------------------------------------------------------------------------
# Bootstrap and Ridge Regression, we perform scaling, bootstrapping and Fitting 
# of the ridge model (adjusting for scaled nutrients)
#
# The function 'relWeight' computes the relative weights of the micronutrients 
# that are used in constructing the Diet Moderation Index (DMI).
#--------------------------------------------------------------------------------

# Inputs: 

# data:  The dataset containing all the necessary variables in the study contained 
#         in the variables set.
# nsim:  The number of simulations to run.
# NRS :  The number of the rows to be selected from the data above for analysis.

# Output:

# Let m be the number of antioxidant nutrients considered in the DMI, then a matrix 
# of order nsim x m is the output, whose each row represents the relative weights 
# of the micro-nutrients considered in the model, i.e. Vitamins A, B12, C, D, E,
# Selenium, Zinc (in an order specified in the data).

relWeight<- function(data, nsim, NRS) {
  
  interaction_terms<- c(
    "PM2.5_Child:scaled_vitaminA",
    "PM2.5_Child:scaled_vitaminB12",
    "PM2.5_Child:scaled_vitaminC",
    "PM2.5_Child:scaled_vitaminD",
    "PM2.5_Child:scaled_vitaminE",
    "PM2.5_Child:scaled_selenium",
    "PM2.5_Child:scaled_zinc"
  )
  
  # Pre-allocate matrix
  est <- matrix(NA, nrow= nsim, ncol= length(interaction_terms))
  
  colnames(est) <- interaction_terms
  
  # Simulation indices for bootstrap
  sim_data <- replicate(nsim, sample(seq_len(nrow(data)), NRS, replace= FALSE))
  
  for(i in 1:nsim){
    
    temp <- sim_data[, i]
    mydata <- data[temp, ]
    
    # Scale nutrients
    mydata <- mydata %>% mutate(across(all_of(nutrients), ~as.numeric(scale(.x)),
          .names = "scaled_{.col}"))
    
    # Design matrix
    X <- model.matrix(hb_child ~ PM2.5_Child* (scaled_vitaminA + scaled_vitaminB12 +
          scaled_vitaminC + scaled_vitaminD + scaled_vitaminE + scaled_selenium +
          scaled_zinc ) + sex_child + WealthIndex + TypeOfLocality + Cooking_fuel + 
        Education + iron, data = mydata)[,-1]
    
    Y <- mydata$hb_child
    
    # Ridge CV
    cv_ridge <- cv.glmnet(X, Y, alpha = 0)
    ridge_final <- glmnet(X, Y, alpha = 0, lambda = cv_ridge$lambda.min)
    coef_vec <- as.matrix(coef(ridge_final))
    coef_names <- rownames(coef_vec)
    
    matched <- match(interaction_terms, coef_names)
    est[i, ] <- coef_vec[matched, 1]
  }
  as.data.frame(est)
}  

#-------------------------------------------------------------------------------
# Function: for the vector x, returns a table summarizing the counts, means and 
# standard deviations of the positive and the negative entries in x. 
#-------------------------------------------------------------------------------

score <- function(x) {
  
  pos<- x[x > 0]  # collect positive entries in x
  neg<- x[x < 0]  # collect negative entries in x
  
  data.frame(# length, mean and standard deviation of positive items
    p.n= length(pos), p.mean= mean(pos, na.rm= TRUE), 
             p.sd= sd(pos, na.rm = TRUE),
    
             # length, mean and standard deviation of negative items
             n.n = length(neg), n.mean = mean(neg, na.rm = TRUE),
             n.sd = sd(neg, na.rm = TRUE))
}

#--------------------------------------------------------------------------------
# Here we provide rationale for the choice of sample size 20000 for ridge 
# regression estimation through a calculation for different sample sizes done 
# over multiple simulation runs.
#--------------------------------------------------------------------------------

set.seed(1571)

# choice of sample sizes 
nrs <- c(1000, 6000, 11000, 16000, 21000, 26000, 31000, 36000)
nsim<- 200  # number of simualation runs

output_list<- vector("list", length(nrs)) # a list to store the output

for(k in seq_along(nrs)){
  i <- nrs[k]
  cat("=")
  temp <- relWeight(data= d4c_cleaned, nsim = nsim, NRS = i)
  output_list[[k]] <- temp
}

weight_matrix<- lapply(output_list, function(temp){
    
    apply(temp, 2, function(x) {
        pos <- x[x > 0]
        if(length(pos) == 0){
          return(0)}
        p.mean <- mean(pos, na.rm= TRUE)
        positive_prop <- length(pos)/length(x)
        p.mean * positive_prop
      })
  })

weight_matrix<- do.call(rbind, weight_matrix)
weight_matrix<- as.data.frame(weight_matrix)

weight_matrix_std <- cbind.data.frame(nrs, 100* sweep(weight_matrix, 1, 
                                        rowSums(weight_matrix), "/"))
colnames(weight_matrix_std)<- c("NRS", "VitaminA", "VitaminB12", "VitaminC",
                                "VitaminD", "VitaminE", "Selenium", "Zinc")

write.csv(weight_matrix_std, "Percentage_weights_DMI(Anemia_Prevalence).csv")

#--------------------------------------------------------------------------------
# Finally, we compute the vector of weights for the optimal choice of the sample   
# size 20,000 in the bootstrapped ridge-regression framework, which is to be used 
# for the validation purposes in the script 'DMI_Anemia_prevalence_validation_codes.R'
#--------------------------------------------------------------------------------

set.seed(6578)

# we compute the relative weights based on the best ridge regression model
est<- relWeight(data= d4c_cleaned, nsim= 200, NRS= 20000)

# Apply the scores on the estimates of the regression
tab <- do.call(rbind, apply(est, 2, score))
tab <- as.data.frame(tab)
rownames(tab) <- c("Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D", "Vitamin E",
                   "Se", "Zn")

# Positive proportion
tab$positive_prop<- tab$p.n/ nsim

# Weight calculation in percentage form
wt<- 100* (tab$p.mean* tab$positive_prop)/ sum(tab$p.mean* tab$positive_prop,
                                                na.rm = TRUE)
tab$weight<- wt
names(wt)<- nutrients

# Export the table
write_xlsx(as.data.frame(tab), "Tablehb_wts.xlsx")
 
# export the relative weights for the next stage of the analysis
write.csv(wt, "DMI_ChU5.csv")

# finally write the NFHS-4 data for women for further analysis
write.csv(d4c_cleaned, "NFHS4_cleandata_children.csv")

