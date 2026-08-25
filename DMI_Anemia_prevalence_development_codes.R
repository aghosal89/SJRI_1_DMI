
#--------------------------------------------------------------------------------
# This script contains the codes that develops DMI from the dataset triangulated 
# from the NFHS-4, NSSO (68th round), Satellite-PM2.5 exposure datasets
#--------------------------------------------------------------------------------

# clear the memory first
rm(list = ls())

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
d4c2<- readRDS("NSSO_nutrients4.rds") # NSSO data 

d4c1<- d4c1 %>%
  # select children aged 6-59 months
  filter(age_child >= 6) %>%  
  # select columns
  select(caseid, psu, StateNFHS,
    PM2.5_Child, sex_child, WealthIndex, TypeOfLocality,
    Cooking_fuel, hh_member_smoke_frequency,
    Education, hb_child, anemia_child)

# select the per-capita nutrients
d4c2<- d4c2 %>% 
  select(caseid= NFHS_ID, vitaminA = pc_VitaminA.RAE, 
         vitaminB12 = pc_VitaminB12, vitaminC = pc_TotalAscorbicAcid.VITC, 
         vitaminD= pc_VitaminD.D2.D3, vitaminE= pc_alphaTocopherolEquivalent.VITE, 
         selenium = pc_SeleniumSe.SE, zinc = pc_ZincZn.ZN, iron = pc_IronFe.FE)

d4c2$caseid <- as.character(d4c2$caseid)

# merge the datasets based on the column 'caseid'
d4c<- inner_join(d4c2, d4c1, by= "caseid")

# gather the names of the categorical variables in the model
cat_vars<- c("sex_child", "WealthIndex", "TypeOfLocality", "Cooking_fuel",
              "hh_member_smoke_frequency", "Education")

# gather the nutrients in the dataset
nutrients <- c("vitaminA", "vitaminB12", "vitaminC", "vitaminD", "vitaminE" ,
               "selenium", "zinc")

# gather the names of the numerical variables in the model
num_vars<- c(nutrients, "iron", "PM2.5_Child","hb_child")

# to check if there are any missing data
anyNA(d4c) #TRUE

#--------------------------------------------------------------------------------
# Before further analysis, percentile trimming is performed removing the outliers 
# below 0.5% and above 99.5% for every variable in the model. Each variable from
# 'num_vars' is considered here.
#--------------------------------------------------------------------------------

clean_trimmed<- d4c

keep<- rep(TRUE, nrow(clean_trimmed))
for(v in c(nutrients, "iron")){
  
  q <- quantile(
    clean_trimmed[[v]],
    c(0.005, 0.995),
    na.rm = TRUE
  )
  
  keep <- keep &
    clean_trimmed[[v]] >= q[1] &
    clean_trimmed[[v]] <= q[2]
}

clean_trimmed <- clean_trimmed[keep, ]

# dimension of the trimmed data
dim(clean_trimmed) #50325, 20

# Since data missingness is problematic, we remove those rows
clean_trimmed <- na.omit(clean_trimmed)
dim(clean_trimmed) 

me<- apply(clean_trimmed[, nutrients], 2, median)
iqr<- apply(clean_trimmed[, nutrients], 2, IQR)
mn<- apply(clean_trimmed[, nutrients], 2, mean)
std<- apply(clean_trimmed[, nutrients], 2, sd)

as.data.frame(bind_rows(me, iqr, mn, std))

#-------------------------------------------------------------
# codes to create the contents of the table 1 in the document
#-------------------------------------------------------------

# descriptive statistics for the numerical variables in the model
tab1<- sapply(clean_trimmed[, num_vars],
            function(x) c(N = sum(!is.na(x)),
                          Mean = mean(x, na.rm=TRUE),
                          Median=median(x,na.rm=TRUE),
                          Sd = sd(x, na.rm=TRUE),
                          Q1     = quantile(x, 0.25, na.rm = TRUE),
                          Q3     = quantile(x, 0.75, na.rm = TRUE),
                          Min =  min(x, na.rm=TRUE),
                          Max =  max(x, na.rm=TRUE)))

tab1=as.data.frame(t(tab1))

# add variable names
numvar_labels <- c(vitaminA = "Vitamin A", vitaminB12 = "Vitamin B",
  vitaminC = "Vitamin C", vitaminD = "Vitamin D",
  selenium = "Selenium", zinc = "Zinc", iron = "Iron", vitaminE = "Vitamin E", 
  PM2.5_Child = "PM2.5", hb_child = "Hb"
)

rownames(tab1) <- numvar_labels
write.csv(tab1,"Table1_nfhs4.csv")

#-------------------------------------------------------------------------------
# codes to create the contents of the table 2 in the document
#-------------------------------------------------------------------------------

tab2 <- lapply(cat_vars, function(v) {
  tbl <- table(clean_trimmed[[v]])
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
write.csv(tab2_combined, "Table2_nfhs4.csv")


# we save this data 'clean_trimmed' for next stages of analysis
#write.csv(clean_trimmed, "NFHS4_cleanData.csv")

# We save the categorical variables as factors for further analysis
clean_trimmed <- clean_trimmed %>%
  mutate(across(all_of(cat_vars), as.factor))

write.csv(clean_trimmed, "NFHS4_cleandata_children.csv")

#--------------------------------------------------------------------------------
# The following function is used to create the contents of the tables in the 
# document. 
#--------------------------------------------------------------------------------

# Inputs:

# 1. data = The triangulated NFHS-4 data on which the model is built.
# 2. nv = nutrient variable
# 3. ov = outcome variable which is Haemoglobin of the observed child

# Output: A table.

analyze_quartiles <- function(data, nv, ov = "hb_child") {
  
  # Remove missing exposure values
  temp_data <- data %>%
    filter(!is.na(.data[[nv]]))
  
  # Quantile cut points
  qdata <- quantile(
    temp_data[[nv]],
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )
  
  qdata_round <- round(qdata, 3)
  
  # Create quartiles
  temp_data <- temp_data %>%
    mutate(quartile = ntile(.data[[nv]], 4))
  
  # Quartile labels
  quartile_labels <- c(
    paste0("[", qdata_round[1], ", ", qdata_round[2], ")"),
    paste0("[", qdata_round[2], ", ", qdata_round[3], ")"),
    paste0("[", qdata_round[3], ", ", qdata_round[4], ")"),
    paste0("[", qdata_round[4], ", ", qdata_round[5], "]")
  )
  
  temp_data$quartile <- factor(temp_data$quartile, levels = 1:4, 
                               labels = quartile_labels)
  
  # Summary table
  result_table <- temp_data %>%
    group_by(quartile) %>%
    summarise(
      N = n(),
      Mean_Hb = mean(.data[[ov]], na.rm = TRUE),
      SD_Hb = sd(.data[[ov]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # ANOVA
  formula_obj <- as.formula(paste(ov, "~ quartile"))
  anova_test <- tryCatch(
    aov(formula_obj, data = temp_data),
    error = function(e) NULL
  )
  
  p_value <- NA
  if (!is.null(anova_test)) {
    p_value <- summary(anova_test)[[1]][["Pr(>F)"]][1]
  }
  
  result_table$P_value <- ""
  result_table$P_value[1] <- ifelse(
    is.na(p_value),
    NA,
    ifelse(p_value < 0.001, "<0.001", round(p_value, 3))
  )
  return(result_table)
}


# Run analysis
results_list <- lapply(nutrients, function(nutrient) {
  analyze_quartiles(clean_trimmed, nutrient)
})

names(results_list) <- numvar_labels[nutrients]

# Combine results
combined_results <- bind_rows(results_list, .id = "Exposure")

# Export
write_xlsx(combined_results, "Hb_vs_Quartiles_Table.xlsx")


# similar table for the categorical predictors in 'cat_vars'
tab_final <- do.call(rbind, lapply(cat_vars, function(v) {
  clean_trimmed %>%
    group_by(.data[[v]]) %>%
    summarise(
      N = n(),
      Mean_Hb = mean(hb_child, na.rm = TRUE),
      SD_Hb = sd(hb_child, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Total_N = sum(N),   # <-- total added beside each category
      Percent = round((N / Total_N) * 100, 1),
      Variable = v,
      Category = .data[[v]]) %>%
    select(Variable, Category, Percent, Mean_Hb, SD_Hb)})
)

pvals <- sapply(cat_vars, function(v) {
  summary(aov(hb_child ~ clean_trimmed[[v]], data = clean_trimmed))[[1]][["Pr(>F)"]][1]
})

pvals_df <- data.frame(
  Variable = names(pvals),
  P_value = pvals
)

tab_final <- merge(tab_final, pvals_df, by = "Variable", all.x = TRUE)
tab_final$Variable <- catvar_labels[tab_final$Variable]
tab_final


tab_final <- tab_final %>%
  group_by(Variable) %>%
  mutate(
    row_id = row_number(),
    #Total_N = ifelse(row_id == 1, Total_N, NA),
    P_value = ifelse(row_id == 1, P_value, NA)
  ) %>%
  ungroup() %>%
  mutate(
    Variable = ifelse(row_id == 1, Variable, ""),
    #Total_N  = ifelse(is.na(Total_N), "", as.character(Total_N)),
    P_value  = ifelse(is.na(P_value), "", as.character(P_value))
  ) %>%
  select(-row_id)
tab_final

writexl::write_xlsx(tab_final, "Categorical_association_test.xlsx")

#--------------------------------------------------------------------------------
# Bootstrap and Ridge Regression, we perform scaling, bootstrapping and Fitting 
# of the ridge model (adjusting for scaled nutrients)

# The function 'relWeight' computes the relative weights of the micronutrients 
# that are used in constructing the Diet Moderation Index (DMI).
#--------------------------------------------------------------------------------

# Inputs: 

# data =  the dataset containing all the necessary variables in the study contained 
#         in the variables set.
# nsim =  The number of simulations to be run
# NRS =   the number of the rows to be selected from 'data' for analysis

# Output:

# A matrix of order nsim x 7, whose each row represents the relative weights of 
# the micro-nutrients considered in the model, i.e. Vitamins A, B12, C, D, E,
# Selenium, Zinc (in an order specified in the data).

relWeight <- function(data, nsim, NRS){
  
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
  est <- matrix(
    NA,
    nrow = nsim,
    ncol = length(interaction_terms)
  )
  
  colnames(est) <- interaction_terms
  
  # Simulation indices for bootstrap
  sim_data <- replicate(nsim, sample(seq_len(nrow(data)), NRS, replace= TRUE))
  
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

# Summary function: returns a table summarizing the positive and 
# the negative items in the vector x
score <- function(x) {
  
  pos <- x[x > 0]  # number of positive items in x
  neg <- x[x < 0]  # number of negative items in x
  
  data.frame(
    # length, mean and standard deviation of positive items
    p.n= length(pos),
    p.mean= mean(pos, na.rm= TRUE),
    p.sd= sd(pos, na.rm = TRUE),
    
    # length, mean and standard deviation of negative items
    n.n = length(neg),
    n.mean = mean(neg, na.rm = TRUE),
    n.sd = sd(neg, na.rm = TRUE))
}


#--------------------------------------------------------------------------------
# Here we provide rationale for the choice of sample size 20000 for ridge 
# estimate a calculation for different sample sizes is given below
#--------------------------------------------------------------------------------

set.seed(1571)
nrs <- c(1000, 6000, 11000, 16000, 21000, 26000, 31000, 36000)
nsim<- 200

output_list <- vector("list", length(nrs))

for(k in seq_along(nrs)){
  
  i <- nrs[k]
  cat("=")
  temp <- relWeight(data= clean_trimmed, nsim = nsim, NRS = i)
  output_list[[k]] <- temp
}


output <- dplyr::bind_rows(output_list)
colnames(output)<- c('Vitamin A', 'Vitamin B', 'Vitamin C', 'Vitamin D', 'Vitamin E',
                     'Se', 'Zn')
output$NRS<- as.factor(rep(nrs, each= nsim))

writexl::write_xlsx(output, "Tableout1.xlsx") 

output_list[[8]]
weight_matrix<- lapply(output_list, function(temp){
    
    apply(temp, 2, function(x) {
        pos <- x[x > 0]
        if(length(pos) == 0){
          return(0)
        }
        p.mean <- mean(pos, na.rm = TRUE)
        positive_prop <- length(pos)/length(x)
        p.mean * positive_prop
      })
  }
)

weight_matrix <- do.call(rbind, weight_matrix)
weight_matrix <- as.data.frame(weight_matrix)

weight_matrix_std <- cbind.data.frame(nrs, 100* sweep(weight_matrix, 1, 
                                        rowSums(weight_matrix), "/"))
colnames(weight_matrix_std)<- c("NRS", "VitaminA", "VitaminB12", "VitaminC",
                                "VitaminD", "VitaminE", "Selenium", "Zinc")

write.csv(weight_matrix_std, "Percentage_weights_DMI(Anemia_Prevalence).csv")

#--------------------------------------------------------------------------------
# Finally, we compute the weights vector for the choice of sample size 20000 for 
# ridge regression, which is to be used for the validation purposes in the script
# 'DMI_Anemia_prevalence_validation_codes.R'
#--------------------------------------------------------------------------------

set.seed(6578)

# we compute the relative weights based on the best ridge regression model
est <- relWeight(data= clean_trimmed, nsim= 200, NRS= 20000)

# Apply the scores on the estimates of the regression
tab <- do.call(rbind, apply(est, 2, score))
tab <- as.data.frame(tab)
rownames(tab) <- c("Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D", "Vitamin E",
                   "Se", "Zn")

# Positive proportion
tab$positive_prop<- tab$p.n/ nsim

# Weight calculation in percentage form
wt <- 100* (tab$p.mean* tab$positive_prop)/ sum(tab$p.mean* tab$positive_prop,
                                                na.rm = TRUE)
tab$weight<- wt
names(wt)<- c("Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D", "Vitamin E",
              "Se", "Zn")

# Export the table
write_xlsx(as.data.frame(tab), "Tablehb_wts.xlsx")
 
# export the relative weights for the next stage of the analysis
write.csv(wt, "Relative_weights_of_nutrients.csv")
