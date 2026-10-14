
#-------------------------------------------------------------------------------
# This script contains the codes that validates the DMI obtained from the 
# triangulated NFHS-4 dataset on the triangulated NFHS-5 dataset.
#-------------------------------------------------------------------------------

# clear the memory first
rm(list= ls())
gc()

# Set the working directory
setwd("~/Documents/Aritra_sjri_files")

# read the coding libraries necessary for our analysis
library(survey)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(readxl)
library(writexl)
library(latex2exp)
library(haven)

d5c_smoke<- read_dta("IAHR7AFL.DTA", col_select= c("hhid", "hv252")) %>%
  transmute(hhid= hhid, hh_frequency_smoke= case_when(hv252 == 0  ~ "never", 
                                       hv252 == 1  ~ "daily",
                                       hv252 == 2  ~ "weekly",
                                       hv252 == 3  ~ "monthly",
                                       hv252 == 4  ~ "less than once a month"))
dim(d5c_smoke) #636699, 2

# select the variables corresponding to the children
d5c<- read_dta("IAKR7AFL.DTA", col_select= c("caseid",  #women id
          "v021",  #psu
          "v024",    # state
          "v005",   #sample Weight
          "b4",    #sex
          "b8",   #current age
          "hw1",  # age in month
          "hw2",  # weight
          "hw3",  # Height
          "hw70", # haz
          "hw71",  # waz
          "hw72", # whz
          "hw73",  # BMZ
          "hw56",  # HB adjusted for altitude
          "hw57",  # Anemia level
          "v190",  # wealth index
          "v025",  # type of residence: urban, rural
          "v161",   # type of cooking fuel
          "v022",   # strata variable
          "v023",   # strata variable (alternative)
          "v106"   # women education
          #"s719"    # second-hand smoke
))

d5c<- d5c %>%
  transmute(cooking_fuel= ifelse(v161 %in% c(1, 2, 3, 4), "Clean", 
                                 ifelse(v161 == 5, "Kerosene",
                                        ifelse(v161 %in% c(6, 7, 8, 9, 10, 11), 
                                               "Unclean", NA))),
                
    ## Child sex
    child_gender = factor(b4, levels = c(1, 2), labels = c("Male", "Female")),
                
    # Age of the child in months
    child_age= hw1,
    
    ## Wealth quintile
    wealthindex = factor(v190, levels = 1:5, labels = c( "Poorest", "Poorer", 
                                                    "Middle", "Richer", "Richest")),
                
    ## Residence
    residence_type= factor(v025, levels = c(1, 2), labels = c("Urban", "Rural")),
                
    ## Mother's education
    women_education = factor(v106, levels = 0:3, labels = c("No education", "Primary",
                                                            "Secondary","Higher")),
                
    ## Child anaemia level
    anemia_levelc = case_when(hw57 == 4 ~ 0, hw57 %in% c(1, 2, 3) ~ 1, 
                              TRUE ~ NA_real_),
    
    caseid= caseid,  #women id
    psu= v021,  #psu
    sample_weight= v005/1000000,   #sample Weight
    strata= v023, 
    #sh_smoke = case_when(s719 == 1 ~ "yes", s719 == 0 ~ "no", TRUE ~ NA_character_),
    
    state = case_when(
      v024 == 1  ~ "Jammu & Kashmir",
      v024 == 2  ~ "Himachal Pradesh",
      v024 == 3  ~ "Punjab",
      v024 == 4  ~ "Chandigarh",
      v024 == 5  ~ "Uttarakhand",
      v024 == 6  ~ "Haryana",
      v024 == 7  ~ "NCT of Delhi",
      v024 == 8  ~ "Rajasthan",
      v024 == 9  ~ "Uttar Pradesh",
      v024 == 10 ~ "Bihar",
      v024 == 11 ~ "Sikkim",
      v024 == 12 ~ "Arunachal Pradesh",
      v024 == 13 ~ "Nagaland",
      v024 == 14 ~ "Manipur",
      v024 == 15 ~ "Mizoram",
      v024 == 16 ~ "Tripura",
      v024 == 17 ~ "Meghalaya",
      v024 == 18 ~ "Assam",
      v024 == 19 ~ "West Bengal",
      v024 == 20 ~ "Jharkhand",
      v024 == 21 ~ "Odisha",
      v024 == 22 ~ "Chhattisgarh",
      v024 == 23 ~ "Madhya Pradesh",
      v024 == 24 ~ "Gujarat",
      v024 == 25 ~ "Dadra & Nagar Haveli and Daman & Diu",
      v024 == 27 ~ "Maharashtra",
      v024 == 28 ~ "Andhra Pradesh",
      v024 == 29 ~ "Karnataka",
      v024 == 30 ~ "Goa",
      v024 == 31 ~ "Lakshadweep",
      v024 == 32 ~ "Kerala",
      v024 == 33 ~ "Tamil Nadu",
      v024 == 34 ~ "Puducherry",
      v024 == 35 ~ "Andaman & Nicobar Islands",
      v024 == 36 ~ "Telangana",
      v024 == 37 
      ~ "Ladakh",
      TRUE ~ NA_character_
    )
) 

dim(d5c) #232920, 12

# Read the data on PM2.5
PM_NFHS_5<- read.csv("Mean_PM_NFHS-5.csv") %>%
  rename(psu= DHSCLUST, PM2.5= MEAN)
dim(PM_NFHS_5) #27576, 2

d5c<- merge.data.frame(PM_NFHS_5, d5c, by="psu", all.y = TRUE)
dim(d5c) #232920, 13

# External Data Validation-Outcome of Anemia
match_id<- readRDS("NFHS5_NSS2223_match.rds")
nss_id<- readRDS("hces22_lvl_01_identification_v1.RDS")

hces2205_dat<- readRDS("hces2205_total_consumption_qty_in_100g_adj_v1_With_FG_Total.RDS")
df_1<- hces2205_dat[, c(1, seq(188, 200, 1))]
df_1[, names(df_1) != "common_id"] <- df_1[, names(df_1) != "common_id"] * 100
hces22_nutrient_v1 <- readRDS("hces22_nutrient_v1.RDS")

NSS = merge(df_1, hces22_nutrient_v1, by="common_id")
NSS = merge.data.frame(nss_id[,c("common_id", "hh_size")], NSS, by="common_id")
dim(NSS) #261696, 217

antioxidant_vars<- c("VA_RAE_mcg","TotalAscorbicAcid.VITC","VB12",
                      "ZincZn.ZN","SeleniumSe.SE",
                      "VITD","alphaTocopherolEquivalent.VITE","IronFe.FE")

# getting the per-capita nutrients
NSS[, c(antioxidant_vars)]=(NSS[,c(antioxidant_vars)]/NSS$hh_size)/30

d5c$hhid <- substr(trimws(d5c$caseid),
                       1, nchar(trimws(d5c$caseid)) - 2)
d5c$hhid<- trimws(d5c$hhid)
d5c_smoke$hhid<- trimws(d5c_smoke$hhid)

# merge the data for household smoking frequency 
d5c<- merge.data.frame(d5c, d5c_smoke, by= "hhid")
dim(d5c) #232920, 15

match_id$NFHS5_hhid<- trimws(match_id$NFHS5_hhid)

nfhs5_child_all= merge.data.frame(match_id, d5c, by.x= "NFHS5_hhid", by.y= "hhid")
# dim(nfhs_child_all)  #232920, 16

NFHS_NSS= merge.data.frame(nfhs5_child_all, NSS[,c("common_id", antioxidant_vars)],
                           by.x= "NSS22_23_hhuid", by.y= "common_id")

# rename the nutrient columns 
NFHS_NSS<- NFHS_NSS %>%
  rename(vitaminA= VA_RAE_mcg,
         vitaminB12= VB12,
         vitaminC= TotalAscorbicAcid.VITC,
         vitaminD= VITD,
         vitaminE= alphaTocopherolEquivalent.VITE,
         selenium= SeleniumSe.SE,
         zinc= ZincZn.ZN,
         iron= IronFe.FE)

# collecting the nutrient variables 
nutrients <- c("vitaminA", "vitaminB12", "vitaminC", "vitaminD", "vitaminE",
          "selenium", "zinc")

# collecting the descriptive variables
cat_vars <- c("child_gender", "wealthindex", "residence_type", "cooking_fuel", 
              "women_education", "anemia_levelc", "hh_frequency_smoke")

# collecting the numerical variables 
num_vars <- c(nutrients, "PM2.5", "iron")

# sampling variables
svars <- c("psu", "state", "strata", "sample_weight")

# to check if there are any missing data
anyNA(NFHS_NSS) #TRUE
dim(NFHS_NSS)  # 232797, 24

# subset the data so that children under age 5 are selected
data1 <- NFHS_NSS %>%
  filter(child_age >= 6)
dim(data1) #198363, 24
#View(data1)

# percentile trimming: we remove the outliers (below 0.5% and above 99.5%)
keep <- rep(TRUE, nrow(data1))

for (v in num_vars) {
  q <- quantile(data1[[v]],
                probs = c(0.005, 0.995),
                na.rm = TRUE)
  
  keep <- keep &
    data1[[v]] >= q[1] &
    data1[[v]] <= q[2]
}

d5c <- data1[keep, ]
dim(d5c) #187194, 24 

# counts of missingness -------------------------
all_vars<- c(num_vars, cat_vars) 

missing_summary<- data.frame(Variable = all_vars,
                  Missing = sapply(d5c[all_vars], \(x) sum(is.na(x))),
                 Non_missing = sapply(d5c[all_vars], \(x) sum(!is.na(x)))
)

missing_summary$Missing_pct<- round(100 * missing_summary$Missing 
                                    / nrow(d5c), 2); missing_summary

# select the data columns for further analysis
d5c<- d5c %>% select(., all_of(c(num_vars, cat_vars, svars)))
d5c<- na.omit(d5c) # remove the missing observations

dim(d5c)  #154,683; 20

# scale the nutrients for further analysis
d5c <- d5c %>%
  mutate(across(all_of(nutrients), ~ as.numeric(scale(.x)),
                .names = "scaled_{.col}"))


# create the DMI
wt<- read.csv("DMI_ChU5.csv")[,-1]
wt <- ifelse(is.na(wt), 0, wt)/100

# Compute DMI for children
d5c$DMI<- as.numeric(as.matrix(d5c[, paste("scaled", 
                                               nutrients, sep ="_")]) %*% wt)
d5c$DMIscaled<- (d5c$DMI- (-3))/ (6-(-3))

# store the NFHS-5 data for children for further use
write.csv(d5c, "NFHS5_cleandata_children.csv")

#Descriptive statistics for categorical variables
tab2 <- lapply(cat_vars, function(v) {
  tbl <- table(d5c[[v]])
  pct <- prop.table(tbl)* 100
  data.frame(
    Frequency  = as.numeric(tbl),
    Percentage = round(as.numeric(pct), 1),
    row.names  = names(tbl)
  )
})

names(tab2) <- cat_vars
catvar_labels <- c(child_gender= "Sex of Child", 
                   wealthindex= "Wealth Index",
                   residence_type= "Type of Locality", 
                   cooking_fuel= "Cooking Fuel",
                   hh_frequency_smoke= "Smoking frequency per household",
                   women_education= "Mothers Education",
                   sh_smoke= "Second-hand smoke",
                   anemia_levelc= "Anemia Status")

tab2_combined <- do.call(rbind, lapply(names(tab2), function(var) {
    df <- tab2[[var]]
    data.frame(Variable = c(catvar_labels[var], rep("", max(nrow(df) - 1, 0))),
      Category = rownames(df), Frequency = df$Frequency, 
      Percentage = df$Percentage)}))

tab2_combined$`N(%)`<- paste(tab2_combined$Frequency, " (", 
                             tab2_combined$Percentage, ")", sep ="")

tab2_combined
write.csv(tab2_combined, "Table2_nfhs5_ChU5.csv")

# to check whether the psu's are nested within the states in the data
all(tapply(d5c$strata,
            d5c$psu,
         function(x) length(unique(x))) == 1)

# summary measures for the nutrients
me<- apply(d5c[, nutrients], 2, median)
iqr<- apply(d5c[, nutrients], 2, IQR)
mn<- apply(d5c[, nutrients], 2, mean)
std<- apply(d5c[, nutrients], 2, sd)

as.data.frame(bind_rows(me, iqr, mn, std))

#----------------------------------
# run the survey glm regression 
#----------------------------------

options(survey.lonely.psu = "adjust")
S.design <- svydesign(id =~ psu, strata= ~strata, weights = ~sample_weight,
                      data = d5c, nest = TRUE)

S.design$variables$women_education <- relevel(factor(S.design$variables$women_education),
  ref = "No education")

S.design$variables$wealthindex <- relevel(factor(S.design$variables$wealthindex),
  ref = "Poorest")

S.design$variables$residence_type<- relevel(factor(S.design$variables$residence_type),
                                          ref = "Urban")

S.design$variables$hh_frequency_smoke<- relevel(factor(S.design$variables$hh_frequency_smoke),
                                            ref= "never")

# fit the model for children
fit_c<- svyglm(anemia_levelc ~ child_gender + hh_frequency_smoke + residence_type+ 
                 cooking_fuel+ wealthindex + iron + women_education + PM2.5+ 
                 (PM2.5:DMIscaled),
              design = S.design,
              family = quasibinomial())

summary(fit_c)

#converting DMI_anemia into quartiles
d5c$DMI_q <- cut(d5c$DMIscaled, breaks = quantile(d5c$DMIscaled,
                    probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
  include.lowest = TRUE, labels = c("Min - Q1", "Q1 - Q2", "Q2 - Q3", "Q3 - Max")
)

# create the design object again
S.design <- svydesign(id= ~psu, strata= ~state, weights= ~woman_sampling_weight, 
                      data= d5c, nest = TRUE)

# run the regression again with the quantile variable for the DMI scores
fit4<- svyglm(anemia_level ~ child_gender + residence_type + cooking_fuel + 
                wealthindex + iron+ women_education + hh_frequency_smoke + 
                PM2.5+ PM2.5:DMI_q, design = S.design, family = quasibinomial())

summary(fit4)

#for binary variables mean=prevalence
svyby(formula=~ anemia_levelc, by=~ child_gender, design= S.design, 
      FUN= svymean, vartype= "ci")
