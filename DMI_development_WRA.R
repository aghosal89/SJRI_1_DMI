rm(list=ls())
### Required Libraries ####
library(dplyr)
library(magrittr)
library(data.table)
library(readr)
#### Setting working directory ####
setwd("/Volumes/Dev/PRJ_BIRACKnit2_WRA")
source("~/Library/CloudStorage/OneDrive-StJohn'sNationalAcademyofHealthSciences/Consulatation/functions.r")
## Getting Matched ID with NSSO ###
load("/Volumes/Dev/PRJ_BIRACKnit2_WRA/Data/Alldata.RData")
Matched_ID<-Alldata$NFHS4
Matched_ID<-Matched_ID[,c("donHH.id","ID","hh_size_L04",names(Matched_ID)[51:535])]
names(Matched_ID)[1:2]<-c("NSSO_HHID","NFHS_ID")
### WRA Data from NFHS with Confounders ######
nfhs4_wm <- readRDS("/Volumes/Dev/PRJ_BIRACKnit2_WRA/Data/nfhs4_wm.rds")
NSSO_nutrients <- readRDS("/Volumes/Dev/PRJ_BIRACKnit2_WRA/Data/NSSO_nutrients.rds")
NSSO_nutrients$NFHS_ID<-as.character(NSSO_nutrients$NFHS_ID)
Matched_ID=Matched_ID[,names(Matched_ID)[1:2]]
nfhs4_wm1<-merge(nfhs4_wm,Matched_ID,by.x = "caseid",by.y = "NFHS_ID",all.x = T)
nfhs4_wm1<-nfhs4_wm1[!is.na(nfhs4_wm1$NSSO_HHID),]
nfhs4_wm1<-as.data.frame(nfhs4_wm1)
nfhs4_wm1<-merge(nfhs4_wm1,NSSO_nutrients,by="NSSO_HHID",all.x = T)
nfhs4_wm1<-nfhs4_wm1[!is.na(nfhs4_wm1$pc_alphaTocopherolEquivalent.VITE),]
## Percapita Nutrient data for HH

percapita_nfhs4<-nfhs4_wm1[!duplicated(nfhs4_wm1$hhid),]
percapita_nfhs4<-percapita_nfhs4[,c("hhid","pc_TotalAscorbicAcid.VITC" #TotalAscorbicAcid.VITC
                                    ,"pc_SeleniumSe.SE" #SeleniumSe.SE
                                    ,"pc_ZincZn.ZN" # Zinc ZN
                                    ,"pc_VitaminA.RAE" # VitaminA.RAE
                                    ,"pc_VitaminD.D2.D3" # VitaminD.D2.D3
                                    ,"pc_IronFe.FE" # IronFe.FE
                                    ,"pc_VitaminB12" # Vitamin B12
                                    ,"pc_alphaTocopherolEquivalent.VITE" )]


saveRDS(percapita_nfhs4,"percapita_nfhs4.rds")
### Exposure Data ####
nfhs4_h_ide<- readRDS("~/Library/CloudStorage/OneDrive-SharedLibraries-StJohn'sNationalAcademyofHealthSciences/TATA NIN - nfhs/nfhs4/working/rds/nfhs4_h_identification.rds")
nfhs4_h_ide<-nfhs4_h_ide[nfhs4_h_ide$hhid %in% nfhs4_wm1$hhid,]
nfhs4_wm_exposure<-merge(nfhs4_wm1[,c("caseid","hhid","month_interview","year_interview","date_interview")],nfhs4_h_ide[,c("hhid","hv001")],by="hhid",all.x = TRUE)
nfhs4_wm_exposure$hv001_=if_else(nchar(nfhs4_wm_exposure$hv001)==5,
                                 paste0("0",nfhs4_wm_exposure$hv001),
                                 as.character(nfhs4_wm_exposure$hv001))
nfhs4_wm_exposure$DHSID=paste0("IA201400",nfhs4_wm_exposure$hv001_)


data_NFHS4GPS <- read.csv("Data/data_NFHS4GPS.csv")

nfhs4_wm_exposure<-merge(nfhs4_wm_exposure,data_NFHS4GPS,by = "DHSID",all.x = TRUE)


pm_years <- 2007:2016

calc_pm25_5yr_avg <- function(interview_year, interview_month, pm_values) {
  
  end_year <- interview_year
  end_month <- interview_month
  
  start_year <- end_year
  start_month <- end_month - 59
  
  while(start_month <= 0) {
    start_month <- start_month + 12
    start_year <- start_year - 1
  }
  
  total <- 0
  total_months <- 0
  
  for (yr in pm_years) {
    
    if (yr < start_year | yr > end_year) next
    
    if (yr == start_year & yr == end_year) {
      months <- end_month - start_month + 1
    } else if (yr == start_year) {
      months <- 12 - start_month + 1
    } else if (yr == end_year) {
      months <- end_month
    } else {
      months <- 12
    }
    
    pm <- pm_values[as.character(yr)]
    
    if (!is.na(pm)) {
      total <- total + pm * months
      total_months <- total_months + months
    }
  }
  
  if (total_months == 0) return(NA)
  total / total_months
}

library(dplyr)

pm_cols <- paste0("X", 2007:2016)

Data <- nfhs4_wm_exposure %>%
  rowwise() %>%
  mutate(
    pm25_5yr_avg = calc_pm25_5yr_avg(
      year_interview,
      month_interview,
      c_across(all_of(pm_cols)) %>% setNames(2007:2016)
    )
  ) %>%
  ungroup()


## Selecting required columns ####

Data<-Data[,c("caseid","pm25_5yr_avg")]
selected_vars<-c("caseid","hhid","age","psu","state","residence","education","wealth_index","sample_weight"
                 ,"BMI_cat","second_hand_smoke","cooking_fuel_group","Hb_level"
                 ,"pc_TotalAscorbicAcid.VITC" #TotalAscorbicAcid.VITC
                 ,"pc_SeleniumSe.SE" #SeleniumSe.SE
                 ,"pc_ZincZn.ZN" # Zinc ZN
                 ,"pc_VitaminA.RAE" # VitaminA.RAE
                 ,"pc_VitaminD.D2.D3" # VitaminD.D2.D3
                 ,"pc_IronFe.FE" # IronFe.FE
                 ,"pc_VitaminB12" # Vitamin B12
                 ,"pc_alphaTocopherolEquivalent.VITE" # alphaTocopherolEquivalent.VITE
)
nfhs4_wm_selected_vars<-nfhs4_wm1[,selected_vars]
### adding exposurer data ###
nfhs4_wm_model_data<-merge(nfhs4_wm_selected_vars,Data,by = "caseid",all.x = TRUE)

sapply(nfhs4_wm_model_data, function(x) sum(is.na(x)))

nfhs4_wm_model_data<-nfhs4_wm_model_data[!is.na(nfhs4_wm_model_data$pm25_5yr_avg),]

### Bootsrapping & Ridge Regression ################

d4w<- nfhs4_wm_model_data %>%
  rename(vitaminA= pc_VitaminA.RAE, vitaminB12= pc_VitaminB12, 
         vitaminC= pc_TotalAscorbicAcid.VITC, vitaminD= pc_VitaminD.D2.D3,
         vitaminE= pc_alphaTocopherolEquivalent.VITE, selenium= pc_SeleniumSe.SE,
         zinc= pc_ZincZn.ZN, iron= pc_IronFe.FE, pm2.5= pm25_5yr_avg) %>% 
  mutate(across(all_of(nutrients), ~as.numeric(scale(.x)), .names = "scaled_{.col}"))

relWeight <- function(data, nsim, NRS){
  
  interaction_terms<- c(
    "pm2.5:scaled_vitaminA",
    "pm2.5:scaled_vitaminB12",
    "pm2.5:scaled_vitaminC",
    "pm2.5:scaled_vitaminD",
    "pm2.5:scaled_vitaminE",
    "pm2.5:scaled_selenium",
    "pm2.5:scaled_zinc"
  )
  
  # Pre-allocate matrix
  est <- matrix(NA, nrow = nsim, ncol = length(interaction_terms))
  
  colnames(est)<- interaction_terms
  
  # Simulation indices for bootstrap
  sim_data<- replicate(nsim, sample(seq_len(nrow(data)), NRS, replace= FALSE))
  
  for(i in 1:nsim){
    
    temp <- sim_data[, i]
    mydata <- data[temp, ]
    
    # Scale nutrients
    mydata <- mydata %>% mutate(across(all_of(nutrients), ~as.numeric(scale(.x)),
                                       .names = "scaled_{.col}"))
    
    # Design matrix
    X <- model.matrix(Hb_level ~ pm2.5* (scaled_vitaminA + scaled_vitaminB12 +
                                           scaled_vitaminC + scaled_vitaminD + scaled_vitaminE + scaled_selenium +
                                           scaled_zinc ) + residence + education + wealth_index + BMI_cat + 
                        second_hand_smoke + cooking_fuel_group + iron, data = mydata)[,-1]
    
    Y <- mydata$Hb_level
    
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

est <- relWeight(data= d4w, nsim= 200, NRS= 20000)

# Apply the scores on the estimates of the regression
tab <- do.call(rbind, apply(est, 2, score))
tab <- as.data.frame(tab)
rownames(tab) <- c("Vitamin A", "Vitamin B", "Vitamin C", "Vitamin D", "Vitamin E",
                   "Se", "Zn")

# Positive proportion
tab$positive_prop<- tab$p.n/ nsim

# Weight calculation in percentage form
wt <-  100*(tab$p.mean* tab$positive_prop)/ sum(tab$p.mean* tab$positive_prop,
                                                na.rm = TRUE)
tab$weight<- wt
names(wt)<- nutrients

# Export the table
#write_xlsx(as.data.frame(tab), "Tablehb_wts.xlsx")

# export the relative weights for the next stage of the analysis
write.csv(wt, "DMI_wra.csv")

wt <- ifelse(is.na(wt), 0, wt)/100

# Compute DMI for children
d4w$DMI<- as.numeric(as.matrix(d4w[, paste("scaled", 
                                           names(wt), sep ="_")]) %*% wt)

d4w$DMIscaled<- (d4w$DMI- (-3))/ (6-(-3))



