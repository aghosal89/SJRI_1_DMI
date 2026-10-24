rm(list = ls())
## Required Libraries ###
library(dplyr)
library(haven)
library(data.table)
library(sjlabelled)
library(survey)
library(readxl)
setwd("/Volumes/Dev/PRJ_BIRACKnit2_WRA")
### Reading Data ####
nfhs5_wm = read_dta("~/Library/CloudStorage/OneDrive-SharedLibraries-StJohn'sNationalAcademyofHealthSciences/TATA NIN - nfhs/nfhs5/raw_data/IAIR7ADT//IAIR7AFL.DTA",
                    col_select = c("caseid"   #women id
                                   ,"v454"   #ipregnant
                                   ,"v005"   #sample Weight
                                   ,"v021"   #psu
                                   , "v024" #state
                                   , "v025" #urban/rural
                                   ,"v106" #educational attainment
                                   , "s116" #caste
                                   , "v190" #wealth index
                                   ,"v006"# end month
                                   ,"v007" # end year
                                   ,"v009"# month_birth
                                   ,"v010" # month_year
                                   ,"v456"# Hb level
                                   ,"v437" # weight
                                   ,"v438" # height
                                   #,"v463a" #smoking 
                                   #,"v463b" # smoke pipe usage
                                   #,"v463c" # tobacco 
                                   #,"v463d" # snuff usage
                                   #,"v463e" # smokes cigars
                                   #,"v463f" # paan masala
                                   # ,"v463g" # paan masala with tobacco
                                   # ,"v463x" # smokes other
                                   ,"v012" # Age
                                   ,"v404" # Breastfeeding
                                   ,"s719" # someone else smoke in hh
                                   ,"v161" ## cooking fuel
                                   ,"v022" # strata
                    ))


names(nfhs5_wm)<-c("caseid","sample_weight","end_month","end_year","month_birth","month_year","Age","psu","strata","state","residence","education","cooking_fuel",
                   "wealth_index","breastfeeding","weight","height","pregnant","Hb","caste","second_hand_smoke")


#---------- Data Cleaning ------------#
nfhs5_wm<-nfhs5_wm[nfhs5_wm$pregnant==0 & nfhs5_wm$breastfeeding==0,] # Non-pregnant & lactating mothers

nfhs5_wm<-nfhs5_wm[nfhs5_wm$Hb>=40 & nfhs5_wm$Hb<180,] # plaussible Hb
nfhs5_wm$Hb<-nfhs5_wm$Hb/10
nfhs5_wm<-nfhs5_wm[!nfhs5_wm$weight%in%c(9994,9995,9996),]
nfhs5_wm<-nfhs5_wm[!nfhs5_wm$height%in%c(9994,9995,9996),]
nfhs5_wm$BMI<-(nfhs5_wm$weight/10)/(nfhs5_wm$height/1000)^2
nfhs5_wm$BMI %>% range(na.rm = TRUE)
nfhs5_wm<-nfhs5_wm[nfhs5_wm$BMI>=12 & nfhs5_wm$BMI<=70,] # BMI outliers

nfhs5_wm$residence<-factor(nfhs5_wm$residence,levels = c(1,2),labels = c("urban","rural"))
nfhs5_wm$education<-factor(nfhs5_wm$education,levels =c(0,1,2,3),labels = c("no education","primary","secondary","higher"))
nfhs5_wm$wealth_index<-factor(nfhs5_wm$wealth_index,levels =c(1,2,3,4,5),labels = c("poorest","poor","middle","richer","richest"))
nfhs5_wm$caste<-ifelse(nfhs5_wm$caste %in% c(1,2),"SC/ST",
                       ifelse(nfhs5_wm$caste == 3, "OBC","Others"))


nfhs5_wm$second_hand_smoke<-factor(nfhs5_wm$second_hand_smoke,levels = c(0,1),labels = c("No","Yes"))
nfhs5_wm$cooking_fuel_group <- dplyr::case_when(
  nfhs5_wm$cooking_fuel %in% c(1,2,3,4) ~ "Clean fuel",
  nfhs5_wm$cooking_fuel == 5 ~ "Kerosene",
  nfhs5_wm$cooking_fuel %in% c(6,7,8,9,10,11) ~ "Solid fuel",
  TRUE ~ NA_character_
)

nfhs5_wm$BMI_cat<-ifelse(nfhs5_wm$BMI<18.5,"underweight",
                         ifelse((nfhs5_wm$BMI>=18.5 & nfhs5_wm$BMI<25),"Normalweight",
                                ifelse((nfhs5_wm$BMI>=25 & nfhs5_wm$BMI<30),"Overweight","Obese")))

saveRDS(nfhs5_wm,"Data/nfhs5_wm_covariates.rds")

##---------- Exposure data -------------##
## We can't use year wise exposure data for NFHS4, that francio will send new one, 5 year avg exposurer data
# Cluster Level PM2.5 -----------------------------------------------------
Mean_PM_NFHS_5 <- read.csv("Data/Mean_PM_NFHS-5.csv")
Mean_PM_NFHS_5=dplyr::rename(Mean_PM_NFHS_5,clust_no=DHSCLUST,PM2.5=MEAN)

NFHS_PM=merge.data.frame(Mean_PM_NFHS_5,nfhs5_wm,by.x="clust_no",by.y = "psu",all.y = TRUE)
NFHS_PM<-NFHS_PM[!is.na(NFHS_PM$PM2.5),]
saveRDS(NFHS_PM,"Data/NFHS_PM.RDS")

##---------- Nutrient data -------------##

# Triangulated ID ---------------------------------------------------------
matchID= readRDS("/Volumes/Dev/PRJ_BIRACKnit2_WRA/Codes/francio/NFHS5_NSS2223_match.rds")
# NSS 22-23 ---------------------------------------------------------------
NSS_path="~/Library/CloudStorage/OneDrive-SharedLibraries-StJohn'sNationalAcademyofHealthSciences/Mr.Jawahar.RM - NSSO 2022-2023/working/rds_files"
NSS_identification=readRDS(file = paste0(NSS_path,"/hces22_lvl_01_identification_v1.RDS"))


hces2205_total_consumption_qty_in_100g_adj_v1_With_FG_Total=readRDS(file = paste0(NSS_path,"/hces2205_total_consumption_qty_in_100g_adj_v1_With_FG_Total.RDS"))
df_1 = hces2205_total_consumption_qty_in_100g_adj_v1_With_FG_Total[,c(1,seq(188,200,1))]
df_1[, names(df_1) != "common_id"] <- df_1[, names(df_1) != "common_id"] * 100
hces22_nutrient_v1=readRDS(file = paste0(NSS_path,"/hces22_nutrient_v1.RDS"))

NSS = merge(df_1,hces22_nutrient_v1,by="common_id")
NSS = merge.data.frame(NSS_identification[,c("common_id","hh_size")],NSS,by="common_id")

antioxidant_vars <- c("VA_RAE_mcg","TotalAscorbicAcid.VITC","VB12",
                      "ZincZn.ZN","SeleniumSe.SE",
                      "VITD","alphaTocopherolEquivalent.VITE","IronFe.FE")
macro_vars <- c("Energy.ENERC.Kcal","Protein.PROTCNT",
                "Carbohydrate.CHOAVLDF","TotalFat.FATCE",
                "TotalDietaryFibre.FIBTG")
vars_for_IIT=c("Fresh_Fruits_T","Vegetables_T")
NSS[,c(antioxidant_vars,macro_vars,vars_for_IIT)]=(NSS[,c(antioxidant_vars,macro_vars,vars_for_IIT)]/NSS$hh_size)/30
summary(NSS[,c(antioxidant_vars,macro_vars,vars_for_IIT)])

trim01 <- function(x, p=0.005){q <- quantile(x, c(p,1-p), na.rm=TRUE); x[x<q[1] | x>q[2]] <- NA; x}
NSS[, c(antioxidant_vars, macro_vars,vars_for_IIT)] <- lapply(NSS[, c(antioxidant_vars, macro_vars,vars_for_IIT)], trim01)
summary(NSS[,c(antioxidant_vars,macro_vars,vars_for_IIT)])

# Merging NFHS and NSS ----------------------------------------------------
NFHS_PM <- readRDS("/Volumes/Dev/PRJ_BIRACKnit2_WRA/Data/NFHS_PM.RDS")
NFHS_PM$hhid <- substr(trimws(NFHS_PM$caseid),
                       1, nchar(trimws(NFHS_PM$caseid)) - 2)
NFHS_PM$hhid=trimws(NFHS_PM$hhid)
#NFHS_PM <- NFHS_PM %>% select(hhid, everything())
matchID$NFHS5_hhid=trimws(matchID$NFHS5_hhid)
nfhs5_woman_all=merge.data.frame(matchID,NFHS_PM,by.x = "NFHS5_hhid",by.y = "hhid")
NFHS_NSS=merge.data.frame(nfhs5_woman_all,NSS[,c("common_id",antioxidant_vars,macro_vars,vars_for_IIT)],
                          by.x="NSS22_23_hhuid",by.y="common_id")

summary(NFHS_NSS[,c(antioxidant_vars,macro_vars,vars_for_IIT)])
saveRDS(NFHS_NSS,"Data/NFHS_NSS_IIT_Delhi.RDS")
### Percapita nutrient for household 
percapita_nfhs5<-NFHS_NSS[!duplicated(NFHS_NSS$NFHS5_hhid),]
percapita_nfhs5<-percapita_nfhs5[,c("NFHS5_hhid",antioxidant_vars,macro_vars)]

## DMI creation ----------------------------------------------------------
NFHS_NSS<-readRDS("Data/NFHS_NSS_IIT_Delhi.RDS")
NFHS_NSS$PM2.5_IQR<-NFHS_NSS$PM2.5/IQR(NFHS_NSS$PM2.5)

nfhs5_nutrients<-NFHS_NSS[,c( "TotalAscorbicAcid.VITC", "SeleniumSe.SE", "ZincZn.ZN", "VA_RAE_mcg"
                              , "VITD","VB12", "alphaTocopherolEquivalent.VITE")]

names(nfhs5_nutrients)<-c("VitC","Selenium","Zinc","VitA","VitD","VitB12","VitE")
nfhs5_nutrients$survey<-rep("NFHS5",nrow(nfhs5_nutrients))
nfhs4_nutients <- readRDS("/Volumes/Dev/PRJ_BIRACKnit2_WRA/nfhs4_nutients.rds")
nutrient_data<-rbind(nfhs4_nutients,nfhs5_nutrients)
#saveRDS(nutrient_data,"nutrient_data.rds")
## Scaling of antioxidants
NFHS_NSS %>% summary()
names(NFHS_NSS)
NFHS_NSS <-NFHS_NSS[,1:35]
NFHS_NSS<-na.omit(NFHS_NSS)
NFHS_NSS[antioxidant_vars[1:7]] <- lapply(NFHS_NSS[antioxidant_vars[1:7]], function(x) {
  (x - mean(x)) / (sd(x))
})

sapply(NFHS_NSS[antioxidant_vars], summary)

## Compute DMI
## -3 to 6 after DMI calculate, then apply it normal
DMI_wra <- read.csv("Data/DMI_wra.csv")
DMI_wra$X<-c("VA_RAE_mcg","VB12","TotalAscorbicAcid.VITC","VITD","alphaTocopherolEquivalent.VITE","SeleniumSe.SE","ZincZn.ZN")
DMI_wra<-DMI_wra[match(antioxidant_vars[1:7],DMI_wra$X),]

DMI<-as.matrix(NFHS_NSS[,antioxidant_vars[1:7]])%*% DMI_wra[,2]

NFHS_NSS$DMI<-DMI[,1]
NFHS_NSS$DMInorm<-(NFHS_NSS$DMI-(-3))/(6-(-3))
## Model creation ----------------------------------------------------------
NFHS_NSS$anemia_level<-ifelse(NFHS_NSS$Hb<12,1,0)
data_nfhs5_for_plot<-NFHS_NSS[,c("caseid","DMInorm","anemia_level","sample_weight","state")]
psu_data<-readRDS("Data/nfhs5_wm_covariates.rds")
psu_data<-psu_data[,c("caseid","psu")]
data_nfhs5_for_plot<-merge(data_nfhs5_for_plot,psu_data,by = "caseid",all.x = T) 
saveRDS(data_nfhs5_for_plot,"data_nfhs5_for_plot.rds")
#NFHS_NSS<-NFHS_NSS[!NFHS_NSS$strata==171,]
NFHS_NSS$sample_weight1<-NFHS_NSS$sample_weight/1000000
options(survey.lonely.psu = "average")
library(dplyr)
NFHS_NSS$state<-as.factor(NFHS_NSS$state)
saveRDS(NFHS_NSS,"Data/NFHS5_validation.rds")
S.design <- svydesign(id= ~clust_no, weights= ~sample_weight1, 
                      data= NFHS_NSS,strata =~ strata)


fit<- svyglm(anemia_level ~ education   + wealth_index+
               BMI_cat + residence+second_hand_smoke+cooking_fuel_group+PM2.5+PM2.5:DMInorm+IronFe.FE,
             design = S.design, 
             family = quasibinomial())


summary(fit)



result <- cbind(
  OR = exp(coef(fit)),
  exp(confint(fit)),
  p_value = summary(fit)$coefficients[,4]
)
result

vcm<- survey:::vcov.svyglm(fit) [c("PM2.5", "PM2.5:DMI"), 
                                 c("PM2.5", "PM2.5:DMI")]
cm<- coef(fit)[c("PM2.5", "PM2.5:DMI")]

conf_band <- function(x, coef_m, vcov_m) {
  
  est <- coef_m[1] + coef_m[2] * x
  se <- sqrt(matrix(c(1, x), nrow = 1) %*% vcov_m %*% matrix(c(1, x), ncol = 1))
  
  lb <- exp(est - 1.96 * se)
  ub <- exp(est + 1.96 * se)
  
  return(c("estimate" = exp(est), "lower" = lb, "upper" = ub))
}

cf <- t(sapply(NFHS_NSS$DMI, FUN = conf_band, coef_m = cm,
               vcov_m = vcm))
colnames(cf)<- c("Odds", "Lower", "Upper")
plot_df<- cbind.data.frame(DMI = NFHS_NSS$DMI, cf)


ggplot(plot_df, aes(x = DMI, y = Odds)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper),
              alpha = 0.2,
              col = "red") +
  theme_classic() +
  geom_boxplot(
    data = plot_df,
    aes(y = 0.98, x = DMI),
    outlier.size = 0.3,
    width = 0.005
  ) +
  labs(x = "DMI", y = TeX(r'($\beta$ coefficient for PM$_{2.5}$)'))


