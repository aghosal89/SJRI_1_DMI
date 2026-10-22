
#=====================================================================
# This script is used to create the figures and tables in the paper
#=====================================================================

# clear the memory first
rm(list = ls())
gc()

# Set the working directory
setwd("~/Documents/Aritra_sjri_files")

# read the coding libraries necessary for our analysis
library(dplyr)
library(ggplot2)
library(ggh4x)
library(gridExtra)
library(patchwork)
library(ggeffects)
library(corrplot)
library(car)
library(latex2exp)
library(sf)
library(geodata)
library(ggrepel)
library(viridis)
library(stringr)
library(tidyr)
library(survey)

#===================================================
# data processing for women and children in NFHS-5
#===================================================

# read the data for children for both the surveys 
d4c<- read.csv("NFHS4_cleandata_children.csv")
d5c<- read.csv("NFHS5_cleandata_children.csv")

# collecting the descriptive variables
cat_vars<- c("child_gender", "wealthindex", "residence_type", "cooking_fuel", 
             "women_education", "anemia_levelc", "state", "hh_frequency_smoke",
             "psu")

# collecting the nutrient variables
nutrients<- c("vitaminA", "vitaminB12", "vitaminC", "vitaminD", "vitaminE", 
              "selenium", "zinc")

# collecting the numerical variables 
num_vars <- c(nutrients, "PM2.5", "iron", "child_age", "sample_weight")

# read the data for women for both the surveys 
d4w<- readRDS("Aswathi_nfhs4_wm_model_data.rds") %>%
  rename(vitaminA= pc_VitaminA.RAE, vitaminB12= pc_VitaminB12, 
         vitaminC= pc_TotalAscorbicAcid.VITC, vitaminD= pc_VitaminD.D2.D3,
         vitaminE= pc_alphaTocopherolEquivalent.VITE, selenium= pc_SeleniumSe.SE,
         zinc= pc_ZincZn.ZN, iron= pc_IronFe.FE, pm2.5= pm25_5yr_avg) %>% 
  mutate(across(all_of(nutrients), ~as.numeric(scale(.x)), 
                .names = "scaled_{.col}")) %>%
  mutate(state = case_when(
    state == 1  ~ "Jammu & Kashmir",
    state == 2  ~ "Himachal Pradesh",
    state == 3  ~ "Punjab",
    state == 4  ~ "Chandigarh",
    state == 5  ~ "Uttarakhand",
    state == 6  ~ "Haryana",
    state == 7  ~ "NCT of Delhi",
    state == 8  ~ "Rajasthan",
    state == 9  ~ "Uttar Pradesh",
    state == 10 ~ "Bihar",
    state == 11 ~ "Sikkim",
    state == 12 ~ "Arunachal Pradesh",
    state == 13 ~ "Nagaland",
    state == 14 ~ "Manipur",
    state == 15 ~ "Mizoram",
    state == 16 ~ "Tripura",
    state == 17 ~ "Meghalaya",
    state == 18 ~ "Assam",
    state == 19 ~ "West Bengal",
    state == 20 ~ "Jharkhand",
    state == 21 ~ "Odisha",
    state == 22 ~ "Chhattisgarh",
    state == 23 ~ "Madhya Pradesh",
    state == 24 ~ "Gujarat",
    state == 25 ~ "Dadra & Nagar Haveli and Daman & Diu",
    state == 27 ~ "Maharashtra",
    state == 28 ~ "Andhra Pradesh",
    state == 29 ~ "Karnataka",
    state == 30 ~ "Goa",
    state == 31 ~ "Lakshadweep",
    state == 32 ~ "Kerala",
    state == 33 ~ "Tamil Nadu",
    state == 34 ~ "Puducherry",
    state == 35 ~ "Andaman & Nicobar Islands",
    state == 36 ~ "Telangana",
    state == 37 ~ "Ladakh",
    TRUE ~ NA_character_
  )) 

d4w$state<-as.factor(d4w$state)

d5w<- readRDS("Aswathi_NFHS5_validation.rds") %>%
  rename(vitaminA= VA_RAE_mcg, vitaminB12= VB12, 
         vitaminC= TotalAscorbicAcid.VITC, vitaminD= VITD,
         vitaminE= alphaTocopherolEquivalent.VITE, selenium= SeleniumSe.SE,
         zinc= ZincZn.ZN, iron= IronFe.FE) %>%
  mutate(state = case_when(
    state == 1  ~ "Jammu & Kashmir",
    state == 2  ~ "Himachal Pradesh",
    state == 3  ~ "Punjab",
    state == 4  ~ "Chandigarh",
    state == 5  ~ "Uttarakhand",
    state == 6  ~ "Haryana",
    state == 7  ~ "NCT of Delhi",
    state == 8  ~ "Rajasthan",
    state == 9  ~ "Uttar Pradesh",
    state == 10 ~ "Bihar",
    state == 11 ~ "Sikkim",
    state == 12 ~ "Arunachal Pradesh",
    state == 13 ~ "Nagaland",
    state == 14 ~ "Manipur",
    state == 15 ~ "Mizoram",
    state == 16 ~ "Tripura",
    state == 17 ~ "Meghalaya",
    state == 18 ~ "Assam",
    state == 19 ~ "West Bengal",
    state == 20 ~ "Jharkhand",
    state == 21 ~ "Odisha",
    state == 22 ~ "Chhattisgarh",
    state == 23 ~ "Madhya Pradesh",
    state == 24 ~ "Gujarat",
    state == 25 ~ "Dadra & Nagar Haveli and Daman & Diu",
    state == 27 ~ "Maharashtra",
    state == 28 ~ "Andhra Pradesh",
    state == 29 ~ "Karnataka",
    state == 30 ~ "Goa",
    state == 31 ~ "Lakshadweep",
    state == 32 ~ "Kerala",
    state == 33 ~ "Tamil Nadu",
    state == 34 ~ "Puducherry",
    state == 35 ~ "Andaman & Nicobar Islands",
    state == 36 ~ "Telangana",
    state == 37 ~ "Ladakh",
    TRUE ~ NA_character_
  )) 

d5w$state<-as.factor(d5w$state)

#----------------------------------------------------------------
# Fit logistic regression on anaemia prevalence for WRA and ChU5
#----------------------------------------------------------------

# for women -------------------------

# construct the survey design
dmi_w<- read.csv("DMI_wra.csv")[,-c(1,2)]

dmi_w <- ifelse(is.na(dmi_w), 0, dmi_w)/100

# Compute DMI for children
d5w$DMI<- as.numeric(as.matrix(d5w[, nutrients]) %*% as.matrix(dmi_w) )

d5w$DMIscaled<- (d5w$DMI- (-3))/(6-(-3))

options(survey.lonely.psu= "average")
S.design_w<- svydesign(id= ~clust_no, weights= ~sample_weight1, 
                      data= d5w, strata=~ strata)

fit_w<- svyglm(anemia_level ~ education+ wealth_index+ iron+
               BMI_cat+ residence+second_hand_smoke+ cooking_fuel_group+ PM2.5+
               PM2.5:DMIscaled, design= S.design_w, family = quasibinomial())

summary(fit_w)

# for children --------------------------

# define the survey structure
options(survey.lonely.psu= "average")
S.design_c <- svydesign(id= ~psu, strata= ~strata, weights= ~sample_weight, 
                        data= d5c, nest = TRUE)

# fit the model for children
S.design_c$variables$women_education <- relevel(factor(S.design_c$variables$women_education),
                                              ref = "No education")

S.design_c$variables$wealthindex <- relevel(factor(S.design_c$variables$wealthindex),
                                          ref = "Poorest")

S.design_c$variables$residence_type<- relevel(factor(S.design_c$variables$residence_type),
                                            ref = "Urban")

S.design_c$variables$hh_frequency_smoke<- relevel(factor(S.design_c$variables$hh_frequency_smoke),
                                              ref = "never")

# fitted model for the children
fit_c<- svyglm(anemia_levelc ~ child_gender + cooking_fuel + hh_frequency_smoke + 
                 residence_type + cooking_fuel + wealthindex + iron + 
                 women_education + PM2.5+ (PM2.5:DMIscaled),
               design = S.design_c, 
               family = quasibinomial())

summary(fit_c)

# create state-wide summary for children (ChU5)
state_data5c<- svyby(~PM2.5 + anemia_levelc + DMIscaled, ~state, design = S.design_c,
                      FUN = svymean, na.rm = TRUE) %>%
  rename(State= state, mean_PM2.5= PM2.5, Prevalence= anemia_levelc, 
         Mean_DMI = DMIscaled) %>%
  mutate(Category= "ChU5 (6-59 months)", Mean_DMI= Mean_DMI*100)

# create state-wide summary data for women (WRA)
state_data5w<- svyby(~PM2.5 + anemia_level + DMIscaled, ~state, design= S.design_w,
                     FUN = svymean, na.rm = TRUE) %>%
  rename(State= state, mean_PM2.5= PM2.5, Prevalence = anemia_level,
         Mean_DMI = DMIscaled) %>%
  mutate(Category= "WRA (15-49 years)", Mean_DMI= Mean_DMI*100)

#state_data5w<- d5w %>%
#  group_by(state) %>%
#  summarise(
#    Prevalence = 100 * mean(anemia_level, na.rm = TRUE),
#    Mean_DMI = mean(DMI, na.rm = TRUE),
#    mean_PM2.5= mean(PM2.5, na.rm= TRUE),
#    n = n(),
#    .groups = "drop"
#  ) %>% rename(State= state)

#state_data5c<- d5c %>%
#  group_by(state) %>%
#  summarise(
#    Prevalence = 100 * mean(anemia_levelc, na.rm = TRUE),
#    Mean_DMI = mean(DMI, na.rm = TRUE),
#    mean_PM2.5= mean(PM2.5, na.rm= TRUE),
#    n = n(),
#    .groups = "drop"
#  ) %>% rename(State= state)

# storing the state names as characters
state_data5c$State <- as.character(state_data5c$State)
state_data5w$State <- as.character(state_data5w$State)

# read the shape files for india's map
india_sf <- st_read("India_State_Boundary.shp")

# some discrepancy in the state names, hence we modify them
# within the shape file
india_sf <- india_sf %>%
  mutate(
    State_Name = dplyr::recode(
      State_Name,
      
      # Spelling differences
      "Tamilnadu" = "Tamil Nadu",
      "Chhattishgarh" = "Chhattisgarh",
      "Telengana" = "Telangana",
      
      # Union Territories
      "Andaman & Nicobar" = "Andaman & Nicobar Islands",
      "Daman and Diu and Dadra and Nagar Haveli" =
        "Dadra & Nagar Haveli And Daman & Diu",
      
      # Jammu & Kashmir
      "Jammu and Kashmir" = "Jammu & Kashmir"
    )
  )


# in the NFHS-5 data for children
state_data5c <- state_data5c %>%
  mutate(
    State = dplyr::recode(
      State,
      "Dadra & Nagar Haveli and Daman & Diu" =
        "Dadra & Nagar Haveli And Daman & Diu",
      "NCT of Delhi" = "Delhi"
    )
  )

# in the NFHS-5 data for women
state_data5w <- state_data5w %>%
  mutate(
    State = dplyr::recode(
      State,
      "Dadra & Nagar Haveli and Daman & Diu" =
        "Dadra & Nagar Haveli And Daman & Diu",
      "NCT of Delhi" = "Delhi"
    )
  )

unique(india_sf$State_Name)
unique(state_data5w$State)
unique(state_data5c$State)

setdiff(unique(india_sf$State_Name), unique(state_data5c$State))
setdiff(unique(state_data5c$State), unique(india_sf$State_Name))
#"Andaman & Nicobar Islands" %in% state_data4$State #FALSE
#"Lakshadweep" %in% state_data4$State #FALSE 

setdiff(unique(india_sf$State_Name), unique(state_data5w$State))
setdiff(unique(state_data5w$State), unique(india_sf$State_Name))

# left join for the NFHS-5 maps
india_map5c <- india_sf %>%
  left_join(
    state_data5c,
    by = c("State_Name" = "State")
  )

# left join for the NFHS-5 maps
india_map5w <- india_sf %>%
  left_join(
    state_data5w,
    by = c("State_Name" = "State")
  )

india_sf_old <- india_sf %>%
  mutate(
    State_2015 = ifelse(
      State_Name %in% c("Jammu & Kashmir", "Ladakh"),
      "Jammu & Kashmir",
      State_Name
    )
  )

india_sf_old <- india_sf_old %>%
  group_by(State_2015) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

india_map5c <- india_map5c %>%
  rename(State = State_Name)

india_map5w <- india_map5w %>%
  rename(State = State_Name)

india_maps <- bind_rows(
  india_map5c %>% mutate(Category = "ChU5 (6-59 months)"),
  india_map5w %>% mutate(Category = "WRA (15-49 years)")
)

#-------------------------------------------------
# Plot of the PM2.5 exposure for NFHS-5 studies
#-------------------------------------------------

lims<- range(india_map5c$mean_PM2.5, na.rm = TRUE)

ggplot(india_map5c) +
  geom_sf(aes(fill = mean_PM2.5),
          colour = "white",
          linewidth = 0.1) +
  scale_fill_viridis_c(
    option = "plasma",
    limits = lims,
    oob = scales::squish,
    direction=-1,
    name = expression("PM"[2.5]~"("*mu*"g/"*m^3*")")
  ) +
  labs(title = expression("Average " * PM[2.5] * " Exposure by State (NFHS-5)")
  ) + theme_classic()+   
  # gives panel borders by default
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1, "lines"),   # optional: adds a little gap so both panels' right borders are visible
    strip.background = element_rect(fill = "white", colour = "black")
  )

ggsave("Rplot_AvePMExposure_bystate.png", width = 7, height = 6, units = "in",
  dpi = 600)

#india_maps_c<- india_maps %>% subset(Category== "Children 6-59 months")
#------------------------------------------------------------------------------
# Plot of the Anemia prevalence for Children under-5 and women in NFHS-5 study
#------------------------------------------------------------------------------


lims <- range(india_maps$Prevalence, na.rm = TRUE)

ggplot(india_maps) +
  geom_sf(aes(fill = Prevalence),
          colour = "white",
          linewidth = 0.1) +
  facet_wrap(~ Category) +
  scale_fill_viridis_c(
    option = "plasma",
    limits = lims,
    direction = -1,
    oob = scales::squish
  ) +
  labs(title = "Average Anaemia Prevalence by State (NFHS-5)") +
  theme_classic()+   
  # gives panel borders by default
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1, "lines"),   # optional: adds a little gap so both panels' right borders are visible
    strip.background = element_rect(fill = "white", colour = "black")
  )

ggsave("Rplot_AvePrevalence_bystate.png", width = 7, height = 5, units = "in",
       dpi = 600)

#---------------------
# For children only
#---------------------

# compute the range of Anemia Prevalence for children in India
lims <- range(india_maps_c$Prevalence, na.rm = TRUE)

ggplot(india_maps_c) +
  geom_sf(aes(fill = Prevalence),
          colour = "white",
          linewidth = 0.1) +
  scale_fill_viridis_c(
    option = "plasma",
    limits = lims,
    direction = -1,
    oob = scales::squish
  ) +
  labs(title = "Average Anaemia Prevalence") +
  theme_classic()+   
  # gives panel borders by default
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    # optional: adds a little gap so both panels' right borders are visible
    panel.spacing = unit(1, "lines"),   
    strip.background = element_rect(fill = "white", colour = "black")
  )

#----------------------------------------------------------------------------------------
# Plot of the DMI score for Children Under 5 years and WRA (15-49 years) in NFHS-5 study
#----------------------------------------------------------------------------------------

# compute the range of DMIs in India
lims <- range(india_maps$Mean_DMI, na.rm = TRUE)

ggplot(india_maps) +
  geom_sf(aes(fill = Mean_DMI),
          colour = "white",
          linewidth = 0.1) +
  facet_wrap(~ Category) +
  scale_fill_viridis_c(
    name = "DMI",
    option = "plasma",
    limits = lims,
    direction = -1,
    oob = scales::squish
  ) +
  labs(title = "Average DMI(%) by State (NFHS-5)") +
  theme_classic() +   # gives panel borders by default
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "white", colour = "black")
  )

ggsave("Rplot_AveDMI_bystate.png", width = 7, height = 5, units = "in",
       dpi = 600)

#-------------------
# for children only
#-------------------

# compute the range of DMIs for children in India
lims <- range(india_maps_c$Mean_DMI, na.rm = TRUE)

ggplot(india_maps_c) +
  geom_sf(aes(fill = Mean_DMI),
          colour = "white",
          linewidth = 0.1) +
  scale_fill_viridis_c(
    name = "DMI",
    option = "plasma",
    limits = lims,
    direction = -1,
    oob = scales::squish
  ) +
  labs(title = "Average DMI(%)") +
  theme_classic() +   # gives panel borders by default
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "white", colour = "black")
  )

ggsave("Rplot_AveDMI_bystate.png", width = 7, height = 5, units = "in",
       dpi = 600)


#================================
# Plot of the correlation matrix
#================================

# rename the data columns
cor_dat4<- cor(data.frame(d4w[,c("vitaminA", "vitaminB12", "vitaminC", "vitaminD",
                                   "vitaminE", "selenium", "zinc", "iron", "pm2.5")]),
               use= "complete.obs")
  
#cor_dat5<- cor(NFHS_NSS[,nutrients], use ="complete.obs") 
colnames(cor_dat4)<- rownames(cor_dat4)<- c("Vitamin A", "Vitamin B12", "Vitamin C", 
                                            "Vitamin D", "Vitamin E", "Selenium", 
                                            "Zinc", "Iron", "PM(2.5)")

png("Rplot_correlation_percapita_nutrients.png", width = 7, height = 6, 
    units = "in", res = 600)

corrplot(cor_dat4, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45, tl.cex = 0.9,
         col = colorRampPalette(c("#4C72B0", "white", "#DD8452"))(200),
         title = "NFHS-4 / NSSO-68 (per-capita nutrients)",
         mar = c(0, 0, 2, 0))

dev.off()


#-------------- 
# For children
#--------------

vcm<- survey:::vcov.svyglm(fit_c)[c("PM2.5", "PM2.5:DMIscaled"), 
                                  c("PM2.5", "PM2.5:DMIscaled")]
cm<- coef(fit_c)[c("PM2.5", "PM2.5:DMIscaled")]

conf_band <- function(x, coef_m, vcov_m) {
  
  est <- coef_m[1] + coef_m[2] * x
  se <- sqrt(matrix(c(1, x), nrow = 1) %*% vcov_m %*% matrix(c(1, x), ncol = 1))
  
  lb <- exp(est - 1.96 * se)
  ub <- exp(est + 1.96 * se)
  
  return(c("estimate" = exp(est), "lower" = lb, "upper" = ub))
}

cf_c <- t(sapply(d5c$DMIscaled, FUN= conf_band, coef_m = cm,
               vcov_m = vcm))
colnames(cf_c)<- c("Odds", "Lower", "Upper")
plot_dfc<- cbind.data.frame(DMI = d5c$DMIscaled, cf_c)

pc<- ggplot(plot_dfc, aes(x = DMI, y = Odds)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper),
              alpha = 0.2,
              col = "red") +
  theme_bw() +
  geom_boxplot(
    data = plot_dfc,
    aes(y = 0.982, x = DMI),
    outlier.size = 0.3,
    width = 0.004
  ) +
  scale_y_continuous(limits = c(0.98, 1.01)) +
  labs(x = "DMI", y = TeX(r'(Coefficient for PM$_{2.5}$)'),
       title = "ChU5 (6-59 months)")

pc

#------------
# For women
#------------

vcm<- survey:::vcov.svyglm(fit_w) [c("PM2.5", "PM2.5:DMIscaled"), 
                                  c("PM2.5", "PM2.5:DMIscaled")]
cm<- coef(fit_w)[c("PM2.5", "PM2.5:DMIscaled")]

conf_band <- function(x, coef_m, vcov_m) {
  
  est <- coef_m[1] + coef_m[2] * x
  se <- sqrt(matrix(c(1, x), nrow = 1) %*% vcov_m %*% matrix(c(1, x), ncol = 1))
  
  lb <- exp(est - 1.96 * se)
  ub <- exp(est + 1.96 * se)
  
  return(c("estimate" = exp(est), "lower" = lb, "upper" = ub))
}

cf <- t(sapply(d5w$DMIscaled, FUN = conf_band, coef_m = cm,
               vcov_m = vcm))
colnames(cf)<- c("Odds", "Lower", "Upper")
plot_df<- cbind.data.frame(DMI = d5w$DMIscaled, cf)

pw<- ggplot(plot_df, aes(x = DMI, y = Odds)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper),
              alpha = 0.2,
              col = "red") +
  theme_bw() +
  geom_boxplot(
    data = plot_df,
    aes(y = 0.982, x = DMI),
    outlier.size = 0.3,
    width = 0.004
  ) +
  scale_y_continuous(limits = c(0.98, 1.01)) +
  labs(x = "DMI", y = TeX(r'(Coefficient for PM$_{2.5}$)'),
       title= "WRA (15-49 years)")
pw
pc+pw
ggsave("Rplot_predodds_vsDMI.png", width= 8, height= 5.5, units= "in", dpi= 300)

#========================
# Densities of nutrients
#========================

# Select only nutrient variables
d4c_long <- d4c %>%
  select(vitaminA, vitaminB12,vitaminC, vitaminD, vitaminE, selenium, zinc) %>%
  
  pivot_longer(everything(), names_to= "Nutrient", values_to= "Value") %>%
  mutate(Survey = "NFHS-4")

d5c_long <- d5c %>%
  select(vitaminA, vitaminB12,vitaminC, vitaminD, vitaminE, selenium, zinc) %>%
  pivot_longer(everything(), names_to= "Nutrient", values_to= "Value") %>%
  mutate(Survey = "NFHS-5")

plot_data <- bind_rows(d4c_long, d5c_long)

# Order of panels
plot_data$Nutrient <- factor(
  plot_data$Nutrient,
  levels = c("vitaminA", "vitaminB12", "vitaminC",
             "vitaminD", "vitaminE", "selenium", "zinc")
)

# Labels for the facet strips
nutrient_labels <- c(
  vitaminA  = "Vitamin~A~(mu*g)",
  vitaminB12 = "Vitamin~B[12]~(mu*g)",
  vitaminC  = "Vitamin~C~(mg)",
  vitaminD  = "Vitamin~D~(mu*g)",
  vitaminE  = "Vitamin~E~(mg)",
  selenium  = "Selenium~(mu*g)",
  zinc      = "Zinc~(mg)"
)

ggplot(plot_data,
       aes(Value,
           colour = Survey,
           fill = Survey)) +
  geom_density(alpha = 0.25, linewidth = 0.3, na.rm = TRUE) +
  facet_wrap(
    ~Nutrient,
    scales = "free",
    ncol = 3,
    labeller = labeller(
      Nutrient = as_labeller(nutrient_labels, label_parsed)
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Per-capita nutrient consumption",
    y = "Density",
    x= NULL
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("Rplot_nutrient_distributions.png", width = 9, height = 7, 
       units = "in", dpi= 600)


#------------------------------------
# Barplot comparing the DMI weights 
#------------------------------------

# DMI weights for women
dmi_w<- read.csv("DMI_wra.csv")[,-c(1,2)]
dmi_w<- ifelse (is.na(dmi_w), 0, dmi_w)

dmi_c<- read.csv("DMI_ChU5.csv")[,-1]
dmi_c<- ifelse(is.na(dmi_c), 0, dmi_c)*100

dmi_w <- dmi_w[dmi_w != 0]
dmi_c <- dmi_c[dmi_c != 0]

# Create a long-format data frame
plot_data <- tibble(
  Nutrient = c("Vitamin A", "Vitamin B12", "Vitamin D",
               "Vitamin E", "Selenium", "Zinc"),
  `WRA (15-49 years)` = dmi_w,
  `ChU5 (6-59 months)` = dmi_c
) %>%
  pivot_longer(cols = c(`WRA (15-49 years)`, `ChU5 (6-59 months)`),
               names_to = "Group",
               values_to = "Weight")


ggplot(plot_data,
       aes(x = Nutrient, y = Weight, fill = Group)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65) +
  labs(x = NULL,
       y = "DMI Weight",
       fill = NULL) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )+ geom_text(aes(label = sprintf("%.0f%%", Weight)),
          position = position_dodge(width = 0.75),
          vjust = -0.3,
          size = 3.5)+
  scale_y_continuous(limits = c(0, 40))
 
ggsave("Rplot_DMInutrient_comparison.png", width = 7, height = 4, units = "in",
             dpi = 600)

#-------------------------------------------------------------
# Table for the ODDS of covariates in the Logistic regression
#-------------------------------------------------------------

OR_svyglm<- function(model, terms, delta= 1) {
  beta <- coef(model)[terms]
  pvalue<- summary(model)$coefficients[terms, "Pr(>|t|)"]
  ci <- confint(model)[terms,]
  OR <- exp(delta* beta)
  ci_OR <- exp(delta* ci)
  pct <- ifelse(OR>= 1, (OR - 1) * 100, (1 - OR) * 100)
  results<- list(OR = OR, CI_low = ci_OR, CI_high = ci_OR, pvalue= pvalue)
  
  results_df <- data.frame(
    OR      = as.numeric(results$OR),
    CI_low  = results$CI_low[, "2.5 %"],
    CI_high = results$CI_high[, "97.5 %"],
    p_value = as.numeric(results$pvalue)
  )
  return(round(results_df,3))
}


results <- OR_svyglm(fit_c, c("PM2.5", "PM2.5:DMIscaled"), delta= 10)
results

results <- OR_svyglm(fit_w, c("PM2.5", "PM2.5:DMIscaled"), delta= 10)
results

results <- OR_svyglm(fit_c, names(fit_c$coefficients), delta= 1)
results

results <- OR_svyglm(fit_w, names(fit_w$coefficients), delta= 1)
results

