#load lavaan and semTools packages to perform confirmatory factor analysis and measure composite reliability
library(lavaan)    
library(semTools)

#load robustbase and car packages to perform univariate and multivariate analyses
library(robustbase)
library(car)

#load the data file
IADQ <- read.csv("IADQ dataset.csv")

#Confirmatory factor analyses of continuous data:

#Model 1: 1-factor solution
cont_mod1 <-'f =~ ajd10 + ajd11 + ajd12 + ajd13 + ajd14 + ajd15'
cont_mod1_fit <- cfa(cont_mod1, data=IADQ, estimator = "MLR")
summary(cont_mod1_fit, standardized = TRUE, fit.measures = TRUE)

#Model 2: 2-factor solution
cont_mod2 <- 'f1 =~ ajd10 + ajd11 + ajd12
              f2 =~ ajd13 + ajd14 + ajd15
              #Correlation between factors
              f1 ~~ f2'
cont_mod2_fit <- cfa(cont_mod2, data=IADQ, estimator = "MLR")
summary(cont_mod2_fit, standardized = TRUE, fit.measures = TRUE)

#Model 3: modified 2-factor solution
cont_mod3 <- 'f1 =~ ajd10 + ajd11
              f2 =~ ajd12 + ajd13 + ajd14 + ajd15
              #Correlation between factors
              f1 ~~ f2'
cont_mod3_fit <- cfa(cont_mod3, data=IADQ, estimator = "MLR")
summary(cont_mod3_fit, standardized = TRUE, fit.measures = TRUE)

#Model 4: overlapping 2-factor solution
cont_mod4 <- 'f1 =~ ajd10 + ajd11 + ajd12
              f2 =~ ajd12 + ajd13 + ajd14 + ajd15
              #Correlation between factors
              f1 ~~ f2'
cont_mod4_fit <- cfa(cont_mod4, data=IADQ, estimator = "MLR")
summary(cont_mod4_fit, standardized = TRUE, fit.measures = TRUE)


#Confirmatory factor analyses for endorsement rates:

#Model 1: 1-factor solution
endo_mod1 <- 'fe =~ ajd10e + ajd11e + ajd12e + ajd13e + ajd14e + ajd15e'
endo_mod1_fit <- cfa(endo_mod1, data=IADQ, estimator = "WLSMV", ordered = c("ajd10e", "ajd11e", "ajd12e", "ajd13e", "ajd14e", "ajd15e"))
summary(endo_mod1_fit, standardized = TRUE, fit.measures = TRUE)

#Model 2: 2-factor solution
endo_mod2 <- 'f1e =~ ajd10e + ajd11e + ajd12e
              f2e =~ ajd13e + ajd14e + ajd15e
              #Correlation between factors
              f1e ~~ f2e'
endo_mod2_fit <- cfa(endo_mod2, data=IADQ, estimator = "WLSMV", ordered = c("ajd10e", "ajd11e", "ajd12e", "ajd13e", "ajd14e", "ajd15e"))
summary(endo_mod2_fit, standardized = TRUE, fit.measures = TRUE)

#Model 3: modified 2-factor solution
endo_mod3 <- 'f1e =~ ajd10e + ajd11e
              f2e =~ ajd12e + ajd13e + ajd14e + ajd15e
              #Correlation between factors
              f1e ~~ f2e'
endo_mod3_fit <- cfa(endo_mod3, data=IADQ, estimator = "WLSMV", ordered = c("ajd10e", "ajd11e", "ajd12e", "ajd13e", "ajd14e", "ajd15e"))
summary(endo_mod3_fit, standardized = TRUE, fit.measures = TRUE)

#Model 4: overlapping 2-factor solution
endo_mod4 <- 'f1e =~ ajd10e + ajd11e + ajd12e
              f2e =~ ajd12e + ajd13e + ajd14e + ajd15e
              #Correlation between factors
              f1e ~~ f2e'
endo_mod4_fit <- cfa(endo_mod4, data=IADQ, estimator = "WLSMV", ordered = c("ajd10e", "ajd11e", "ajd12e", "ajd13e", "ajd14e", "ajd15e"))
summary(endo_mod4_fit, standardized = TRUE, fit.measures = TRUE)


#Calculate composite reliability for preoccupation and failure to adapt in model 2
comp_rel <- compRelSEM(cont_mod2_fit)
print(comp_rel, digits = 3)


#Intrinsic concurrent validity: Stressor scale ---> preoccupation, failure to adapt
se_po <- lmrob(po ~ se, data = IADQ)
summary(se_po)
se_fta <- lmrob(fta ~ se, data = IADQ)
summary(se_fta)
se_multi <- manova(cbind(po, fta) ~ se, data = IADQ)
Anova(se_multi, type = "III", robust = TRUE)

#Extrinsic concurrent validity: Preoccupation ---> PHQ-9, GAD-7
po_phq9 <- lmrob(phq9 ~ po, data = IADQ)
summary(po_phq9)
po_gad7 <- lmrob(gad7 ~ po, data = IADQ)
summary(po_gad7)
po_multi <- manova(cbind(phq9, gad7) ~ po, data = IADQ)
Anova(po_multi, type = "III", robust = TRUE)

#Extrinsic concurrent validity: Failure to adapt ---> PHQ-9, GAD-7
fta_phq9 <- lmrob(phq9 ~ fta, data = IADQ)
summary(fta_phq9)
fta_gad7 <- lmrob(gad7 ~ fta, data = IADQ)
summary(fta_gad7)
fta_multi <- manova(cbind(phq9, gad7) ~ fta, data = IADQ)
Anova(fta_multi, type = "III", robust = TRUE)

#Extrinsic concurrent validity: Preoccupation, Failure to adapt ---> PHQ-9, GAD-7
iadq_phq9 <- lmrob(phq9 ~ po + fta, data = IADQ)
summary(iadq_phq9)
iadq_gad7 <- lmrob(gad7 ~ po + fta, data = IADQ)
summary(iadq_gad7)
iadq_multi <- manova(cbind(phq9, gad7) ~ po + fta, data = IADQ)
Anova(iadq_multi, type = "III", robust = TRUE)

#Extrinsic concurrent validity: Preoccupation + Failure to adapt ---> PHQ-9, GAD-7
total_phq9 <- lmrob(phq9 ~ total, data = IADQ)
summary(total_phq9)
total_gad7 <- lmrob(gad7 ~ total, data = IADQ)
summary(total_gad7)
total_multi <- manova(cbind(phq9, gad7) ~ total, data = IADQ)
Anova(total_multi, type = "III", robust = TRUE)
