#the following code is provided in a general format (e.g., predictor variables, response variables, columns #s etc.,)
#that is readily adaptable for specific analyses
library(stats)
library(nlme)
library(performance)
library(survival)
library(Seurat)
library(writexl)

#BLSA
#load predictor/outcome datasets (PET, SomaScan proteomic, pQTL, ADRD biomarker etc.,)
df1<-readRDS(file = "filename")

#calculate co-morbidity index
df1$cindex_freq <-rowSums(df1[,c("obesity0", "bp_risk0","diabetes0",
"cancer0", "hid0", "chf0", "copd0", "ckd0")])
df1$cindex_denom <- 8-(is.na(df1$obesity0) + is.na(df1$bp_risk0) + 
is.na(df1$diabetes0) + is.na(df1$cancer0) + 
is.na(df1$hid0) + is.na(df1$chf0) + is.na(df1$copd0)+ is.na(df1$ckd0))
df1$cindex_percent <- (df1$cindex_freq/df1$cindex_denom)*100

#generate slopes (e.g., PiB)
outcome <- lme(outcome ~ Time, data = df1, random = ~ Time | ID)
check_model(outcome)
outcome_slope <- coef(outcome)
colnames(outcome_slope )[colnames(outcome_slope) =="Time"] <- "outcome"
outcome_slope <- as.data.frame(outcome_slope)[c(2)]
outcome_slope <- tibble::rownames_to_column(outcome_slope, "ID")
outcome_slope <- sapply( outcome_slope, as.numeric )
df1<-left_join(df1,outcome_slope, by=c("ID"))

#format for analyses
length(df1$ID) # observations
length(unique(df1$ID)) #  participants 
df1<-as.data.frame(df1)
df1$race <- as.factor(df1$race)
df1$race <- relevel(df1$race, ref = "0")
df1$sex <- as.factor(df1$sex)
df1$sex <- relevel(df1$sex, ref = "0")
df1$apoe <- as.factor(df1$apoe)
df1$apoe <- relevel(df1$apoe, ref = "0")
#select predictors
df_vars <- as.data.frame(df1)[c(00:00)]

#predictors x dichotomous response variable (e.g., PIB status)
output <- list()  # Create empty list
for (i in 1:length(df_vars)){
  var=colnames(df_vars)[i]
  output[[var]] = list() # Create one entry for each variable
  vars <- df_vars[,i]
  fit <- glm(outcome ~ vars + age_covary + Sex + Race + educ_years + apoe4 + cindex_percent+batch
             , data = df1, family = "binomial")
  output[[var]]=summary(fit)$coefficients[c(2), ]
}
output<-as.data.frame(output)
write_xlsx(output,"results.xlsx")

#predictors x continuous response variable (e.g., PiB slopes)
output <- list()
for (i in 1:length(df_vars)){
  var=colnames(df_vars)[i]
  output[[var]] = list() 
  vars <- df_vars[,i]
  fit <- lm(outcome ~ vars + age + sex + race + educ_years + cindex_percent + PIB+apoe4+batch,
            data = df1, na.action=na.omit)
  output[[var]]=summary(fit)$coefficients[ 2,]
}
output<-as.data.frame(output)
write_xlsx(output,"results.xlsx")

#for PiB status stratified analyses, filter PiB==1 or PiB==0
#for sex stratified analyses, filter sex==1 or sex==0
#for ADRD biomarker analyses, include "eGFR"
#for minimally adjusted models, remove "race" and "educ_years"



#ARIC
#load predictor/outcome datasets
df1<-readRDS(file = "ARIC_filename")
#predictors x dichotomous response variable (e.g., PIB status)
output <- list()  # Create empty list
for (i in 1:length(df_vars)){
  var=colnames(df_vars)[i]
  output[[var]] = list() # Create one entry for each variable
  vars <- df_vars[,i]
  fit <- glm(outcome ~ vars + age + sex + race_center + education + E4 + BMI+hyperten+smk+diabetes+eGFR
             , data = df1, family = "binomial")
  output[[var]]=summary(fit)$coefficients[c(2), ]
}
output<-as.data.frame(output)
write_xlsx(output,"results.xlsx")

#predictors x dichotomous response variable (e.g., incident dementia risk)
output <- list()  # Create empty list
for (i in 1:length(df_vars)){
  var=colnames(df_vars)[i]
  output[[var]] = list() # Create one entry for each variable
  vars <- df_vars[,i]
  fit <- coxph(Surv(time_toDx) ~ vars + age + sex + race_center + education + E4 + BMI+hyperten+smk+diabetes+eGFR
             , data = df1)
  output[[var]]=summary(fit)$coefficients[c(2), ]
}
output<-as.data.frame(output)
write_xlsx(output,"results.xlsx")



#Microglia
#load predictor/outcome datasets
df1<-readRDS(file = "microglia_filename")
results <- FindMarkers(df1, ident.1 = 1, assay=SCT, test.use="wilcox", min.pct=0.01, logfc.threshold=0.1)
write_xlsx(results,"results.xlsx")