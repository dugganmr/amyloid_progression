#primary analyses (linear regression)
#adaptable for specific variables (e.g., columns #s) and file names
library(tidyverse)       
library(haven)           
library(readxl)          
library(data.table)
library(nlme)


#load predictor data (SomaScan proteins)
df_a<-readRDS(file = "filename")

#WGCNA
powers = c(seq(1,20,by=1))
sft = pickSoftThreshold(df_a[,7:948], powerVector = powers, networkType="signed", verbose = 2) 
par(mfrow = c(1,2))
cex1 = 0.9
#  Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
#export 4.5x6.5
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
mergingThresh<-0.10 
net<-blockwiseModules(df_a[,7:948], 
                      corType = "bicor", 
                      power = 14, 
                      networkType = "signed", 
                      TOMType = "signed", 
                      minModuleSize = 5,
                      mergeCutHeight = mergingThresh,
                      deepSplit = 2, 
                      pamRespectsDendro = FALSE)
MEsAutomatic<-net$MEs
rownames(MEsAutomatic) <- 1:nrow(MEsAutomatic)
df_a<-cbind(df_a, MEsAutomatic)
#load outcome data set (PET, ADRD biomarkers etc.,)
df_b <- read_sas("filename")
#merge datasets
df1 <-left_join(df_b, df_a, by = c("ID", "visit"))
#exclude participants with neurovascular conditions that could affect brain structure or function (e.g., strokes)
neuro_comp<-read_excel("filename")
df2 <- left_join(df1, neuro_comp, by = c("ID", "visit"))
df2 <- filter(df2, exclude == "0")
#exclude participants with cognitive impairment
cog_impair <-read_excel("filename")
df2 <- left_join(df2, cog_impair, by = c("ID", "visit"))
df2 <- filter(df2, dx == "0")
#load covariates
covar <- read_sas("filename")
df2 <-left_join(df2, covar, by = c("ID", "visit"))

#calculate co-morbidity index
df2$cindex_freq <-rowSums(df2[,c("obesity0", "bp_risk0","diabetes0",
"cancer0", "hid0", "chf0", "copd0", "ckd0")])
df2$cindex_denom <- 8-(is.na(df2$obesity0) + is.na(df2$bp_risk0) + 
is.na(df2$diabetes0) + is.na(df2$cancer0) + 
is.na(df2$hid0) + is.na(df2$chf0) + is.na(df2$copd0)+ is.na(df2$ckd0))
df2$cindex_percent <- (df2$cindex_freq/df2$cindex_denom)*100

#Generate Slopes
outcome <- lme(outcome ~ Time, data = df2, random = ~ Time | ID)
outcome_slope <- coef(outcome) #
outcome_slope <- sapply( outcome_slope, as.numeric )
outcome_slope[,"outcome"] <- as.numeric(scale(outcome_slope[,"outcome"]))
soma2<-left_join(soma2,outcome_slope, by=c("ID"))


#format for analyses
length(df2$ID) # observations
length(unique(df2$ID)) #  participants 
df2<-as.data.frame(df2)
df2$race <- as.factor(df2$race)
df2$race <- relevel(df2$race, ref = "0")
df2$sex <- as.factor(df2$sex)
df2$sex <- relevel(df2$sex, ref = "0")
df2$apoe <- as.factor(df2$apoe)
df2$apoe <- relevel(df2$apoe, ref = "0")
df2$infection <- as.factor(df2$infection )
df2$infection  <- relevel(df2$infection, ref = "0")
#select predictors
df_vars <- as.data.frame(df2)[c(00:00)]

#linear regression models (loop)
output <- list()
for (i in 1:length(df_vars)){
  var=colnames(df_vars)[i]
  output[[var]] = list() 
  vars <- df_vars[,i]
  fit <- lm(outcome ~ vars + age + sex + race + educ_years + cindex_percent + PIB+apoe4+batch,
            data = soma2, na.action=na.omit)
  output[[var]]=summary(fit)$coefficients[ 2,]
}
output<-as.data.frame(output)
write_xlsx(output,"results.xlsx")