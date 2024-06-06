library(survival)
library(survminer)

library(dplyr)

source("./signature_functions.r")

## get the 50% survival time
get_50_survival_time <- function(fit){

df <- surv_summary((fit))

high_df <- df[df$Feature == "high",]
low_df <- df[df$Feature == "low",]

high_05sample = which(high_df$surv <= 0.5)
low_05sample = which(low_df$surv <= 0.5)

if (length(high_05sample) == 0 ){
    time_at_50_percent_high = NA
}
else {
   time_at_50_percent_high = high_df[high_05sample[1],"time"]
}
if (length(low_05sample) == 0 ){
    time_at_50_percent_low = NA
}
else {
   time_at_50_percent_low = low_df[low_05sample[1],"time"]
}

    print(paste("50% survival time for 'high' group:", time_at_50_percent_high))
    print(paste("50% survival time for 'low' group:", time_at_50_percent_low))

    return(NULL)
}

## validation
val_df <- read.csv(file = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/val_all_classes_df.csv", header = T)
val_clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")

val_df$pred_1_0 <- val_df$Class_1_Prob / val_df$Class_0_Prob

val_df$RFS_time <- val_clinical[match(val_df$PID,val_clinical$PID),"RFS_time"]
val_df$RFS_status <- val_clinical[match(val_df$PID,val_clinical$PID),"RFS_status"]
head(val_df)

rownames(val_df) <- val_df[,1]

# Iterating over each column
save_path = "./survival_results/"
if(!dir.exists(save_path)){
    dir.create(save_path,recursive = T)
}

variable <- "pred_1_0"

val_df_ <- val_df[,c(variable,"RFS_time","RFS_status")]
colnames(val_df_) <- c("variable","RFS_time","RFS_status")
cutpoint <- surv_cutpoint(data = val_df_, time = "RFS_time", event = "RFS_status", variables = "variable")
cutpoint_val <- summary(cutpoint)$cutpoint[[1]]

val_df_$Feature <- ifelse(val_df_$variable > cutpoint_val, "high", "low")
val_df_[order(val_df_$variable,decreasing = T),]

fit <- survfit(Surv(RFS_time, RFS_status) ~ Feature, data = val_df_)
cox.fit <- coxph(Surv(RFS_time, RFS_status) ~ Feature, data = val_df_)

c_index = summary(cox.fit)$concordance ## c-index
print(c_index)

# Extract survival times and probabilities
get_50_survival_time(fit)

p <- ggsurvplot(fit,
                data = val_df_,
                linetype = c("solid", "solid"),
                surv.median.line = "hv", 
                surv.scale = "percent",
                pval = TRUE, 
                risk.table = TRUE,
                conf.int = TRUE, 
                conf.int.alpha = 0.1, 
                conf.int.style = "ribbon",
                risk.table.y.text = TRUE,
                palette = c("#3300CC", "#CC3300"),
                xlab = "Recurrence time"
                )
pdf(paste0(save_path,"Survival analysis on validation via ",variable,".pdf"),width=8,height=6)
print(p)
dev.off()

## test
test_df <- read.csv(file = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/test_all_classes_df.csv", header = T)
test_clinical <- read.csv("/mnt/data/lyx/IMC/WSI_HE/WSI_Testing/WSI_Test_clinical.csv")

test_df$pred_1_0 <- test_df$Class_1_Prob / test_df$Class_0_Prob

test_df$PID <- paste0("P",gsub("-","_",test_df$PID))
test_clinical$PID <- paste0("P",test_clinical$PID)
test_clinical$PID <- gsub(" ","",test_clinical$PID)

test_df <- left_join(test_df,test_clinical,by = "PID")
test_df <- na.omit(test_df)

head(test_df,20)
rownames(test_df) <- test_df[,1]

# Iterating over each column
save_path = "./survival_results/"
if(!dir.exists(save_path)){
    dir.create(save_path,recursive = T)
}

variable <- "pred_1_0"

test_df_ <- test_df[,c(variable,"RFS_time","RFS_status")]
write.csv(test_df_, file = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/test_result.csv")
test_df_ <- na.omit(test_df_)

# test_df_ <- read.csv("/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/test_result.csv",row.names=1)
colnames(test_df_) <- c("variable","RFS_time","RFS_status")

cutpoint <- surv_cutpoint(data = test_df_, time = "RFS_time", event = "RFS_status", variables = "variable")
cutpoint_test <- summary(cutpoint)$cutpoint[[1]]

test_df_$Feature <- ifelse(test_df_$variable > cutpoint_test, "high", "low")
test_df_[order(test_df_$variable,decreasing = T),]

fit <- survfit(Surv(RFS_time, RFS_status) ~ Feature, data = test_df_)
cox.fit <- coxph(Surv(RFS_time, RFS_status) ~ Feature, data = test_df_)

c_index = summary(cox.fit)$concordance ## c-index
print(c_index)

get_50_survival_time(fit)

p <- ggsurvplot(fit,
                data = test_df_,
                linetype = c("solid", "solid"),
                surv.median.line = "hv", 
                surv.scale = "percent",
                pval = TRUE, 
                risk.table = TRUE,
                conf.int = TRUE, 
                conf.int.alpha = 0.1, 
                conf.int.style = "ribbon",
                risk.table.y.text = TRUE,
                palette = c("#3300CC", "#CC3300"),
                xlab = "Recurrence time"
                )
pdf(paste0(save_path,"Survival analysis on test via ",variable,".pdf"),width=8,height=6)
print(p)
dev.off()

## Plus with other clinical information
### val 
if(T){
val_df <- read.csv(file = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/val_all_classes_df.csv", header = T)
val_clinical <- read.csv("/mnt/data/lyx/IMC/clinical.csv")
save_path = "./survival_results/"
if(!dir.exists(save_path)){
    dir.create(save_path,recursive = T)
}

val_df$pred_1_0 <- val_df$Class_1_Prob / val_df$Class_0_Prob
val_df <- val_df[,c("PID","pred_1_0")]

clinical_features <- c(
    "PID", "Recurrence_site", 
    "fong_score", "Gender", "Age", "KRAS_mutation",
    "TBS", "CRLM_number", "CRLM_size", "Liver_involvement_num", 
    "CRC_site", "CEA", "CA199", "Differential_grade", "T_stage", "Lymph_positive",
    "RFS_status", "RFS_time"
)

val_clinical <- val_clinical[,clinical_features]

val_df <- left_join(val_df,val_clinical,by="PID")
head(val_df)

rownames(val_df) <- val_df[,1]
val_df <- val_df[,-1]

#### COX
unicox_result <- MultipleUniCOX(val_df,UniCOX = T)
ForestPlot(unicox_result,savePath = paste0(save_path,"trainset-uniCox.pdf"))
multicox_result <- MultipleUniCOX(val_df,UniCOX = F)
ForestPlot(multicox_result,savePath = paste0(save_path,"trainset-mulCox.pdf"))

}

### test 
if(T){
test_df <- read.csv(file = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/test_all_classes_df.csv", header = T)
test_clinical <- read.csv("/mnt/data/lyx/IMC/WSI_HE/WSI_Testing/WSI_Test_clinical.csv")
save_path = "./survitest_results/"
if(!dir.exists(save_path)){
    dir.create(save_path,recursive = T)
}

test_df$pred_1_0 <- test_df$Class_1_Prob / test_df$Class_0_Prob
test_df <- test_df[,c("PID","pred_1_0")]
test_df$PID <- paste0("P",gsub("-","_",test_df$PID))
test_clinical$PID <- paste0("P",test_clinical$PID)
test_clinical$PID <- gsub(" ", "", test_clinical$PID)

clinical_features <- c(
    "PID",
    "fong_score", "Gender", "Age", "KRAS_mutation",
    "TBS", "CRLM_number", "CRLM_size", "Liver_involvement_num", 
    "CRC_site",   "Differential_grade", "T_stage", "Lymph_positive", #"CEA","CA199",
    "RFS_status", "RFS_time"
)

test_clinical <- test_clinical[,clinical_features]

test_df <- left_join(test_df, test_clinical, by = "PID")
test_df <- na.omit(test_df)
head(test_df)

rownames(test_df) <- test_df[,1]
test_df <- test_df[,-1]

test_df <- test_df[-match("P18S19932_4",rownames(test_df)),]

# write.csv(test_df,"test_df.csv")

#### COX
unicox_result <- MultipleUniCOX(test_df,UniCOX = T)
ForestPlot(unicox_result,savePath = paste0(save_path,"testset-uniCox.pdf"))
multicox_result <- MultipleUniCOX(test_df,UniCOX = F)
ForestPlot(multicox_result,savePath = paste0(save_path,"testset-mulCox.pdf"))

}


##########################
## the difference between TAR and chemotherapy only in BDME score high group
library(readxl)
library(survminer)

## train
train_ <- read_excel("/home/lenislin/mnt_16T/ProjectData/IMC_CRC/WSI_HE/training.xlsx")
head(train_)
train_ <- as.data.frame(train_)
train_$Type <- as.factor(train_$Type)
colnames(train_) <- c("RFS_status","RFS_time","Type")

val_ <- read_excel("/home/lenislin/mnt_16T/ProjectData/IMC_CRC/WSI_HE/validating.xlsx")
val_ <- as.data.frame(val_)
val_$Type <- as.factor(val_$Type)
head(val_)

train_fit <- survfit(Surv(RFS_time, RFS_status) ~ Type, data = train_)
val_fit <- survfit(Surv(RFS_time, RFS_status) ~ Type, data = val_)

save_path <- "/home/lenislin/mnt_16T/Project/IMC-CRC/analysis/HE/"

p <- ggsurvplot(train_fit,
                data = train_,
                linetype = c(1, 1, 2, 2),
                surv.median.line = "hv", 
                surv.scale = "percent",
                pval = TRUE, 
                risk.table = TRUE,
                conf.int = F, 
                conf.int.alpha = 0.1, 
                conf.int.style = "ribbon",
                risk.table.y.text = TRUE,
                palette = c("#CC0C00FF", "#5C88DAFF", "#95CC5EFF","#F79D1EFF"),
                xlab = "Recurrence time"
                )

pdf(paste0(save_path,"TAR vs. Che on train set.pdf"),width=10,height=7.5)
print(p)
dev.off()

p <- ggsurvplot(val_fit,
                data = val_,
                linetype = c(1, 1, 2, 2),
                surv.median.line = "hv", 
                surv.scale = "percent",
                pval = TRUE, 
                risk.table = TRUE,
                conf.int = F, 
                conf.int.alpha = 0.1, 
                conf.int.style = "ribbon",
                risk.table.y.text = TRUE,
                palette = c("#CC0C00FF", "#5C88DAFF", "#95CC5EFF","#F79D1EFF"),
                xlab = "Recurrence time"
                )

pdf(paste0(save_path,"TAR vs. Che on val set.pdf"),width=10,height=7.5)
print(p)
dev.off()
