
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(broom)
library(MendelianRandomization)
library(psych)
library(openxlsx)


## df.smk: dataframe for SS
## df.cpd: dataframe for cpd


### Main MR analysis:
### Analysis using MendelianRandomization package: https://cran.r-project.org/web/packages/MendelianRandomization/index.html

#### SS

##### Not account for pairwise correlation of IVs
MRsmk<- mr_allmethods(mr_input(bx = df.smk$BETA_exp, bxse=df.smk$SE_exp ,
                               by = df.smk$BETA_out, byse = df.smk$SE_out ), method = "all") #
MRsmk

##### Only perform MR (MR-IVW, MR-Egger) with generalized weighted regression model:
MRsmk_cor<- mr_allmethods(mr_input(bx = df.smk$BETA_exp, bxse=df.smk$SE_exp ,
                                   by = df.smk$BETA_out, byse = df.smk$SE_out, correlation = as.matrix(smk.ldmx) ), method = "all")
MRsmk_cor

##### Only perform MR-IVW with generalized weighted regression model:
MRsmk_cor_IVW<- mr_ivw(mr_input(bx = df.smk$BETA_exp, bxse=df.smk$SE_exp ,
                                   by = df.smk$BETA_out, byse = df.smk$SE_out, correlation = as.matrix(smk.ldmx) ) )
MRsmk_cor_IVW


#### CPD
##### Not account for pairwise correlation of IVs
MRcpd<- mr_allmethods(mr_input(bx = df.cpd$BETA_exp, bxse=df.cpd$SE_exp ,
                               by = df.cpd$BETA_out, byse = df.cpd$SE_out ), method = "all") #
MRcpd

##### Only perform MR (MR-IVW, MR-Egger) with generalized weighted regression model:
MRcpd_cor<- mr_allmethods(mr_input(bx = df.cpd$BETA_exp, bxse=df.cpd$SE_exp ,
                                   by = df.cpd$BETA_out, byse = df.cpd$SE_out, correlation = as.matrix(cpd.ldmx) ), method = "all")
MRcpd_cor


##### Only perform MR-IVW with generalized weighted regression model:
MRcpd_cor_IVW<- mr_ivw(mr_input(bx = df.cpd$BETA_exp, bxse=df.cpd$SE_exp ,
                                   by = df.cpd$BETA_out, byse = df.cpd$SE_out, correlation = as.matrix(cpd.ldmx) ) )
MRcpd_cor_IVW




### Sensitivity analysis:

### Leave-one-out

#### SS
MRsmk_loo=mr_loo(mr_input(bx = df.smk.1$BETA_exp, bxse=df.smk.1$SE_exp,
                          by = df.smk.1$BETA_out, byse = df.smk.1$SE_out  ) ) 
MRsmk_loo


#### CPD
MRcpd_loo=mr_loo(mr_input(bx = df.cpd$BETA_exp, bxse=df.cpd$SE_exp ,
                          by = df.cpd$BETA_out, byse = df.cpd$SE_out) )
MRcpd_loo



### contamination mixture method ###

#### SS
mr_conmix(mr_input(bx = df.smk$BETA_exp, bxse=df.smk$SE_exp ,
                   by = df.smk$BETA_out, byse = df.smk$SE_out) )

#### CPD
mr_conmix(mr_input(bx = df.cpd$BETA_exp, bxse=df.cpd$SE_exp ,
                   by = df.cpd$BETA_out, byse = df.cpd$SE_out) )




### MR-PRESSO ###
### Package: https://github.com/rondolab/MR-PRESSO 

if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

#### SS
SMKpresso=mr_presso(BetaOutcome = "BETA_out", BetaExposure = "BETA_exp", SdOutcome = "SE_out", SdExposure = "SE_exp", 
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = df.smk.1, NbDistribution = 1000,  SignifThreshold = 0.05)
#### CPD
CPDpresso=mr_presso(BetaOutcome = "BETA_out", BetaExposure = "BETA_exp", SdOutcome = "SE_out", SdExposure = "SE_exp", 
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = df.cpd, NbDistribution = 1000,  SignifThreshold = 0.05)




### MR-Mix ###
### Package: https://github.com/gqi/MRMix ###

devtools::install_github("gqi/MRMix")
library(MRMix)
library(dplyr)

## Obtain the MAF for IVs:


### SMK:

smk_std = with(df.smk.2, standardize(BETA_exp, BETA_out, SE_exp, SE_out, xtype="binary", ytype = "continuous", nx = NMISS_exp, ny = NMISS_out, MAF = MAF.x))

MRMix.smk = MRMix(smk_std$betahat_x_std, smk_std$betahat_y_std, smk_std$sx_std, smk_std$sy_std, profile = TRUE)
str(MRMix.smk)



### CPD:
cpd_std = with(df.cpd.1, standardize(BETA_exp, BETA_out, SE_exp, SE_out, xtype="continuous", ytype = "continuous", nx = NMISS_exp, ny = NMISS_out, MAF = MAF.x))

MRMix.cpd = MRMix(cpd_std$betahat_x_std, cpd_std$betahat_y_std, cpd_std$sx_std, cpd_std$sy_std, profile = TRUE)
str(MRMix.cpd)




