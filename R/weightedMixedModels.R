
library(lme4)
library(car)
library(RCurl)

data=read.csv(paste0(object$home_dir, "for Mplus/TSST_allTxseq_LONG3.csv"))
data[data==-9999]=NA

weights=data$CTSETA1_lnCORTR

temp=data%>%dplyr::select(s_id, WAVE, contains(c("CTSETA1", "CORT")))

test=lmer(lnCORTR ~ WAVE + I(WAVE^2) +(1| s_id), data=temp, weights=CTSETA1_lnCORTR)

summary(test)
