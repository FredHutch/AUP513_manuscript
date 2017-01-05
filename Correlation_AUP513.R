
library(devtools)
library(ggplot2)
library(grid)
library(gtools)
library(reshape2)
library(plyr)
devtools::install_github("VIDD-VISC/AUP513", host="github.fhcrc.org/api/v3")
rm(list=ls(all=TRUE))
# Load data and R packages
library(AUP513)
data(AUP513_ics)
data(AUP513_ics_marginal)
data(AUP513_ics_marginal_pooled_across_antigen)
data(AUP513_ics_pooled_across_cytokine)
data(AUP513_ics_pooled_across_cytokine_antigen)

data(AUP513_bama)
data(AUP513_adcc_gtl)
data(AUP513_nab)
#data(AUP513_elispot)

# # ICS data
# datics <- AUP513_ics
# table(datics$week)
# datics <- datics[datics$week=="WK26",]
# datics4 <- datics[datics$tcellsub=="CD4",]
# datics4 <- datics4[datics4$antigen=="ENV1",]

# ICS marginal data
daticsm <- AUP513_ics_marginal
table(daticsm$week)
# daticsm <- daticsm[daticsm$week=="WK26" & daticsm$tcellsub=="CD4" & daticsm$antigen=="ENV1" & daticsm$cytokine=="IFNg",]
daticsm <- daticsm[daticsm$week=="WK26",]
daticsm <- daticsm[,c("ptid","group","week","antigen","tcellsub","cytokine","pctpos_adj")]
colnames(daticsm) <- interaction("ics", colnames(daticsm))

# BAMA data
datbama <- AUP513_bama
table(datbama$visitno)
datbama <- datbama[datbama$visitno=="26" & datbama$isotype=="IgG" & datbama$antigen=="gp120 TV1",]
datbama <- datbama[, c("isotype","ptid","visitno","antigen","dilution","auc_calculated")]
colnames(datbama) <- interaction("bama", colnames(datbama))

# GTL data
datgtl <- AUP513_adcc_gtl
table(datgtl$week)
datgtl <- datgtl[datgtl$week=="26" & datgtl$antigen=="TV-1",]
datgtl$log10titer <- log10(datgtl$titer)
datgtl <- datgtl[!is.na(datgtl$log10titer),]
datgtl <- datgtl[,c("antigen","ptid","group","week","log10titer")]
colnames(datgtl) <- interaction("gtl", colnames(datgtl))

# NAB data
datnab <- AUP513_nab
table(datnab$visitno)
datnab <- datnab[datnab$visitno=="26" & datnab$celltype=="TZM-bl" & datnab$isolate=="MW965.26",]
datnab <- datnab[,c("celltype","ptid","isolate","visitno","titer_mod")]
colnames(datnab) <- interaction("nab", colnames(datnab))

# Merge all
dat <- merge(daticsm, datbama, by.x="ics.ptid", by.y="bama.ptid", all=T)
dat <- merge(dat, datgtl, by.x="ics.ptid", by.y="gtl.ptid", all=T)
dat <- merge(dat, datnab, by.x="ics.ptid", by.y="nab.ptid", all=T)
dat <- dat[!is.na(dat$ics.tcellsub) & !is.na(dat$bama.auc_calculated) & !is.na(dat$gtl.log10titer) & !is.na(dat$nab.titer_mod),]

datTot <- dat
#dat <- datTot[datTot$ics.group==6 & datTot$ics.tcellsub=="CD4" & datTot$ics.cytokine=="IFNg",]
#Make Correlation Scatter Plot
#dat <- dat[,c("ics.ptid","ics.group","ics.week","ics.antigen","ics.tcellsub", "ics.cytokine",  "ics.pctpos_adj","bama.auc_calculated","gtl.log10titer","nab.titer_mod")]
#datlong <- melt(dat, id.vars=c("ics.ptid", "ics.group", "ics.week", "ics.antigen", "ics.tcellsub", "ics.cytokine", "ics.pctpos_adj"))
# # Only CD4
# dat <- datTot[datTot$ics.group==6 & datTot$ics.tcellsub=="CD4",]
# cors <- ddply(dat, c("ics.antigen","ics.cytokine"), summarise, cor = round(cor(ics.pctpos_adj, bama.auc_calculated, method="spearman"), 2))
# ggplot(dat, aes(x=bama.auc_calculated, y=ics.pctpos_adj)) + theme_bw() +
#   xlab("bama.auc_calculated") +
#   ylab("ics.pctpos_adj") +
#   facet_grid(ics.cytokine ~ ics.antigen, scales="free") + 
#   geom_point() +    
#   #scale_colour_manual(name="Group", values=c("CHRONIC Treated"="red", "chronic untreated"="blue", "elite controller"="green","VAX004"="orange", "Placebo"="black")) +
#   geom_text(data=cors, size=3, aes(label=paste("r=", cor, sep="")), x=48000, y=0.4) + 
#   #geom_text(data=cors, size=3, aes(label=paste("r=", cor, sep="")), x=48000, y=0.5, color="red") + 
#   geom_smooth(method=lm, se=FALSE)  +  # Don't add shaded confidence region
#   theme(#plot.margin = unit(c(1,1,0,1), "cm"),
#     legend.position="top",
#     legend.title=element_text(size=7,vjust=-1.6),
#     legend.text=element_text(size=7,vjust=-1.6))

# CD4 and CD8
dat <- datTot[datTot$ics.group==6,]
dat$response <- dat$bama.auc_calculated
corsCD4 <- ddply(dat[dat$ics.tcellsub=="CD4",], c("ics.antigen","ics.cytokine","ics.tcellsub"), summarise, cor = round(cor(ics.pctpos_adj, response, method="spearman"), 2))
corsCD8 <- ddply(dat[dat$ics.tcellsub=="CD8",], c("ics.antigen","ics.cytokine","ics.tcellsub"), summarise, cor = round(cor(ics.pctpos_adj, response, method="spearman"), 2))

ggplot(dat, aes(x=response, y=ics.pctpos_adj, colour=ics.tcellsub)) + theme_bw() +
  xlab("bama.auc_calculated") +
  ylab("ics.pctpos_adj") +
  facet_grid(ics.cytokine ~ ics.antigen) + 
  geom_point() +    
  scale_colour_manual(name="ics.tcellsub", values=c("CD4"="red", "CD8"="blue")) +
  geom_text(data=corsCD4, size=4, aes(label=paste("r =", cor, sep=" ")), x=52000, y=0.7, color="red") + 
  geom_text(data=corsCD8, size=4, aes(label=paste("r =", cor, sep=" ")), x=52000, y=0.65, color="blue") + 
  geom_smooth(method=lm, se=FALSE)  +  # Don't add shaded confidence region
  theme(#plot.margin = unit(c(1,1,0,1), "cm"),
    legend.position="top",
    legend.title=element_text(size=7,vjust=-1.6),
    legend.text=element_text(size=7,vjust=-1.6))










# 
# 
# ##############################
# dat <- dat[,c("ics.ptid","ics.group","ics.week","ics.antigen","ics.tcellsub", "ics.cytokine",  "ics.pctpos_adj","bama.auc_calculated","gtl.log10titer","nab.titer_mod")]
# 
# datlong <- melt(dat, id.vars=c("ics.ptid", "ics.group", "ics.week", "ics.antigen", "ics.tcellsub", "ics.cytokine", "ics.pctpos_adj"))
# cors <- ddply(datlong, c("ics.antigen","variable"), summarise, cor = round(cor(ics.pctpos_adj, value, method="spearman"), 2))
# 
# ggplot(datlong, aes(x=value, y=ics.pctpos_adj)) + theme_bw() +
#   xlab("bama.auc/gtl.log10titer/nab.titer") +
#   ylab("ics.pctpos_adj") +
#   facet_grid(ics.antigen ~ variable, scales="free") + 
#   geom_point() +    
#   #scale_colour_manual(name="Group", values=c("CHRONIC Treated"="red", "chronic untreated"="blue", "elite controller"="green","VAX004"="orange", "Placebo"="black")) +
#   #geom_text(data=cors, size=3, aes(label=paste("r=", cor, sep=""))) + 
#   geom_smooth(method=lm, se=FALSE)  +  # Don't add shaded confidence region
#   theme(#plot.margin = unit(c(1,1,0,1), "cm"),
#     legend.position="top",
#     legend.title=element_text(size=7,vjust=-1.6),
#     legend.text=element_text(size=7,vjust=-1.6))
# 
# ###########################
# 
# ggplot(gtl.long,aes(x = factor(group),y = measurement)) + 
#   facet_grid(condition~week,scales="free") +
#   theme_bw() + xlab("Group")  + 
#   ylab("     AUC                                                Peak                                             Log10_titer       ") +
#   scale_colour_manual(name="Response", values = c("0"="blue","1"="red"),labels=c("Non-Responders","Responders")) +
#   geom_point(size=2,position = position_jitter(width = 0.2, height=0),aes(colour = factor(response)))  +
#   geom_boxplot(data=gtl.long[gtl.long$response!=0,],outlier.colour = "NA",alpha=0) +
#   geom_hline(aes(yintercept = z), col="blue", lty=3) +
#   theme(legend.position="top",
#         legend.title=element_text(size=7,vjust=-1.6), 
#         legend.text=element_text(size=7,vjust=-1.6), 
#         axis.text.x = element_text(size=8),
#         axis.text.y = element_text(size=8),
#         axis.title.x = element_text(size=10,vjust=-0.5), 
#         axis.title.y = element_text(angle=90, size=10,vjust=0.5),
#         strip.text.x = element_text(size = 9, face="bold"), strip.text.y = element_blank(),strip.background = element_blank())
# 




