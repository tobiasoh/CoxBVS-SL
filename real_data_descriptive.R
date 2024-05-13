library(survival)
library(ggplot2)
library(survminer)


library(gtsummary)
load("./RealData/dataset10p_basal.RData")
load("./RealData/brca_full_meta.RData")

# Kaplan-Meier (KM) all data
meta$status = as.integer(meta$vital_status == "Dead")

alive_subjects = which(meta$vital_status=="Alive")
dead_subjects = which(meta$vital_status == "Dead")
meta$time = rep(0,nrow(meta))
meta[alive_subjects,]$time = meta[alive_subjects,]$days_to_last_follow_up
meta[dead_subjects,]$time = meta[dead_subjects,]$days_to_death

meta = meta[-which(is.na(meta$time)),]

sfit <- survival::survfit(Surv(time/365.25, status) ~ 1, data = meta)

# calculate survival probability at 1-, 3- and 5-year time points

#theme_set(theme_bw())
theme_set(theme_classic())
ggsurv <- survminer::ggsurvplot(sfit,
                                conf.int = T, risk.table = TRUE,
                                xlab = "Time since diagnosis (years)",
                                legend = "none", surv.median.line = "hv",
                                break.x.by = 365*3/365,
                                align = "v",
                                legend.title="",
                                legend.labs=c("")
)
ggsurv$plot <- ggsurv$plot +  annotate("text", x = 20, y = 0.9, label = "+  Censor") +
  annotate("text", x=12, y=1, label="A)", fontface = "bold", size=10)

ggsurv



#KM strat based on treatment


sfit <- survival::survfit(Surv(time/365.25, status) ~ as.integer(treatments), data = meta)

# calculate survival probability at 1-, 3- and 5-year time points

theme_set(theme_bw())
ggsurv_treat <- survminer::ggsurvplot(sfit,
                                conf.int = T, risk.table = TRUE,
                                xlab = "Time since diagnosis (years)",
                                legend.title = "Treatment",
                                #legend = "none",#, surv.median.line = "hv"
                                legend.labs = c("No", "Yes"),
                                break.x.by = 365*3/365,
                                risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
                                align = "v"
)
ggsurv_treat$plot <- ggsurv_treat$plot + annotate("text", x = 10, y = 1, label = "+  Censor") +
  annotate("text", x = 20, y =1.2, label = paste0("Log-rank test:\n", surv_pvalue(sfit)$pval.txt))

                            

ggsurv_treat

library(ggpubr)
library(cowplot)
ggarrange(ggsurv$plot, ggsurv_treat$plot, 
          labels = c("a)", "b)"),
          ncol = 1, nrow = 2)

plot_grid(ggsurv$plot, ggsurv_treat$plot, 
          labels = c("a)", "b)"),
          ncol = 1, nrow = 2)


lst = list()
lst[[1]] = ggsurv
lst[[2]] = ggsurv_treat
arrange_ggsurvplots(lst, ncol=1, nrow=2,
                    risk.table.height = 0.4)



#KM cancer subtype specific
meta$Subtype = meta$paper_BRCA_Subtype_PAM50
sfit_subtype = survival::survfit(Surv(time/365.25, status) ~ Subtype, data = meta)

ggsurv_subtype <- survminer::ggsurvplot(sfit_subtype,
                                conf.int = F, risk.table = TRUE,
                                xlab = "Time since diagnosis (years)",
                                legend = "none",#, surv.median.line = "hv"
                                break.x.by = 365*3/365,
                                legend.title="Subtype",
                                legend.labs=c("Basal", "Her2", "Luminal A", "Luminal B", "Normal")
)
ggsurv_subtype$plot <- ggsurv_subtype$plot + annotate("text", x = 21, y = 1, label = "+  Censor") + 
  annotate("text", x = 12, y = 1, label = "B)", fontface="bold", size=10)# +
  #annotate("text", x = 21.5, y =.88, label = paste0("Log-rank test:\n", surv_pvalue(sfit_basal)$pval.txt))


ggsurv_subtype




#
meta_basal = meta[which(meta$paper_BRCA_Subtype_PAM50 == "Basal"),]

sfit_basal = survival::survfit(Surv(time/365.25, status) ~ as.integer(treatments), data=meta_basal)

ggsurv <- survminer::ggsurvplot(sfit_basal,
                                conf.int = T, risk.table = TRUE,
                                xlab = "Time since diagnosis (years)",
                                legend = "none",#, surv.median.line = "hv"
                                break.x.by = 365*3/365
)

ggsurv


#censoring prop total:
sum(meta$status == 0)/nrow(meta)

#censoring prop basal
sum(meta_basal$status == 0)/nrow(meta_basal)

#survival prob 1, 3, 5 years
summ = summary(sfit, times = c(1, 3, 5,10))
summ_df = data.frame(time=summ$time, n.risk = summ$n.risk, 
                     survival=summ$surv, std.err = summ$std.err, 
                     "lower 95% CI"=summ$lower,
                     "Upper 95% CI"=summ$upper)

knitr::kable(summ_df, format="latex")#, digits=3)

#survival probs subtype specific
summary(sfit_subtype, times=c(1,3,5, 10))

#median survival time
tbl_survfit(sfit, probs=0.5, 
            label_header = "**Median survival (95% CI)**"
)

tbl_survfit(sfit_subtype, probs=0.5, 
            label_header = "**Median survival (95% CI)**"
)
#overall and subtype-specific together
tbl_median_surv = tbl_survfit(list(sfit, sfit_subtype), probs=0.5, 
            label_header = "**Median survival (95% CI)**"
)

knitr::kable(tbl_median_surv, format="latex")

#unsupervised learning?
M3C::pca(t(X), labels = meta$Subtype, dotsize = 3, legendtitle = "Subtype")
M3C::tsne(t(X), labels = meta$Subtype, dotsize = 3, legendtitle = "Subtype")
M3C::umap(t(X), labels = meta$Subtype, dotsize = 3, legendtitle = "Subtype")



load("./RealData/result_10p_basal_full.RData")




#check PH assumption for age and treatment:

cox_clin = coxph(Surv(meta$time, meta$status) ~ treatments + age_at_diagnosis, data=meta)
test.ph = cox.zph(cox_clin)
