# Code for creating survivorship curves, survivorship summary statistics, log-rank/Mantel-Haenszel tests

# Required packages:
library(tidyverse)
library(reshape2)
library(survival)
library(survminer)
library(car)
library(grid)
library(ggpubr)

# Loading and reformatting data--------------------------------------------------------------------

surv_data <- read_csv("survival.csv")

# Changing the order of the factors
surv_data <- surv_data[with(surv_data, order(generation, maternal_age, strain)),]

# Pulling out the censor column, which we will need later
censor <- surv_data$censor
surv_data <- surv_data[,-34]

# Modifying the dataset from wide to long format (needed for plotting and analysis):
surv_long <- melt(data = surv_data, 
                  id.vars = c("generation", "maternal_age", "strain", 
                              "plate", "well"),
                  variable.name = "day",
                  value.name = "dead_alive")

# Survival analysis---------------------------------------------------------------------------------

# Calculate survival times for each individual
survival_data_per_ind <- surv_long %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(surv_time = sum(dead_alive, na.rm = TRUE))

# Adding the censor column back and renaming it
survival_data_per_ind <- cbind(survival_data_per_ind, censor)

survival_data_per_ind <- rename(survival_data_per_ind, "censor" = "...7")

# Making a column with combined full sample information: generation, maternal age, and strain (useful for analysis and plotting later)
survival_data_per_ind$final_line <- paste(survival_data_per_ind$generation,
                                          survival_data_per_ind$maternal_age,
                                          survival_data_per_ind$strain)

# Creating a dataset that does not include censored individuals (for summary statistics)
surv_data_no_censor <-  filter(survival_data_per_ind, censor != 0)

# Calculating summary statistics
surv_summary <- surv_data_no_censor %>%
  group_by(generation, maternal_age, strain) %>%
  summarize(
    mean=mean(surv_time, na.rm = TRUE), 
    SE=sd(surv_time, na.rm =TRUE)/sqrt(length(!is.na(surv_time))),
    n=length(surv_time))
surv_summary

# Conducting log-rank/Mantel-Haenszel tests on the dataset that includes censored individuals
# Comparing Y and 0 and Y and M cohorts within each strain

# BmanL5 strain: -------------------------------------------------------
surv_L5 <-  filter(survival_data_per_ind, strain == "L5")

# Y vs. M cohort:
surv_L5_M <- filter(surv_L5, maternal_age != "old") # excluding 0 cohort

L5_M <- survdiff(Surv(surv_time, censor) ~ maternal_age, data = surv_L5_M)
1 - pchisq(L5_M$chisq, length(L5_M$n) - 1) # calculating p-value using the chi-squared distribution

# Y vs. O cohort:
surv_L5_O <- filter(surv_L5, maternal_age != "mid") # excluding M cohort

L5_O <- survdiff(Surv(surv_time, censor) ~ maternal_age, data = surv_L5_O)
1 - pchisq(L5_O$chisq, length(L5_O$n) - 1)

# BmanRUS strain: -------------------------------------------------------
surv_RUS <-  filter(survival_data_per_ind, strain == "RUS")

# There was no M cohort for this strain
# Y vs. O cohort:
surv_RUS_O <- filter(surv_RUS, generation != "F0") # excluding F0 generation data

RUS_O <- survdiff(Surv(surv_time, censor) ~ maternal_age, data = surv_RUS_O)
1 - pchisq(RUS_O$chisq, length(RUS_O$n) - 1)

# BmanRUS-RE strain: -------------------------------------------------------
surv_RUS_RE <- filter(survival_data_per_ind, strain == "RUS-RE")

# Y vs. M cohort:
surv_RUS_RE_M <- filter(surv_RUS_RE, maternal_age != "old") # excluding O cohort

RUS_RE_M <- survdiff(Surv(surv_time, censor) ~ maternal_age, data = surv_RUS_RE_M)
1 - pchisq(RUS_RE_M$chisq, length(RUS_RE_M$n) - 1)

# Y vs. O cohort:
surv_RUS_RE_O <- filter(surv_RUS_RE, maternal_age != "mid") # excluding M cohort

RUS_RE_O <- survdiff(Surv(surv_time, censor) ~ maternal_age, data = surv_RUS_RE_O)
1 - pchisq(RUS_RE_O$chisq, length(RUS_RE_O$n) - 1)

# BpL1 strain: -------------------------------------------------------
surv_BpL1 <- filter(survival_data_per_ind, strain == "BpL1")

# Y vs. M cohort:
surv_Bp_M <- filter(surv_BpL1, maternal_age != "old") # excluding O cohort

Bp_M <- survdiff(Surv(surv_time, censor) ~ maternal_age, data = surv_Bp_M)
1 - pchisq(Bp_M$chisq, length(Bp_M$n) - 1)

# Y vs. O cohort:
surv_Bp_O <- filter(surv_BpL1, maternal_age != "mid") # excluding M cohort

Bp_O <- survdiff(Surv(surv_time, censor) ~ maternal_age, data = surv_Bp_O)
1 - pchisq(Bp_O$chisq, length(Bp_O$n) - 1)

# Plotting survivorship curves----------------------------------------------------------------------
# The following code should fully reproduce the manuscript figure (multipanel plot of survivorship curves)

# BmanL5 strain:-----------------------------------------------------

surv_L5$final_line <- factor(surv_L5$final_line, 
                             levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"))

fit_L5 <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_L5)

print(fit_L5)
summary(fit_L5)

surv_L5 <- ggsurvplot(
  fit_L5,
  data = surv_L5,
  linetype = c("dotted","solid","solid", "solid"),
  palette = c("black", "black", "#999999", "#0072B2"),
  size = 0.75,
  conf.int = FALSE,        
  pval = FALSE,             
  ggtheme = theme_bw(),      
  xlab = element_blank(),
  xlim = c(0,25),
  ylab = "",
  font.tickslab = c(14))

surv_L5$plot <- surv_L5$plot + theme(plot.margin = unit(c(0,0.1,0,0.5), "lines"))

surv_L5$plot <- surv_L5$plot +  geom_label(
  label="BmanL5", 
  size = 6,
  x=22,
  y=0.9,
  label.size = 0,
  color = "black")
surv_L5

# BmanRUS strain:-----------------------------------------------------

surv_RUS$final_line <- factor(surv_RUS$final_line, 
                             levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"))

fit_RUS <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_RUS)

print(fit_RUS)
summary(fit_RUS)

surv_RUS <- ggsurvplot(
  fit_RUS,
  data = surv_RUS,
  linetype = c("dotted","solid", "solid"),
  palette = c("black", "black", "#0072B2"),
  size = 0.75,        
  conf.int = FALSE,       
  pval = FALSE,           
  ggtheme = theme_bw(), 
  xlab = element_blank(),
  xlim = c(0,25),
  ylab = "",
  font.tickslab = c(14))

surv_RUS$plot <- surv_RUS$plot + theme(plot.margin = unit(c(0,0.1,0,0.5), "lines"))

surv_RUS$plot <- surv_RUS$plot +  geom_label(
  label="BmanRUS", 
  size = 6,
  x=21,
  y=0.9,
  label.size = 0,
  color = "black")
surv_RUS

# BmanRUS-RE strain:-----------------------------------------------------

surv_RUS_RE$final_line <- factor(surv_RUS_RE$final_line, 
                             levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"))

fit_RUS_RE <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_RUS_RE)

print(fit_RUS_RE)
summary(fit_RUS_RE)

surv_RUS_RE <- ggsurvplot(
  fit_RUS_RE,
  data = surv_RUS_RE,
  linetype = c("dotted","solid","solid", "solid"),
  palette = c("black", "black","#999999", "#0072B2"),
  size = 0.75,     
  conf.int = FALSE, 
  pval = FALSE,        
  ggtheme = theme_bw(),   
  xlab = element_blank(),
  xlim = c(0,25),
  ylab = "",
  font.tickslab = c(14))

surv_RUS_RE$plot <- surv_RUS_RE$plot + theme(plot.margin = unit(c(0,0.1,0,0.5), "lines"))

surv_RUS_RE$plot <- surv_RUS_RE$plot +  geom_label(
  label="BmanRUS-RE", 
  size = 6,
  x=20,
  y=0.9,
  label.size = 0,
  color = "black")
surv_RUS_RE

# BpL1 strain:-----------------------------------------------------

surv_BpL1$final_line <- factor(surv_BpL1$final_line, 
                             levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"))

fit_BpL1 <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_BpL1)

print(fit_BpL1)
summary(fit_BpL1)

surv_BpL1 <- ggsurvplot(
  fit_BpL1,
  data = surv_BpL1,
  linetype = c("dotted","solid","solid", "solid"),
  palette = c("black", "black","#999999", "#0072B2"),
  size = 0.75,             
  conf.int = FALSE,     
  pval = FALSE,         
  ggtheme = theme_bw(),
  xlab = element_blank(),
  xlim = c(0,25),
  ylab = "",
  legend.title = "",
  legend.labs = c("F0", "Y", "M", "O"),
  font.tickslab = c(14))

surv_BpL1$plot <- surv_BpL1$plot + theme(legend.key.width = unit(3, "line"),
                                         legend.direction = "horizontal",
                                         plot.margin = unit(c(0,0.1,0,0.5), "lines"),
                                         legend.text=element_text(size=14))

surv_BpL1$plot <- surv_BpL1$plot +  geom_label(
  label="BpL1", 
  size = 6,
  x=23,
  y=0.9,
  label.size = 0,
  color = "black")
surv_BpL1

# Creating the multipanel plot:------------------------------

tiff("surv_curves.tiff", units="in", width=5, height=11, res=300)

surv_curves <- ggarrange(surv_BpL1$plot, surv_RUS$plot, surv_RUS_RE$plot, surv_L5$plot, 
                         ncol = 1, nrow = 4,
                         font.label = list(size = 18),
                         common.legend = TRUE, legend = "top")

annotate_figure(surv_curves, left = textGrob("Survival probability", rot = 90, gp = gpar(fontsize = 20)),
                bottom = textGrob("Age (d)", gp = gpar(fontsize = 20)))

dev.off()

