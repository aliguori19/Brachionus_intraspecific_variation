# Required packages:
library(tidyverse)
library(reshape2)
library(survival)
library(survminer)
library(car)
library(grid)
library(ggpubr)
library(flexsurv)
library(MASS)
library(emmeans)

# Lifespan--------------------------------------------------------------------

# Loading and reformatting data
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

# Creating a dataset that does not include censored individuals
surv_data_no_censor <-  filter(survival_data_per_ind, censor != 0)

# Changing order of maternal age levels
surv_data_no_censor$maternal_age <- factor(surv_data_no_censor$maternal_age, levels = c("young","mid","old"))

# Calculating summary statistics
surv_summary <- surv_data_no_censor %>%
  group_by(generation, maternal_age, strain) %>%
  summarize(
    mean=mean(surv_time, na.rm = TRUE), 
    SE=sd(surv_time, na.rm =TRUE)/sqrt(length(!is.na(surv_time))),
    n=length(surv_time))
surv_summary

# GLM on lifespan (count data)

surv_data_no_censor_F1 <- filter(surv_data_no_censor, generation == "F1")

# Checking data for assumptions of normality and homogeneity of variance

surv_aov <- aov((surv_time) ~ maternal_age * strain, surv_data_no_censor_F1)
surv_resid <- residuals(surv_aov)

shapiro.test(surv_resid)
leveneTest((surv_time) ~ maternal_age * strain, surv_data_no_censor_F1)

# Assumptions are not met, can't find any useful transformations (count data)
# Will use GLM instead
# Using Poisson for count data

pois.glm_surv <- glm(surv_time ~ maternal_age * strain, family = poisson, data = surv_data_no_censor_F1)
summary(pois.glm_surv)

# (null deviance - residual deviance)/null deviance
(1291.5-850.16)/1291.5

drop1(pois.glm_surv, test = "Chi") # significant interaction

# Model validation (plotting residuals)
plot(pois.glm_surv)

# Pairwise comparisons: which maternal age groups differ within strains?

EMM <- emmeans(pois.glm_surv, ~ maternal_age * strain)
EMM    # displays the cell means

# Simple pairwise comparisons
pairs(EMM, simple = "maternal_age") 

# Plotting survivorship curves

# BmanL5

surv_L5_nc <- filter(surv_data_no_censor, strain == "L5")

surv_L5_nc$final_line <- factor(surv_L5_nc$final_line, 
                                levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"))

fit_L5 <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_L5_nc)

print(fit_L5)
summary(fit_L5)

surv_L5 <- ggsurvplot(
  fit_L5,
  data = surv_L5_nc,
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

# BmanRUS

surv_RUS_nc <- filter(surv_data_no_censor, strain == "RUS")

surv_RUS_nc$final_line <- factor(surv_RUS_nc$final_line, 
                                 levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"))

fit_RUS <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_RUS_nc)

print(fit_RUS)
summary(fit_RUS)

surv_RUS <- ggsurvplot(
  fit_RUS,
  data = surv_RUS_nc,
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

# BmanRUS-RE

surv_RUS_RE_nc <- filter(surv_data_no_censor, strain == "RUS-RE")

surv_RUS_RE_nc$final_line <- factor(surv_RUS_RE_nc$final_line, 
                                    levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"))

fit_RUS_RE <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_RUS_RE_nc)

print(fit_RUS_RE)
summary(fit_RUS_RE)

surv_RUS_RE <- ggsurvplot(
  fit_RUS_RE,
  data = surv_RUS_RE_nc,
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

# BpL1

surv_BpL1_nc <- filter(surv_data_no_censor, strain == "BpL1")

surv_BpL1_nc$final_line <- factor(surv_BpL1_nc$final_line, 
                                  levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"))

fit_BpL1 <- survfit(Surv(surv_time, censor) ~ final_line, data = surv_BpL1_nc)

print(fit_BpL1)
summary(fit_BpL1)

surv_BpL1 <- ggsurvplot(
  fit_BpL1,
  data = surv_BpL1_nc,
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

# Creating a multipanel plot

tiff("surv_curves_nc.tiff", units="in", width=5, height=11, res=300)

surv_curves <- ggarrange(surv_BpL1$plot, surv_RUS$plot, surv_RUS_RE$plot, surv_L5$plot, 
                         ncol = 1, nrow = 4,
                         font.label = list(size = 18),
                         common.legend = TRUE, legend = "top")

annotate_figure(surv_curves, left = textGrob("Survival probability", rot = 90, gp = gpar(fontsize = 20)),
                bottom = textGrob("Age (d)", gp = gpar(fontsize = 20)))

dev.off()

# LRO--------------------------

# Loading and reformatting data

neonate_data <- read_csv("neonates.csv")

# Changing the order of the factors
neonate_data <- neonate_data[with(neonate_data, order(generation, maternal_age, strain)),]

# Removing censored individuals (we were not able to quantify their LRO, lost before natural death)
neonate_data_no_censor <- filter(neonate_data, censor != 0)

# Modifying the dataset from wide to long format (needed for plotting and analysis):
neo_long_nc <- melt(data = neonate_data_no_censor,
                    id.vars = c("generation", "maternal_age", "strain", "plate", "well"),
                    variable.name = "day",
                    value.name = "num_neonates")

# Removing 'censor' altogether, not needed after censored individuals are removed
neo_long_nc <- filter(neo_long_nc, day != "censor")

# Changing order of maternal age levels
neo_long_nc$maternal_age <- factor(neo_long_nc$maternal_age, levels = c("young","mid","old"))

# Calculating and analyzing lifetime reproductive output (LRO) per individual

# Calculating LRO, maximum LRO, and 50% LRO (needed for downstream analyses)
LRO_per_ind <- neo_long_nc %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(LRO = sum(num_neonates, na.rm = TRUE),
            max_RO = max(num_neonates, na.rm = TRUE),
            fifty_percent_LRO = (sum(num_neonates, na.rm = TRUE))*0.5)

LRO_per_ind$final_line <- paste(LRO_per_ind$generation,
                                LRO_per_ind$maternal_age,
                                LRO_per_ind$strain)

LRO_per_ind$maternal_age <- factor(LRO_per_ind$maternal_age)

# Summary statistics for LRO and maximum LRO
LRO_summary <- LRO_per_ind %>%
  group_by(generation, maternal_age, strain) %>%
  summarize(
    mean=mean(LRO, na.rm = TRUE), 
    med=median(LRO, na.rm = TRUE),
    SE=sd(LRO, na.rm =TRUE)/sqrt(length(!is.na(LRO))),
    n=length(LRO),
    mean_max=mean(max_RO, na.rm=TRUE),
    SE_max=sd(max_RO, na.rm =TRUE)/sqrt(length(!is.na(max_RO))),
    med_max=median(max_RO, na.rm = TRUE),
    min=min(LRO, na.rm=TRUE),
    max=max(LRO, na.rm=TRUE),
    var=var(LRO, na.rm = TRUE))
LRO_summary

LRO_summary$final_line <- paste(LRO_summary$generation, # needed for plotting later
                                LRO_summary$maternal_age,
                                LRO_summary$strain)

LRO_summary$maternal_age <- factor(LRO_summary$maternal_age)

# Statistical analyses for LRO (GLM)

# Excluding the F0 generation

LRO_per_ind_F1 <- filter(LRO_per_ind, generation == "F1")

# Checking data for assumptions of normality and homogeneity of variance

LRO_aov <- aov((LRO) ~ maternal_age * strain, LRO_per_ind_F1)
LRO_resid <- residuals(LRO_aov)

shapiro.test(LRO_resid)
leveneTest((LRO) ~ maternal_age * strain, LRO_per_ind_F1)

# Assumptions are not met, can't find any useful transformations (count data)
# Will use GLM
# LRO data are also overdispersed- using negative binomial GLM

nb.glm <- glm.nb(LRO ~ maternal_age * strain, link = "log", data = LRO_per_ind_F1)
summary(nb.glm, cor = FALSE)

# (null deviance - residual deviance)/null deviance
(2247.16-988.16)/2247.16

drop1(nb.glm, test = "Chi")

# Model validation- (plotting residuals)
plot(nb.glm)

# Pairwise comparisons: which maternal age groups differ within strains?

EMM <- emmeans(nb.glm, ~ maternal_age * strain)
EMM    # displays the cell means

# Simple pairwise comparisons
pairs(EMM, simple = "maternal_age")    # compare maternal age cohorts for each strain

# LRO Figures

# BmanL5:
wilcox_L5 <- c("","a","b","b") # from pairwise analyses

LRO_sum_L5 <- filter(LRO_summary, strain == "L5")

LRO_sum_L5$final_line <- factor(LRO_sum_L5$final_line, levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"), labels = c("F0", "Y", "M","O"))

LRO_per_ind_L5 <- filter(LRO_per_ind, strain == "L5")

LRO_per_ind_L5$final_line <- factor(LRO_per_ind_L5$final_line, levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"), labels = c("F0", "Y", "M","O"))

# Plotting means and SE
LRO_raw <- ggplot(data = LRO_sum_L5, mapping = aes(x = final_line, y = mean))+
  geom_point(size = 3) +
  geom_errorbar(data = LRO_sum_L5, 
                aes(x = final_line, ymin=mean-SE, ymax=mean+SE), width=0.1) +
  geom_text(label=wilcox_L5, 
            nudge_x = 0.3,
            size = 8)

# Plotting the raw data points
LRO_combo <- LRO_raw +
  geom_jitter(data = LRO_per_ind_L5, aes(x = final_line, y = LRO),
              width=.2, alpha=0.2, size = 2)

LRO_L5_plot <- LRO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(-0.5,38)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        #axis.title=element_text(size=24), 
        axis.text=element_text(size=14))
LRO_L5_plot

# BmanRUS:
LRO_sum_RUS <- filter(LRO_summary, strain == "RUS")

LRO_sum_RUS$final_line <- factor(LRO_sum_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"), labels = c("F0", "Y","O"))

LRO_per_ind_RUS <- filter(LRO_per_ind, strain == "RUS")

LRO_per_ind_RUS$final_line <- factor(LRO_per_ind_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"), labels = c("F0", "Y","O"))

# Plotting means and SE
LRO_raw <- ggplot(data = LRO_sum_RUS, mapping = aes(x = final_line, y = mean))+
  geom_point(size = 3) +
  geom_errorbar(data = LRO_sum_RUS, 
                aes(x = final_line, ymin=mean-SE, ymax=mean+SE), width=0.1)

# Plotting the raw data points
LRO_combo <- LRO_raw +
  geom_jitter(data = LRO_per_ind_RUS, aes(x = final_line, y = LRO),
              width=.2, alpha=0.2, size = 2)

LRO_RUS_plot <- LRO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(-0.5,38)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=14))
LRO_RUS_plot

# BmanRUS-RE:
wilcox_RUS_RE <- c("","a","ab","b") # from pairwise analyses

LRO_sum_RUS_RE <- filter(LRO_summary, strain == "RUS-RE")

LRO_sum_RUS_RE$final_line <- factor(LRO_sum_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"), labels = c("F0", "Y", "M","O"))

LRO_per_ind_RUS_RE <- filter(LRO_per_ind, strain == "RUS-RE")

LRO_per_ind_RUS_RE$final_line <- factor(LRO_per_ind_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"), labels = c("F0", "Y", "M","O"))

# Plotting means and SE
LRO_raw <- ggplot(data = LRO_sum_RUS_RE, mapping = aes(x = final_line, y = mean))+
  geom_point(size = 3) +
  geom_errorbar(data = LRO_sum_RUS_RE, 
                aes(x = final_line, ymin=mean-SE, ymax=mean+SE), width=0.1) +
  geom_text(label=wilcox_RUS_RE, 
            nudge_x = 0.3,
            size = 8)

# Plotting the raw data points
LRO_combo <- LRO_raw +
  geom_jitter(data = LRO_per_ind_RUS_RE, aes(x = final_line, y = LRO),
              width=.2, alpha=0.2, size = 2)

LRO_RUS_RE_plot <- LRO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(-0.5,38)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=14))
LRO_RUS_RE_plot

#BpL1:
wilcox_bp <- c("","a","b","b") # from pairwise analyses

LRO_sum_BpL1 <- filter(LRO_summary, strain == "BpL1")

LRO_sum_BpL1$final_line <- factor(LRO_sum_BpL1$final_line, levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"), labels = c("F0", "Y", "M","O"))

LRO_per_ind_BpL1 <- filter(LRO_per_ind, strain == "BpL1")

LRO_per_ind_BpL1$final_line <- factor(LRO_per_ind_BpL1$final_line, levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"), labels = c("F0", "Y", "M","O"))

# Plotting means and SE
LRO_raw <- ggplot(data = LRO_sum_BpL1, mapping = aes(x = final_line, y = mean))+
  geom_point(size = 3) +
  geom_errorbar(data = LRO_sum_BpL1, 
                aes(x = final_line, ymin=mean-SE, ymax=mean+SE), width=0.1) +
  geom_text(label=wilcox_bp, 
            nudge_x = 0.3,
            size = 8)

# Plotting the raw data points
LRO_combo <- LRO_raw +
  geom_jitter(data = LRO_per_ind_BpL1, aes(x = final_line, y = LRO),
              width=.2, alpha=0.2, size = 2)

LRO_Bp_plot <- LRO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(-0.5,38)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=14))
LRO_Bp_plot


# Plotting daily reproduction over time-----------------

# Calculating the mean number of neonates produced across all individuals, per day
RO_per_day <- neo_long_nc %>%
  group_by(generation, maternal_age, strain, day) %>%
  summarize(total_offspring_per_day = sum(num_neonates, na.rm = TRUE),
            mean = mean(num_neonates, na.rm = TRUE),
            SE = sd(num_neonates, na.rm =   
                      TRUE)/sqrt(length(!is.na(num_neonates))))

RO_per_day$final_line <- paste(RO_per_day$generation,
                               RO_per_day$maternal_age,
                               RO_per_day$strain)

# Plotting for each strain:

# BmanL5:
RO_L5 <- filter(RO_per_day, strain == "L5")

RO_L5$final_line <- factor(RO_L5$final_line, levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"))

RO_raw <- ggplot(data=RO_L5, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw+
  geom_line(data=RO_L5, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black", "#999999", "#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid", "solid"),labels = c("F0","Y","M","O")) #+  geom_label(
   # label="BmanL5", 
   # size = 8,
   # x=24,
   # y=4.5,
  #  label.size = 0,
   # color = "black")

RO_L5_plot_time <- RO_combo + 
  ylab("") +
  ylim(c(0,5)) +
  xlab(element_blank()) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = "none")
RO_L5_plot_time

# BmanRUS:
RO_RUS <- filter(RO_per_day, strain == "RUS")

RO_RUS$final_line <- factor(RO_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"))

RO_raw <- ggplot(data=RO_RUS, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw+
  geom_line(data=RO_RUS, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black","#0072B2"), labels = c("F0","Y", "O")) +
  scale_linetype_manual(values = c("dotted","solid","solid"),labels = c("F0","Y", "O")) # +
 # geom_label(
   # label="BmanRUS", 
   # size = 8,
   # x=23,
   # y=4.5,
   # label.size = 0,
   # color = "black")

RO_RUS_plot_time <- RO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(0,5)) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = "none")
RO_RUS_plot_time

# BmanRUS-RE:
RO_RUS_RE <- filter(RO_per_day, strain == "RUS-RE")

RO_RUS_RE$final_line <- factor(RO_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"))

RO_raw <- ggplot(data=RO_RUS_RE, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw+
  geom_line(data=RO_RUS_RE, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black","#999999","#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid","solid"),labels = c("F0","Y","M","O")) #+ geom_label(
  #  label="BmanRUS-RE", 
   # size = 8,
   # x=22,
   # y=4.5,
   # label.size = 0,
   # color = "black")

RO_RUS_RE_time <- RO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(0,5)) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = "none")
RO_RUS_RE_time

# BpL1:
RO_BpL1 <- filter(RO_per_day, strain == "BpL1")

RO_BpL1$final_line <- factor(RO_BpL1$final_line, levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"))

RO_raw <- ggplot(data=RO_BpL1, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw +
  geom_line(data=RO_BpL1, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black","#999999","#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid","solid"),labels = c("F0","Y","M","O"))# + geom_label(
   # label="BpL1", 
   # size = 8,
   # x=25,
   # y=4.5,
   # label.size = 0,
   # color = "black")

RO_Bp_plot_time <- RO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(0,5)) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = c(.85,.5),
        legend.title = element_blank(),
        legend.key.width = unit(3, "line"),
        legend.text=element_text(size=16))
RO_Bp_plot_time

# Maximum daily reproduction (MDR) analyses---------------------

# Checking data for assumptions of normality and homogeneity of variance

MDR_aov <- aov((max_RO) ~ maternal_age * strain, LRO_per_ind_F1)
MDR_resid <- residuals(MDR_aov)

shapiro.test(MDR_resid)
leveneTest((max_RO) ~ maternal_age * strain, LRO_per_ind_F1)

# Assumptions are not met, can't find any useful transformations (count data)
# Will use GLM instead
# Using Poisson for count data

pois.glm_MDR <- glm(max_RO ~ maternal_age * strain, family = poisson, data = LRO_per_ind_F1)
summary(pois.glm_MDR)

# (null deviance - residual deviance)/null deviance
(291.71-137.7)/291.71

drop1(pois.glm_MDR, test = "Chi") # no significant interaction

# Trying model without the interaction term:
pois.glm_MDR2 <- glm(max_RO ~ maternal_age + strain, family = poisson, data = LRO_per_ind_F1)
summary(pois.glm_MDR2) # slightly lower AIC than model with interaction term

# Model validation (plotting residuals)
plot(pois.glm_MDR2)

# Pairwise comparisons: which levels are different from each other within each factor?

EMM1 <- emmeans(pois.glm_MDR2, ~ maternal_age)
pairs(EMM1)

EMM2 <- emmeans(pois.glm_MDR2, ~ strain)
pairs(EMM2)

# Timing of MDR and onset of reproduction-------------------

# Calculating cumulative reproduction PER INDIVIDUAL:

neo_long_nc <- neo_long_nc[with(neo_long_nc, order(generation, maternal_age, strain, plate, well)),]

neo_long_nc <- neo_long_nc %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  mutate(cum_RO = cumsum(num_neonates))

# Removing individuals with a LRO of zero (not relevant here)
LRO_per_ind <- filter(LRO_per_ind, LRO != 0)

neo_merge <- merge(LRO_per_ind, neo_long_nc, by = c("generation","maternal_age","strain","plate","well"))

neo_merge$max_match <- neo_merge$num_neonates == neo_merge$max_RO

MDR_time <- filter(neo_merge, max_match == "TRUE")

MDR_time$day <- as.numeric(MDR_time$day)

earliest_max_day <- MDR_time %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(earliest = min(day))

# Summarizing time of MDR:
MDR_time_summary <- earliest_max_day %>%
  group_by(generation, maternal_age, strain) %>%
  summarize(
    mean=mean(earliest, na.rm = TRUE), 
    med=median(earliest, na.rm = TRUE),
    n=length(earliest),
    SE=sd(earliest, na.rm =TRUE)/sqrt(length(!is.na(earliest))))
MDR_time_summary

# Timing of MDR analyses

# Checking data for assumptions of normality and homogeneity of variance

MDR_time_aov <- aov((earliest) ~ maternal_age * strain, earliest_max_day)
MDR_time_resid <- residuals(MDR_time_aov)

shapiro.test(MDR_time_resid)
leveneTest((earliest) ~ maternal_age * strain, earliest_max_day)

# Assumptions are not met, can't find any useful transformations (count data)
# Will use GLM instead
# Using Poisson for count data

pois.glmMDRtime <- glm(earliest ~ maternal_age * strain, family = poisson, data = earliest_max_day)
summary(pois.glmMDRtime)

drop1(pois.glmMDRtime, test = "Chi") # no significant interaction

# Trying model without the interaction term:
pois.glmMDRtime2 <- glm(earliest ~ maternal_age + strain, family = poisson, data = earliest_max_day)
summary(pois.glmMDRtime2)

# Pairwise comparisons: which levels are different from each other within each factor?

EMMA <- emmeans(pois.glmMDRtime2, ~ maternal_age)
pairs(EMMA)

EMMB <- emmeans(pois.glmMDRtime2, ~ strain)
pairs(EMMB)

# Onset of reproduction calculations

neo_merge$onset_match <- neo_merge$num_neonates == 0

onset <- filter(neo_merge, onset_match == "FALSE")

onset$day <- as.numeric(onset$day)

rep_onset <- onset %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(earliest = min(day))

# Summarizing time of rep onset:
onset_summary <- rep_onset %>%
  group_by(generation, maternal_age, strain) %>%
  summarize(
    mean=mean(earliest, na.rm = TRUE), 
    med=median(earliest, na.rm = TRUE),
    n=length(earliest),
    SE=sd(earliest, na.rm =TRUE)/sqrt(length(!is.na(earliest))))
onset_summary

# Onset of reproduction analyses

# Checking data for assumptions of normality and homogeneity of variance

onset_aov <- aov((earliest) ~ maternal_age * strain, rep_onset)
onset_resid <- residuals(onset_aov)

shapiro.test(onset_resid)
leveneTest((earliest) ~ maternal_age * strain, rep_onset)

# Assumptions are not met, can't find any useful transformations (count data)
# Will use GLM instead
# Using Poisson for count data

pois.glm.onset <- glm(earliest ~ maternal_age * strain, family = poisson, data = rep_onset)
summary(pois.glm.onset)

drop1(pois.glm.onset, test = "Chi") # no significant interaction

# Trying model without the interaction term:
pois.glm.onset2 <- glm(earliest ~ maternal_age + strain, family = poisson, data = rep_onset)
summary(pois.glm.onset)

# Pairwise comparisons: which levels are different from each other within each factor?

EMMA <- emmeans(pois.glm.onset2, ~ maternal_age)
pairs(EMMA)

EMMB <- emmeans(pois.glm.onset2, ~ strain)
pairs(EMMB)

# Age of 50% LRO-------------------------

# Calculating the age of 50% LRO, using cumulative reproduction data

sub <- subset(neo_merge, cum_RO >= fifty_percent_LRO , select = c(generation,maternal_age,strain,plate,well,LRO,fifty_percent_LRO,day,cum_RO))

sub$day <- as.numeric(sub$day)

# Contains the time of 50% LRO production (d) for each individual
sub_summ <- sub %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(time_fifty = min(day))

# Summarizing time of 50% LRO:
LRO_50_summary <- sub_summ %>%
  group_by(generation, maternal_age, strain) %>%
  summarize(
    mean=mean(time_fifty, na.rm = TRUE), 
    med=median(time_fifty, na.rm = TRUE),
    SE=sd(time_fifty, na.rm =TRUE)/sqrt(length(!is.na(time_fifty))))
LRO_50_summary

# Checking data for assumptions of normality and homogeneity of variance

LRO50_aov <- aov((time_fifty) ~ maternal_age * strain, sub_summ)
LRO50_resid <- residuals(LRO50_aov)

shapiro.test(LRO50_resid)
leveneTest((time_fifty) ~ maternal_age * strain, sub_summ)

# Assumptions are not met, can't find any useful transformations (count data)
# Will use GLM instead
# Using Poisson for count data

pois.glm50 <- glm(time_fifty ~ maternal_age * strain, family = poisson, data = sub_summ)
summary(pois.glm50)

# (null deviance - residual deviance)/null deviance
(177.67-122.7)/177.67

drop1(pois.glm50, test = "Chi") # no significant interaction

# Trying model without the interaction term:
pois.glm502 <- glm(time_fifty ~ maternal_age + strain, family = poisson, data = sub_summ)
summary(pois.glm502) # slightly lower AIC than model with interaction term

# Model validation (plotting residuals)
plot(pois.glm502)

# Pairwise comparisons: which levels are different from each other within each factor?

EMMA <- emmeans(pois.glm502, ~ maternal_age)
pairs(EMMA)

EMMB <- emmeans(pois.glm502, ~ strain)
pairs(EMMB)




# Intrinsic value------------

intrinsic_value <- neo_merge %>%
  group_by(generation, maternal_age, strain, plate, well, day) %>%
  summarize(remaining_LRO = LRO-cum_RO,
            int_val = remaining_LRO/LRO)

# Plotting over time

# Calculating the mean number of neonates produced across all individuals, per day
intval_per_day <- intrinsic_value %>%
  group_by(generation, maternal_age, strain, day) %>%
  summarize(mean = mean(int_val, na.rm = TRUE),
            SE = sd(int_val, na.rm =   
                      TRUE)/sqrt(length(!is.na(int_val))))

intval_per_day$final_line <- paste(intval_per_day$generation,
                                   intval_per_day$maternal_age,
                                   intval_per_day$strain)

# Plotting for each strain:
# BmanL5:
intval_L5 <- filter(intval_per_day, strain == "L5")

intval_L5$final_line <- factor(intval_L5$final_line, levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"))

intval_raw <- ggplot(data=intval_L5, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

intval_combo <- intval_raw+
  geom_line(data=intval_L5, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black", "#999999", "#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid", "solid"),labels = c("F0","Y","M","O")) +  geom_label(
    label="BmanL5", 
    size = 8,
    x=24,
    y=0.9,
    label.size = 0,
    color = "black")

intval_L5_plot_time <- intval_combo + 
  ylab("") +
  #ylim(c(0,5)) +
  xlab(element_blank()) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = "none")
intval_L5_plot_time

# BmanRUS:
intval_RUS <- filter(intval_per_day, strain == "RUS")

intval_RUS$final_line <- factor(intval_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"))

intval_raw <- ggplot(data=intval_RUS, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

intval_combo <- intval_raw+
  geom_line(data=intval_RUS, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black", "#0072B2"), labels = c("F0","Y","O")) +
  scale_linetype_manual(values = c("dotted","solid", "solid"),labels = c("F0","Y","O")) +  geom_label(
    label="BmanRUS", 
    size = 8,
    x=23,
    y=0.9,
    label.size = 0,
    color = "black")

intval_RUS_plot_time <- intval_combo + 
  ylab("") +
  #ylim(c(0,5)) +
  xlab(element_blank()) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = "none")
intval_RUS_plot_time

# BmanRUS-RE:
intval_RUS_RE <- filter(intval_per_day, strain == "RUS-RE")

intval_RUS_RE$final_line <- factor(intval_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"))

intval_raw <- ggplot(data=intval_RUS_RE, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

intval_combo <- intval_raw+
  geom_line(data=intval_RUS_RE, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black", "#999999", "#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid", "solid"),labels = c("F0","Y","M","O")) +  geom_label(
    label="BmanRUS-RE", 
    size = 8,
    x=21,
    y=0.9,
    label.size = 0,
    color = "black")

intval_RUS_RE_plot_time <- intval_combo + 
  ylab("") +
  #ylim(c(0,5)) +
  xlab(element_blank()) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = "none")
intval_RUS_RE_plot_time

# BpL1:
intval_Bp <- filter(intval_per_day, strain == "BpL1")

intval_Bp$final_line <- factor(intval_Bp$final_line, levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"))

intval_raw <- ggplot(data=intval_Bp, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

intval_combo <- intval_raw+
  geom_line(data=intval_Bp, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black", "#999999", "#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid", "solid"),labels = c("F0","Y","M","O")) +  geom_label(
    label="BpL1", 
    size = 8,
    x=25,
    y=0.9,
    label.size = 0,
    color = "black")

intval_BpL1_plot_time <- intval_combo + 
  ylab("") +
  #ylim(c(0,5)) +
  xlab(element_blank()) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.text=element_text(size=14),
        legend.position = c(.85,.5),
        legend.title = element_blank(),
        legend.key.width = unit(3, "line"),
        legend.text=element_text(size=16))
intval_BpL1_plot_time


# Multipanel reproduction plot---------------------

tiff("rep.tiff", units="in", width=16, height=12, res=300)

RO_curves <- ggarrange(RO_Bp_plot_time, RO_RUS_plot_time,
                       RO_RUS_RE_time, RO_L5_plot_time, 
                       ncol = 1, nrow = 4,
                       labels = c("B","E","H","K"),
                       font.label = list(size = 18),
                       hjust = -0.25)

RO_curves <- annotate_figure(RO_curves, bottom = textGrob("Age (d)", gp = gpar(fontsize = 24)), 
                             left = textGrob(expression(paste("Offspring ind"^"-1"*"day"^"-1")), 
                               rot = 90, gp = gpar(fontsize = 24)))

LRO_summary <- ggarrange(LRO_Bp_plot, LRO_RUS_plot,
                         LRO_RUS_RE_plot, LRO_L5_plot, 
                         ncol = 1, nrow = 4,
                         labels = c("A","D","G","J"),
                         font.label = list(size = 18),
                         hjust = -0.25)

LRO_summary <- annotate_figure(LRO_summary, bottom = textGrob("Cohort", gp = gpar(fontsize = 24)), 
                               left = textGrob(expression(paste("Lifetime reproductive output (offspring ind"^"-1"*")")),
                                rot = 90, gp = gpar(fontsize = 24)))

int_value_combo <- ggarrange(intval_BpL1_plot_time, intval_RUS_plot_time,
                             intval_RUS_RE_plot_time, intval_L5_plot_time, 
                             ncol = 1, nrow = 4,
                             labels = c("C","F","I","L"),
                             font.label = list(size = 18),
                             hjust = -0.25)

int_value_combo <- annotate_figure(int_value_combo, bottom = textGrob("Age (d)", gp = gpar(fontsize = 24)),
                                left = textGrob("Intrinsic value", rot = 90, gp = gpar(fontsize = 24)))


big_boi <- ggarrange(LRO_summary, RO_curves, int_value_combo,
                     ncol = 3, nrow = 1)
big_boi
dev.off()

# Reproductive period---------------------

# Loading and reformatting data

rep_data <- read_csv("rep_period.csv")

# Changing the order of the factors
rep_data <- rep_data[with(rep_data, order(generation, maternal_age, strain)),]

# Removing censored individuals (we did not measure their full reproductive period)
rep_data_no_censor <- filter(rep_data, censor != 0)

# Modifying the dataset from wide to long format (needed for plotting and analysis):
rep_long_nc <- melt(data = rep_data_no_censor,
                    id.vars = c("generation", "maternal_age", "strain", "plate", "well"),
                    variable.name = "day",
                    value.name = "num_rep_days")
rep_long_nc <- filter(rep_long_nc, day != "censor")

# Changing order of maternal age levels
rep_long_nc$maternal_age <- factor(rep_long_nc$maternal_age, levels = c("young","mid","old"))

# Calculating the reproductive period in days and as a percent of the lifespan, per individual
rep_per_ind <- rep_long_nc %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(rep_time = sum(num_rep_days, na.rm = TRUE),
            length = length(num_rep_days[!is.na(num_rep_days)]),
            prop_rep = (rep_time/length),
            percent_rep = (rep_time/length)*100)

# Calculating summary statistics across individuals
rep_summary <- rep_per_ind %>%
  group_by(generation, maternal_age, strain) %>%
  summarize(
    mean_perc=mean(percent_rep, na.rm = TRUE), 
    med_perc=median(percent_rep, na.rm = TRUE),
    SE_perc=sd(percent_rep, na.rm =TRUE)/sqrt(length(!is.na(percent_rep))),
    mean_days=mean(rep_time, na.rm=TRUE),
    SE_days=sd(rep_time, na.rm =TRUE)/sqrt(length(!is.na(rep_time))),
    med_days=median(rep_time, na.rm=TRUE),
    n = length(percent_rep))
rep_summary

rep_per_ind$maternal_age <- factor(rep_per_ind$maternal_age)
rep_summary$maternal_age <- factor(rep_summary$maternal_age)

rep_per_ind$final_line <- paste(rep_per_ind$generation,
                                rep_per_ind$maternal_age,
                                rep_per_ind$strain)

rep_summary$final_line <- paste(rep_summary$generation,
                                rep_summary$maternal_age,
                                rep_summary$strain)

# Analysis: GLM

# Excluding the F0 generation
rep_per_ind_F1 <- filter(rep_per_ind, generation != "F0")

# Checking data for assumptions of normality and homogeneity of variance
rep_aov <- aov((prop_rep) ~ maternal_age * strain, rep_per_ind_F1)
rep_resid <- residuals(rep_aov)

shapiro.test(rep_resid)
leveneTest((prop_rep) ~ maternal_age * strain, rep_per_ind_F1)

# Assumptions are not met, can't find any useful transformations (proportion data)
# Will use GLM instead
# Rep. period data are also overdispersed- using quasibinomial GLM

qbin.glm <- glm(prop_rep ~ maternal_age * strain, family = quasibinomial, data = rep_per_ind_F1, weights = length)
summary(qbin.glm)

drop1(qbin.glm, test = "F") # interaction is significant- keep in model

# Model validation (plotting residuals)
plot(qbin.glm)

# Pairwise comparisons: which maternal age groups differ within strains?

EMM <- emmeans(qbin.glm, ~ maternal_age * strain)
EMM    # display the cell means

### Simple pairwise comparisons...
pairs(EMM, simple = "maternal_age")    # compare maternal age cohorts for each strain

# Figure

# BmanL5:

wilcox_L5 <- c("","a","b","b") # from pairwise comparisons

rep_per_ind_L5 <- filter(rep_per_ind, strain == "L5")

rep_sum_L5 <- filter(rep_summary, strain == "L5")

rep_sum_L5$final_line <- factor(rep_sum_L5$final_line, levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"), labels = c("F0", "Y", "M","O"))

rep_per_ind_L5$final_line <- factor(rep_per_ind_L5$final_line, levels = c("F0 NA L5", "F1 young L5", "F1 mid L5", "F1 old L5"), labels = c("F0", "Y", "M","O"))

# Plotting the mean and SE
rep_raw <- ggplot(data = rep_sum_L5, mapping = aes(x = final_line, y = mean_perc))+
  geom_point(size = 3) +
  geom_errorbar(data = rep_sum_L5, 
                aes(x = final_line, ymin=mean_perc-SE_perc, ymax=mean_perc+SE_perc), width=0.1) +
  geom_text(label=wilcox_L5, 
            nudge_x = 0.3,
            size = 8)

# Plotting the raw data points
rep_combo <- rep_raw +
  geom_jitter(data = rep_per_ind_L5, aes(x = final_line, y = percent_rep),
              width=.2, alpha=0.2, size = 2) +
  geom_label(
    label="BmanL5", 
    size = 8,
    x=1,
    y=4,
    label.size = 0,
    color = "black")

rep_L5_plot <- rep_combo + 
  ylab(element_blank()) +
  xlab(element_blank()) +
  ylim(c(-0.5,100.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=14))
rep_L5_plot

# BmanRUS:

wilcox_RUS <- c("","a","b") # from pairwise comparisons

rep_per_ind_RUS <- filter(rep_per_ind, strain == "RUS")

rep_sum_RUS <- filter(rep_summary, strain == "RUS")

rep_sum_RUS$final_line <- factor(rep_sum_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"), labels = c("F0", "Y","O"))

rep_per_ind_RUS$final_line <- factor(rep_per_ind_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"), labels = c("F0", "Y","O"))

# Plotting the mean and SE
rep_raw <- ggplot(data = rep_sum_RUS, mapping = aes(x = final_line, y = mean_perc))+
  geom_point(size = 3) +
  geom_errorbar(data = rep_sum_RUS, 
                aes(x = final_line, ymin=mean_perc-SE_perc, ymax=mean_perc+SE_perc), width=0.1) +
  geom_text(label=wilcox_RUS, 
            nudge_x = 0.3,
            size = 8)

# Plotting the raw data points
rep_combo <- rep_raw +
  geom_jitter(data = rep_per_ind_RUS, aes(x = final_line, y = percent_rep),
              width=.2, alpha=0.2, size = 2) +
  geom_label(
    label="BmanRUS", 
    size = 8,
    x=0.95,
    y=4,
    label.size = 0,
    color = "black")

rep_RUS_plot <- rep_combo + 
  ylab(element_blank()) +
  xlab(element_blank()) +
  ylim(c(-0.5,100.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=14))
rep_RUS_plot

# BmanRUS-RE:

rep_per_ind_RUS_RE <- filter(rep_per_ind, strain == "RUS-RE")

rep_sum_RUS_RE <- filter(rep_summary, strain == "RUS-RE")

rep_sum_RUS_RE$final_line <- factor(rep_sum_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"), labels = c("F0", "Y", "M","O"))

rep_per_ind_RUS_RE$final_line <- factor(rep_per_ind_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 young RUS-RE", "F1 mid RUS-RE", "F1 old RUS-RE"), labels = c("F0", "Y", "M","O"))

# Plotting the mean and SE
rep_raw <- ggplot(data = rep_sum_RUS_RE, mapping = aes(x = final_line, y = mean_perc))+
  geom_point(size = 3) +
  geom_errorbar(data = rep_sum_RUS_RE, 
                aes(x = final_line, ymin=mean_perc-SE_perc, ymax=mean_perc+SE_perc), width=0.1)

# Plotting the raw data points
rep_combo <- rep_raw +
  geom_jitter(data = rep_per_ind_RUS_RE, aes(x = final_line, y = percent_rep),
              width=.2, alpha=0.2, size = 2) +
  geom_label(
    label="BmanRUS-RE", 
    size = 8,
    x=1.4,
    y=4,
    label.size = 0,
    color = "black")

rep_RUS_RE_plot <- rep_combo + 
  ylab(element_blank()) +
  xlab(element_blank()) +
  ylim(c(-0.5,100.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=14))
rep_RUS_RE_plot

# BpL1:

wilcox_bp <- c("","a","b","ab") # from pairwise comparisons

rep_per_ind_BpL1 <- filter(rep_per_ind, strain == "BpL1")

rep_sum_BpL1 <- filter(rep_summary, strain == "BpL1")

rep_sum_BpL1$final_line <- factor(rep_sum_BpL1$final_line, levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"), labels = c("F0", "Y", "M","O"))

rep_per_ind_BpL1$final_line <- factor(rep_per_ind_BpL1$final_line, levels = c("F0 NA BpL1", "F1 young BpL1", "F1 mid BpL1", "F1 old BpL1"), labels = c("F0", "Y", "M","O"))

# Plotting the mean and SE
rep_raw <- ggplot(data = rep_sum_BpL1, mapping = aes(x = final_line, y = mean_perc))+
  geom_point(size = 3) +
  geom_errorbar(data = rep_sum_BpL1, 
                aes(x = final_line, ymin=mean_perc-SE_perc, ymax=mean_perc+SE_perc), width=0.1) +
  geom_text(label=wilcox_bp, 
            nudge_x = 0.3,
            size = 8)

#plotting the raw data points
rep_combo <- rep_raw +
  geom_jitter(data = rep_per_ind_BpL1, aes(x = final_line, y = percent_rep),
              width=.2, alpha=0.2, size = 2) +
  geom_label(
    label="BpL1", 
    size = 8,
    x=0.85,
    y=5,
    label.size = 0,
    color = "black")

rep_Bp_plot <- rep_combo + 
  ylab(element_blank()) +
  xlab(element_blank()) +
  ylim(c(-0.5,100.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text=element_text(size=14))
rep_Bp_plot

# Reproducing the manuscript figure:

tiff("rep_period.tiff", units="in", width=6, height=12, res=300)

rep_curves <- ggarrange(rep_Bp_plot, rep_RUS_plot,
                        rep_RUS_RE_plot, rep_L5_plot,
                        ncol = 1, nrow = 4)

annotate_figure(rep_curves, bottom = textGrob("Cohort", gp = gpar(fontsize = 24)), left = textGrob("Reproductive period (% of lifespan)", 
                                                                                                   rot = 90, gp = gpar(fontsize = 24)))

dev.off()

