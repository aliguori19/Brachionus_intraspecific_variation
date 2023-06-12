# Code for visualizing and analyzing reproductive period data: measured as days carrying eggs (and as a percent of the lifespan)

# Required packages:

library(tidyverse)
library(reshape2)
library(car)
library(grid)
library(ggpubr)
library(emmeans)

# Loading and reformatting data--------------------------------------------------------------------

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

# Statistics: -------------------------------------------------------------------

# Excluding the F0 generation
rep_per_ind_F1 <- filter(rep_per_ind, generation != "F0")

# Checking n for each group
balance <- rep_per_ind_F1 %>% count(maternal_age, strain)

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

# Figure: ----------------------------------------------------------------------

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

rep_per_ind_RUS <- filter(rep_per_ind, strain == "RUS")

rep_sum_RUS <- filter(rep_summary, strain == "RUS")

rep_sum_RUS$final_line <- factor(rep_sum_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"), labels = c("F0", "Y","O"))

rep_per_ind_RUS$final_line <- factor(rep_per_ind_RUS$final_line, levels = c("F0 NA RUS", "F1 young RUS", "F1 old RUS"), labels = c("F0", "Y","O"))

# Plotting the mean and SE
rep_raw <- ggplot(data = rep_sum_RUS, mapping = aes(x = final_line, y = mean_perc))+
  geom_point(size = 3) +
  geom_errorbar(data = rep_sum_RUS, 
                aes(x = final_line, ymin=mean_perc-SE_perc, ymax=mean_perc+SE_perc), width=0.1)

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
