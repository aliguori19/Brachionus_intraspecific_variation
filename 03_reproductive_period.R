# Code for visualizing and analyzing reproductive period data: measured as days carrying eggs (and as a percent of the lifespan)

# Required packages:

library(tidyverse)
library(reshape2)
library(car)
library(grid)
library(ggpubr)

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

# Calculating the reproductive period in days and as a percent of the lifespan, per individual
rep_per_ind <- rep_long_nc %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(rep_time = sum(num_rep_days, na.rm = TRUE),
            length = length(num_rep_days[!is.na(num_rep_days)]),
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

# BmanL5:
rep_per_ind_L5 <- filter(rep_per_ind, strain == "L5")

kruskal.test(percent_rep ~ final_line, data = rep_per_ind_L5)
pairwise.wilcox.test(rep_per_ind_L5$percent_rep, rep_per_ind_L5$final_line,
                     p.adjust.method = "BH")

# BmanRUS:
rep_per_ind_RUS <- filter(rep_per_ind, strain == "RUS")

kruskal.test(percent_rep ~ final_line, data = rep_per_ind_RUS)
pairwise.wilcox.test(rep_per_ind_RUS$percent_rep, rep_per_ind_RUS$final_line,
                     p.adjust.method = "BH")

# BmanRUS-RE:
rep_per_ind_RUS_RE <- filter(rep_per_ind, strain == "RUS-RE")

kruskal.test(percent_rep ~ final_line, data = rep_per_ind_RUS_RE)
pairwise.wilcox.test(rep_per_ind_RUS_RE$percent_rep, rep_per_ind_RUS_RE$final_line,
                     p.adjust.method = "BH")

# BpL1:
rep_per_ind_BpL1 <- filter(rep_per_ind, strain == "BpL1")

kruskal.test(percent_rep ~ final_line, data = rep_per_ind_BpL1)
pairwise.wilcox.test(rep_per_ind_BpL1$percent_rep, rep_per_ind_BpL1$final_line,
                     p.adjust.method = "BH")

# Figures: ----------------------------------------------------------------------

# BmanL5:

wilcox_L5 <- c("a","a","b","b")

rep_sum_L5 <- filter(rep_summary, strain == "L5")

rep_sum_L5$final_line <- factor(rep_sum_L5$final_line, levels = c("F0 NA L5", "F1 3 L5", "F1 6 L5", "F1 10 L5"), labels = c("F0", "Y", "M","O"))

rep_per_ind_L5$final_line <- factor(rep_per_ind_L5$final_line, levels = c("F0 NA L5", "F1 3 L5", "F1 6 L5", "F1 10 L5"), labels = c("F0", "Y", "M","O"))

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
  xlab("Cohort") +
  ylim(c(-0.5,100.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.title=element_text(size=24), 
        axis.text=element_text(size=14))
rep_L5_plot

# BmanRUS:

wilcox_RUS <- c("a","a","b")

rep_sum_RUS <- filter(rep_summary, strain == "RUS")

rep_sum_RUS$final_line <- factor(rep_sum_RUS$final_line, levels = c("F0 NA RUS", "F1 3 RUS", "F1 11 RUS"), labels = c("F0", "Y","O"))

rep_per_ind_RUS$final_line <- factor(rep_per_ind_RUS$final_line, levels = c("F0 NA RUS", "F1 3 RUS", "F1 11 RUS"), labels = c("F0", "Y","O"))

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

rep_sum_RUS_RE <- filter(rep_summary, strain == "RUS-RE")

rep_sum_RUS_RE$final_line <- factor(rep_sum_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 3 RUS-RE", "F1 6 RUS-RE", "F1 10 RUS-RE"), labels = c("F0", "Y", "M","O"))

rep_per_ind_RUS_RE$final_line <- factor(rep_per_ind_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 3 RUS-RE", "F1 6 RUS-RE", "F1 10 RUS-RE"), labels = c("F0", "Y", "M","O"))

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

wilcox_bp <- c("ab","a","b","a")

rep_sum_BpL1 <- filter(rep_summary, strain == "BpL1")

rep_sum_BpL1$final_line <- factor(rep_sum_BpL1$final_line, levels = c("F0 NA BpL1", "F1 3 BpL1", "F1 6 BpL1", "F1 9 BpL1"), labels = c("F0", "Y", "M","O"))

rep_per_ind_BpL1$final_line <- factor(rep_per_ind_BpL1$final_line, levels = c("F0 NA BpL1", "F1 3 BpL1", "F1 6 BpL1", "F1 9 BpL1"), labels = c("F0", "Y", "M","O"))

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

tiff("rep_period.tiff", units="in", width=6, height=12, res=1200)

rep_curves <- ggarrange(rep_Bp_plot, rep_RUS_plot,
                        rep_RUS_RE_plot, rep_L5_plot,
                        ncol = 1, nrow = 4,
                        labels = c())#,
                        #hjust = -0.25)

annotate_figure(rep_curves, left = textGrob("Reproductive period (% of lifespan)", 
                                            rot = 90, gp = gpar(fontsize = 24)))

dev.off()
