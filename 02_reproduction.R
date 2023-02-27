# Code for analyzing lifetime reproductive output (LRO), maximum daily reproduction (MDR), time of 50% LRO

# Required packages:
library(tidyverse)
library(reshape2)
library(car)

# Loading and reformatting data--------------------------------------------------------------------

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

# Calculating and analyzing lifetime reproductive output (LRO) per individual-----------------------

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

LRO_summary$final_line <- paste(LRO_summary$generation, # needed for plotting later
                                LRO_summary$maternal_age,
                                LRO_summary$strain)

LRO_summary$maternal_age <- factor(LRO_summary$maternal_age)

# Statistical analyses for LRO - Kruskal-Wallis rank sum tests and pairwise Wilcoxon rank sum tests, conducted within each strain:

# BmanL5:
LRO_L5 <- filter(LRO_per_ind, strain == "L5")

kruskal.test(LRO ~ final_line, data = LRO_L5)
pairwise.wilcox.test(LRO_L5$LRO, LRO_L5$final_line,
                     p.adjust.method = "BH")

# BmanRUS:
LRO_RUS <- filter(LRO_per_ind, strain == "RUS")

kruskal.test(LRO ~ final_line, data = LRO_RUS)
pairwise.wilcox.test(LRO_RUS$LRO, LRO_RUS$final_line,
                     p.adjust.method = "BH")

# BmanRUS-RE:
LRO_RUS_RE <- filter(LRO_per_ind, strain == "RUS-RE")

kruskal.test(LRO ~ final_line, data = LRO_RUS_RE)
pairwise.wilcox.test(LRO_RUS_RE$LRO, LRO_RUS_RE$final_line,
                     p.adjust.method = "BH")

# BpL1:
LRO_BpL1 <- filter(LRO_per_ind, strain == "BpL1")

kruskal.test(LRO ~ final_line, data = LRO_BpL1)
pairwise.wilcox.test(LRO_BpL1$LRO, LRO_BpL1$final_line,
                     p.adjust.method = "BH")

# LRO Figures--------------------------------------------------------------

# BmanL5:
wilcox_L5 <- c("a","a","b","b") # from pairwise analyses

LRO_sum_L5 <- filter(LRO_summary, strain == "L5")

LRO_sum_L5$final_line <- factor(LRO_sum_L5$final_line, levels = c("F0 NA L5", "F1 3 L5", "F1 6 L5", "F1 10 L5"), labels = c("F0", "Y", "M","O"))

LRO_per_ind_L5 <- filter(LRO_per_ind, strain == "L5")

LRO_per_ind_L5$final_line <- factor(LRO_per_ind_L5$final_line, levels = c("F0 NA L5", "F1 3 L5", "F1 6 L5", "F1 10 L5"), labels = c("F0", "Y", "M","O"))

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
  xlab("Cohort") +
  ylim(c(-0.5,38)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.title=element_text(size=24), 
        axis.text=element_text(size=14))
LRO_L5_plot

# BmanRUS:
LRO_sum_RUS <- filter(LRO_summary, strain == "RUS")

LRO_sum_RUS$final_line <- factor(LRO_sum_RUS$final_line, levels = c("F0 NA RUS", "F1 3 RUS", "F1 11 RUS"), labels = c("F0", "Y","O"))

LRO_per_ind_RUS <- filter(LRO_per_ind, strain == "RUS")

LRO_per_ind_RUS$final_line <- factor(LRO_per_ind_RUS$final_line, levels = c("F0 NA RUS", "F1 3 RUS", "F1 11 RUS"), labels = c("F0", "Y","O"))

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
LRO_sum_RUS_RE <- filter(LRO_summary, strain == "RUS-RE")

LRO_sum_RUS_RE$final_line <- factor(LRO_sum_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 3 RUS-RE", "F1 6 RUS-RE", "F1 10 RUS-RE"), labels = c("F0", "Y", "M","O"))

LRO_per_ind_RUS_RE <- filter(LRO_per_ind, strain == "RUS-RE")

LRO_per_ind_RUS_RE$final_line <- factor(LRO_per_ind_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 3 RUS-RE", "F1 6 RUS-RE", "F1 10 RUS-RE"), labels = c("F0", "Y", "M","O"))

# Plotting means and SE
LRO_raw <- ggplot(data = LRO_sum_RUS_RE, mapping = aes(x = final_line, y = mean))+
  geom_point(size = 3) +
  geom_errorbar(data = LRO_sum_RUS_RE, 
                aes(x = final_line, ymin=mean-SE, ymax=mean+SE), width=0.1)

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
wilcox_bp <- c("a","b","a","ab") # from pairwise analyses

LRO_sum_BpL1 <- filter(LRO_summary, strain == "BpL1")

LRO_sum_BpL1$final_line <- factor(LRO_sum_BpL1$final_line, levels = c("F0 NA BpL1", "F1 3 BpL1", "F1 6 BpL1", "F1 9 BpL1"), labels = c("F0", "Y", "M","O"))

LRO_per_ind_BpL1 <- filter(LRO_per_ind, strain == "BpL1")

LRO_per_ind_BpL1$final_line <- factor(LRO_per_ind_BpL1$final_line, levels = c("F0 NA BpL1", "F1 3 BpL1", "F1 6 BpL1", "F1 9 BpL1"), labels = c("F0", "Y", "M","O"))

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

# Plotting daily reproduction over time---------------------------------------------

# Summarizing reproductive output across individuals, per day - includes censored individuals (their daily reproduction before loss is informative here)

# Changing the original dataset (including censored individuals) from wide to long format
neo_long <- melt(data = neonate_data,
                    id.vars = c("generation", "maternal_age", "strain", "plate", "well"),
                    variable.name = "day",
                    value.name = "num_neonates")

neo_long <- filter(neo_long, day != "censor")

# Calculating the mean number of neonates produced across all individuals, per day
RO_per_day <- neo_long %>%
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

RO_L5$final_line <- factor(RO_L5$final_line, levels = c("F0 NA L5", "F1 3 L5", "F1 6 L5", "F1 10 L5"))

RO_raw <- ggplot(data=RO_L5, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw+
  geom_line(data=RO_L5, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black", "#999999", "#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid", "solid"),labels = c("F0","Y","M","O")) +  geom_label(
    label="BmanL5", 
    size = 8,
    x=24,
    y=4.5,
    label.size = 0,
    color = "black")

RO_L5_plot_time <- RO_combo + 
  ylab("") +
  ylim(c(0,5)) +
  xlab("Age (d)") +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=14),
        legend.position = "none")
RO_L5_plot_time

# BmanRUS:
RO_RUS <- filter(RO_per_day, strain == "RUS")

RO_RUS$final_line <- factor(RO_RUS$final_line, levels = c("F0 NA RUS", "F1 3 RUS", "F1 11 RUS"))

RO_raw <- ggplot(data=RO_RUS, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw+
  geom_line(data=RO_RUS, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black","#0072B2"), labels = c("F0","Y", "O")) +
  scale_linetype_manual(values = c("dotted","solid","solid"),labels = c("F0","Y", "O")) +
  geom_label(
    label="BmanRUS", 
    size = 8,
    x=23,
    y=4.5,
    label.size = 0,
    color = "black")

RO_RUS_plot_time <- RO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(0,5)) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=16, face="bold"), 
        axis.text=element_text(size=14),
        legend.position = "none")
RO_RUS_plot_time

# BmanRUS-RE:
RO_RUS_RE <- filter(RO_per_day, strain == "RUS-RE")

RO_RUS_RE$final_line <- factor(RO_RUS_RE$final_line, levels = c("F0 NA RUS-RE", "F1 3 RUS-RE", "F1 6 RUS-RE", "F1 10 RUS-RE"))

RO_raw <- ggplot(data=RO_RUS_RE, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw+
  geom_line(data=RO_RUS_RE, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black","#999999","#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid","solid"),labels = c("F0","Y","M","O")) + geom_label(
    label="BmanRUS-RE", 
    size = 8,
    x=22,
    y=4.5,
    label.size = 0,
    color = "black")

RO_RUS_RE_time <- RO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(0,5)) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=16, face="bold"), 
        axis.text=element_text(size=14),
        legend.position = "none")
RO_RUS_RE_time

# BpL1:
RO_BpL1 <- filter(RO_per_day, strain == "BpL1")

RO_BpL1$final_line <- factor(RO_BpL1$final_line, levels = c("F0 NA BpL1", "F1 3 BpL1", "F1 6 BpL1", "F1 9 BpL1"))

RO_raw <- ggplot(data=RO_BpL1, aes(x=day, y=mean)) + 
  geom_point(aes(group = final_line), size = 1.5) +
  geom_errorbar(aes(x=day, ymin=mean-SE, ymax=mean+SE), width=0.1)

RO_combo <- RO_raw +
  geom_line(data=RO_BpL1, aes(x=day, y=mean, color = final_line, linetype = final_line, group=final_line), size = 1) +
  scale_color_manual(values = c("black", "black","#999999","#0072B2"), labels = c("F0","Y","M","O")) +
  scale_linetype_manual(values = c("dotted","solid","solid","solid"),labels = c("F0","Y","M","O")) + geom_label(
    label="BpL1", 
    size = 8,
    x=25,
    y=4.5,
    label.size = 0,
    color = "black")

RO_Bp_plot_time <- RO_combo + 
  ylab("") +
  xlab(element_blank()) +
  ylim(c(0,5)) +
  scale_x_discrete(breaks=seq(0,28,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=18, face="bold"), 
        axis.text=element_text(size=14),
        legend.position = c(.85,.5),
        legend.title = element_blank(),
        legend.key.width = unit(3, "line"),
        legend.text=element_text(size=16))
RO_Bp_plot_time

# Combined reproductive output multipanel figure (reproducing figure in the manuscript)--------------

tiff("LRO.tiff", units="in", width=12, height=12, res=1200)

RO_curves <- ggarrange(RO_Bp_plot_time, RO_RUS_plot_time,
                        RO_RUS_RE_time, RO_L5_plot_time, 
                        ncol = 1, nrow = 4,
                        labels = c("B","D","F","H"),
                        font.label = list(size = 18),
                        hjust = -0.25)

RO_curves <- annotate_figure(RO_curves, left = textGrob(expression(paste("Offspring ind"^"-1"*"day"^"-1")), 
                                                        rot = 90, gp = gpar(fontsize = 24)))

LRO_summary <- ggarrange(LRO_Bp_plot, LRO_RUS_plot,
                            LRO_RUS_RE_plot, LRO_L5_plot, 
                            ncol = 1, nrow = 4,
                            labels = c("A","C","E","G"),
                            font.label = list(size = 18),
                            hjust = -0.25,
                            common.legend = TRUE, legend = "top")

LRO_summary <- annotate_figure(LRO_summary, left = textGrob(expression(paste("Lifetime reproductive output (offspring ind"^"-1"*")")), 
                                                            rot = 90, gp = gpar(fontsize = 24)))

big_boi <- ggarrange(LRO_summary, RO_curves,
                     ncol = 2, nrow = 1)
big_boi
dev.off()

# Maximum daily reproduction (MDR) analyses--------------------------------------------------

# BmanL5:
kruskal.test(max_RO ~ final_line, data = LRO_L5)
pairwise.wilcox.test(LRO_L5$max_RO, LRO_L5$final_line,
                     p.adjust.method = "BH")

# BmanRUS:
kruskal.test(max_RO ~ final_line, data = LRO_RUS)
pairwise.wilcox.test(LRO_RUS$max_RO, LRO_RUS$final_line,
                     p.adjust.method = "BH")

# BmanRUS-RE:
kruskal.test(max_RO ~ final_line, data = LRO_RUS_RE)
pairwise.wilcox.test(LRO_RUS_RE$max_RO, LRO_RUS_RE$final_line,
                     p.adjust.method = "BH")

# BpL1:
kruskal.test(max_RO ~ final_line, data = LRO_BpL1)
pairwise.wilcox.test(LRO_BpL1$max_RO, LRO_BpL1$final_line,
                     p.adjust.method = "BH")

# Calculating and analyzing timing of the production of 50% of LRO-----------------------------

# Calculating cumulative reproduction PER INDIVIDUAL:
# Then getting the day when 50% of lifetime fecundity was reached:
  
neo_long_nc <- neo_long_nc[with(neo_long_nc, order(generation, maternal_age, strain, plate, well)),]

neo_long_nc <- neo_long_nc %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  mutate(cum_RO = cumsum(num_neonates))

# Removing individuals with a LRO of zero (not relevant here)
LRO_per_ind <- filter(LRO_per_ind, LRO != 0)

neo_merge <- merge(LRO_per_ind, neo_long_nc, by = c("generation","maternal_age","strain","plate","well"))

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

# Statistics:

sub_summ$final_line <- paste(sub_summ$generation,
                                sub_summ$maternal_age,
                                sub_summ$strain)
# BmanL5:
fifty_L5 <- filter(sub_summ, strain == "L5")

kruskal.test(time_fifty ~ final_line, data = fifty_L5)
pairwise.wilcox.test(fifty_L5$time_fifty, fifty_L5$final_line,
                     p.adjust.method = "BH")

# BmanRUS:
fifty_RUS <- filter(sub_summ, strain == "RUS")

kruskal.test(time_fifty ~ final_line, data = fifty_RUS)
pairwise.wilcox.test(fifty_RUS$time_fifty, fifty_RUS$final_line,
                     p.adjust.method = "BH")

# BmanRUS-RE:
fifty_RUS_RE <- filter(sub_summ, strain == "RUS-RE")

kruskal.test(time_fifty ~ final_line, data = fifty_RUS_RE)
pairwise.wilcox.test(fifty_RUS_RE$time_fifty, fifty_RUS_RE$final_line,
                     p.adjust.method = "BH")

# BpL1:
fifty_BpL1 <- filter(sub_summ, strain == "BpL1")

kruskal.test(time_fifty ~ final_line, data = fifty_BpL1)
pairwise.wilcox.test(fifty_BpL1$time_fifty, fifty_BpL1$final_line,
                     p.adjust.method = "BH")
# F0 differs from the F1 generation, but no differences within the F1 cohorts

