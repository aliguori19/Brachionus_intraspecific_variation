# Code for analyzing and visualizing relationships between lifespan and reproduction (trade-offs present?)

# Required packages:
library(tidyverse)
library(reshape2)
library(survival)
library(survminer)
library(lsmeans)
library(car)
library(grid)
library(ggpubr)

# Loading and reformatting data--------------------------------------------------------------------

surv_data <- read_csv("survival.csv")

# Changing the order of the factors and removing censored individuals
surv_data <- surv_data[with(surv_data, order(generation, maternal_age, strain)),]
surv_data_no_censor <- filter(surv_data, censor != 0)

neonate_data <- read_csv("neonates.csv")
neonate_data <- neonate_data[with(neonate_data, order(generation, maternal_age, strain)),]
neonate_data_no_censor <- filter(neonate_data, censor != 0)

rep_data <- read_csv("rep_period.csv")
rep_data <- rep_data[with(rep_data, order(generation, maternal_age, strain)),]
rep_data_no_censor <- filter(rep_data, censor != 0)

# Modifying the datasets from wide to long format (needed for plotting and analysis):
surv_long <- melt(data = surv_data_no_censor, 
                  id.vars = c("generation", "maternal_age", "strain", 
                              "plate", "well"),
                  variable.name = "day",
                  value.name = "dead_alive")
surv_long <- filter(surv_long, day != "censor")

neo_long <- melt(data = neonate_data_no_censor,
                 id.vars = c("generation", "maternal_age", "strain", "plate", "well"),
                 variable.name = "day",
                 value.name = "num_neonates")
neo_long <- filter(neo_long, day != "censor")

rep_long <- melt(data = rep_data_no_censor,
                    id.vars = c("generation", "maternal_age", "strain", "plate", "well"),
                    variable.name = "day",
                    value.name = "num_rep_days")
rep_long <- filter(rep_long, day != "censor")

# Per individual calculations: ------------------------------------------------------------------

# Calculating lifespan per individual:
survival_data_per_ind <- surv_long %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(surv_time = sum(dead_alive, na.rm = TRUE))

# Calculating LRO and MDR per individual:
LRO_per_ind <- neo_long %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(LRO = sum(num_neonates, na.rm = TRUE),
            max_RO = max(num_neonates, na.rm = TRUE),
            fifty_percent_LT = (sum(num_neonates, na.rm = TRUE))*0.5)

# Calculating time of 50% LRO production
# Excluding individuals with an LRO of zero:
LRO_per_ind_excl0 <- filter(LRO_per_ind, LRO != 0)

neo_long <- neo_long[with(neo_long, order(generation, maternal_age, strain, plate, well)),]

neo_long <- neo_long %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  mutate(cum_LRO = cumsum(num_neonates))

neo_merge <- merge(LRO_per_ind_excl0, neo_long, by = c("generation","maternal_age","strain","plate","well"))

sub <- subset(neo_merge, cum_LRO >= fifty_percent_LT , select = c(generation,maternal_age,strain,plate,well,LRO,fifty_percent_LT,day,cum_LRO))

sub$day <- as.numeric(sub$day)
sub_summ <- sub %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(time_fifty = min(day))

# Calculating reproductive period per individual:
rep_per_ind <- rep_long %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(rep_time = sum(num_rep_days, na.rm = TRUE),
            length = length(num_rep_days[!is.na(num_rep_days)]),
            percent_rep = (rep_time/length)*100)

# Exploring tradeoffs: LRO versus lifespan: ---------------------------------------------------

# Making a dataframe with lifespan and LRO data:
tradeoff_df <- cbind(survival_data_per_ind, LRO_per_ind$LRO)
tradeoff_df <- rename(tradeoff_df, "LRO" = "...7")
tradeoff_df_F1 <- filter(tradeoff_df, generation != "F0")

tradeoff_df_F1$strain <- factor(tradeoff_df_F1$strain, levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

# Exploring relationships among maternal age cohorts, within each strain:
tradeoff_df_F1_L5 <- filter(tradeoff_df_F1, strain == "BmanL5")
tradeoff_df_F1_RUS <- filter(tradeoff_df_F1, strain == "BmanRUS")
tradeoff_df_F1_RUS_RE <- filter(tradeoff_df_F1, strain == "BmanRUS-RE")
tradeoff_df_F1_BpL1 <- filter(tradeoff_df_F1, strain == "BpL1")

# BmanL5:

# Models for each maternal age:
tradeoff_df_F1_L5_Y <- filter(tradeoff_df_F1_L5, maternal_age == "young")
tradeoff_df_F1_L5_M <- filter(tradeoff_df_F1_L5, maternal_age == "mid")
tradeoff_df_F1_L5_O <- filter(tradeoff_df_F1_L5, maternal_age == "old")

lm_L5_Y <- lm(LRO ~ surv_time, data = tradeoff_df_F1_L5_Y)
summary(lm_L5_Y)
lm_L5_M <- lm(LRO ~ surv_time, data = tradeoff_df_F1_L5_M)
summary(lm_L5_M)
lm_L5_O <- lm(LRO ~ surv_time, data = tradeoff_df_F1_L5_O)
summary(lm_L5_O)

# Statistics (ANCOVA):
tradeoff_df_F1_L5$maternal_age <- factor(tradeoff_df_F1_L5$maternal_age, levels = c("young", "mid", "old"), labels = c(
  "Y", "M","O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_L5)
anova(m.interaction)

# Visualization:
LRO_L5 <- ggplot(tradeoff_df_F1_L5, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  geom_label(
    label="L5", 
    size = 8.5,
    x=6,
    y=30,
    label.size = 0,
    color = "black") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_L5

# BmanRUS:

# Models for each maternal age:
tradeoff_df_F1_RUS_Y <- filter(tradeoff_df_F1_RUS, maternal_age == "young")
tradeoff_df_F1_RUS_O <- filter(tradeoff_df_F1_RUS, maternal_age == "old")

lm_RUS_Y <- lm(LRO ~ surv_time, data = tradeoff_df_F1_RUS_Y)
summary(lm_RUS_Y)
lm_RUS_O <- lm(LRO ~ surv_time, data = tradeoff_df_F1_RUS_O)
summary(lm_RUS_O)

# Statistics (ANCOVA):
tradeoff_df_F1_RUS$maternal_age <- factor(tradeoff_df_F1_RUS$maternal_age, levels = c("young", "old"), labels = c(
  "Y", "O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_RUS)
anova(m.interaction)

# Visualization:
LRO_RUS <- ggplot(tradeoff_df_F1_RUS, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  geom_label(
    label="RUS", 
    size = 8.5,
    x=6,
    y=32,
    label.size = 0,
    color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_RUS

# BmanRUS-RE:

# Models for each maternal age:
tradeoff_df_F1_RUS_RE_Y <- filter(tradeoff_df_F1_RUS_RE, maternal_age == "young")
tradeoff_df_F1_RUS_RE_M <- filter(tradeoff_df_F1_RUS_RE, maternal_age == "mid")
tradeoff_df_F1_RUS_RE_O <- filter(tradeoff_df_F1_RUS_RE, maternal_age == "old")

lm_RUS_RE_Y <- lm(LRO ~ surv_time, data = tradeoff_df_F1_RUS_RE_Y)
summary(lm_RUS_RE_Y)
lm_RUS_RE_M <- lm(LRO ~ surv_time, data = tradeoff_df_F1_RUS_RE_M)
summary(lm_RUS_RE_M)
lm_RUS_RE_O <- lm(LRO ~ surv_time, data = tradeoff_df_F1_RUS_RE_O)
summary(lm_RUS_RE_O)

# Statistics (ANCOVA):
tradeoff_df_F1_RUS_RE$maternal_age <- factor(tradeoff_df_F1_RUS_RE$maternal_age, 
                                             levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_RUS_RE)
anova(m.interaction)

# Visualization:
LRO_RUS_RE <- ggplot(tradeoff_df_F1_RUS_RE, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  geom_label(
    label="RUS-RE", 
    size = 8.5,
    x=9,
    y=39,
    label.size = 0,
    color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_RUS_RE

# BpL1:

# Models for each maternal age:
tradeoff_df_F1_Bp_Y <- filter(tradeoff_df_F1_BpL1, maternal_age == "young")
tradeoff_df_F1_Bp_M <- filter(tradeoff_df_F1_BpL1, maternal_age == "mid")
tradeoff_df_F1_Bp_O <- filter(tradeoff_df_F1_BpL1, maternal_age == "old")

lm_Bp_Y <- lm(LRO ~ surv_time, data = tradeoff_df_F1_Bp_Y)
summary(lm_Bp_Y)
lm_Bp_M <- lm(LRO ~ surv_time, data = tradeoff_df_F1_Bp_M)
summary(lm_Bp_M)
lm_Bp_O <- lm(LRO ~ surv_time, data = tradeoff_df_F1_Bp_O)
summary(lm_Bp_O)

# Statistics (ANCOVA):
tradeoff_df_F1_BpL1$maternal_age <- factor(tradeoff_df_F1_BpL1$maternal_age, 
                                           levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_BpL1)
anova(m.interaction)

# Visualization:
LRO_BpL1 <- ggplot(tradeoff_df_F1_BpL1, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  geom_label(
    label="BpL1", 
    size = 8.5,
    x=6,
    y=25,
    label.size = 0,
    color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_BpL1

# Exploring tradeoffs: Reproductive period (percentage of lifespan) versus lifespan: ---------------------------------------------------

# Making a dataframe with lifespan and reproductive period data:
tradeoff_df_rep <- cbind(survival_data_per_ind, rep_per_ind$percent_rep)
tradeoff_df_rep <- rename(tradeoff_df_rep, "rep_percent" = "...7")
tradeoff_df_rep_F1 <- filter(tradeoff_df_rep, generation != "F0")

tradeoff_df_rep_F1$strain <- factor(tradeoff_df_rep_F1$strain, levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

# Exploring relationships among maternal age cohorts, within each strain:
tradeoff_rep_L5 <- filter(tradeoff_df_rep_F1, strain == "BmanL5")
tradeoff_rep_RUS <- filter(tradeoff_df_rep_F1, strain == "BmanRUS")
tradeoff_rep_RUS_RE <- filter(tradeoff_df_rep_F1, strain == "BmanRUS-RE")
tradeoff_rep_BpL1 <- filter(tradeoff_df_rep_F1, strain == "BpL1")

# BmanL5:

# Models for each maternal age:
tradeoff_rep_L5_Y <- filter(tradeoff_rep_L5, maternal_age == "young")
tradeoff_rep_L5_M <- filter(tradeoff_rep_L5, maternal_age == "mid")
tradeoff_rep_L5_O <- filter(tradeoff_rep_L5, maternal_age == "old")

lm_rep_L5_Y <- lm(rep_percent ~ surv_time, data = tradeoff_rep_L5_Y)
summary(lm_rep_L5_Y)
lm_rep_L5_M <- lm(rep_percent ~ surv_time, data = tradeoff_rep_L5_M)
summary(lm_rep_L5_M)
lm_rep_L5_O <- lm(rep_percent ~ surv_time, data = tradeoff_rep_L5_O)
summary(lm_rep_L5_O)

# Statistics (ANCOVA):
tradeoff_rep_L5$maternal_age <- factor(tradeoff_rep_L5$maternal_age, 
                                       levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_L5)
anova(m.interaction)

# Visualization:
RPP_L5 <- ggplot(tradeoff_rep_L5, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_L5

# BmanRUS:

# Models for each maternal age:
tradeoff_rep_RUS_Y <- filter(tradeoff_rep_RUS, maternal_age == "young")
tradeoff_rep_RUS_O <- filter(tradeoff_rep_RUS, maternal_age == "old")

lm_rep_RUS_Y <- lm(rep_percent ~ surv_time, data = tradeoff_rep_RUS_Y)
summary(lm_rep_RUS_Y)
lm_rep_RUS_O <- lm(rep_percent ~ surv_time, data = tradeoff_rep_RUS_O)
summary(lm_rep_RUS_O)

# Statistics (ANCOVA):
tradeoff_rep_RUS$maternal_age <- factor(tradeoff_rep_RUS$maternal_age, 
                                        levels = c("young", "old"), labels = c("Y","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_RUS)
anova(m.interaction)

# Visualization:
RPP_RUS <- ggplot(tradeoff_rep_RUS, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_RUS

# BmanRUS-RE:

# Models for each maternal age:
tradeoff_rep_RUS_RE_Y <- filter(tradeoff_rep_RUS_RE, maternal_age == "young")
tradeoff_rep_RUS_RE_M <- filter(tradeoff_rep_RUS_RE, maternal_age == "mid")
tradeoff_rep_RUS_RE_O <- filter(tradeoff_rep_RUS_RE, maternal_age == "old")

lm_rep_RUS_RE_Y <- lm(rep_percent ~ surv_time, data = tradeoff_rep_RUS_RE_Y)
summary(lm_rep_RUS_RE_Y)
lm_rep_RUS_RE_M <- lm(rep_percent ~ surv_time, data = tradeoff_rep_RUS_RE_M)
summary(lm_rep_RUS_RE_M)
lm_rep_RUS_RE_O <- lm(rep_percent ~ surv_time, data = tradeoff_rep_RUS_RE_O)
summary(lm_rep_RUS_RE_O)

# Statistics (ANCOVA):
tradeoff_rep_RUS_RE$maternal_age <- factor(tradeoff_rep_RUS_RE$maternal_age, 
                                           levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_RUS_RE)
anova(m.interaction)

# Visualization:
RPP_RUS_RE <- ggplot(tradeoff_rep_RUS_RE, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_RUS_RE

# BpL1:

# Models for each maternal age:
tradeoff_rep_BpL1_Y <- filter(tradeoff_rep_BpL1, maternal_age == "young")
tradeoff_rep_BpL1_M <- filter(tradeoff_rep_BpL1, maternal_age == "mid")
tradeoff_rep_BpL1_O <- filter(tradeoff_rep_BpL1, maternal_age == "old")

lm_rep_BpL1_Y <- lm(rep_percent ~ surv_time, data = tradeoff_rep_BpL1_Y)
summary(lm_rep_BpL1_Y)
lm_rep_BpL1_M <- lm(rep_percent ~ surv_time, data = tradeoff_rep_BpL1_M)
summary(lm_rep_BpL1_M)
lm_rep_BpL1_O <- lm(rep_percent ~ surv_time, data = tradeoff_rep_BpL1_O)
summary(lm_rep_BpL1_O)

# Statistics (ANCOVA):
tradeoff_rep_BpL1$maternal_age <- factor(tradeoff_rep_BpL1$maternal_age, 
                                         levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_BpL1)
anova(m.interaction)

# Visualization:
RPP_BpL1 <- ggplot(tradeoff_rep_BpL1, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_BpL1

# Exploring tradeoffs: Timing of 50% LRO production versus lifespan: ---------------------------------------------------

# Making a dataframe with lifespan and 50% LRO data:
tradeoff_merge <- merge(survival_data_per_ind, sub_summ, by = c("generation","maternal_age","strain","plate","well"))
tradeoff_merge_F1 <- filter(tradeoff_merge, generation != "F0")

tradeoff_merge_F1$strain <- factor(tradeoff_merge_F1$strain, 
                                   levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

# Exploring relationships among maternal age cohorts, within each strain:
tradeoff_merge_L5 <- filter(tradeoff_merge_F1, strain == "BmanL5")
tradeoff_merge_RUS <- filter(tradeoff_merge_F1, strain == "BmanRUS")
tradeoff_merge_RUS_RE <- filter(tradeoff_merge_F1, strain == "BmanRUS-RE")
tradeoff_merge_BpL1 <- filter(tradeoff_merge_F1, strain == "BpL1")

# BmanL5:

# Models for each maternal age:
tradeoff_merge_L5_Y <- filter(tradeoff_merge_L5, maternal_age == "young")
tradeoff_merge_L5_M <- filter(tradeoff_merge_L5, maternal_age == "mid")
tradeoff_merge_L5_O <- filter(tradeoff_merge_L5, maternal_age == "old")

lm_50_L5_Y <- lm(time_fifty ~ surv_time, data = tradeoff_merge_L5_Y)
summary(lm_50_L5_Y)
lm_50_L5_M <- lm(time_fifty ~ surv_time, data = tradeoff_merge_L5_M)
summary(lm_50_L5_M)
lm_50_L5_O <- lm(time_fifty ~ surv_time, data = tradeoff_merge_L5_O)
summary(lm_50_L5_O)

# Statistics (ANCOVA):
tradeoff_merge_L5$maternal_age <- factor(tradeoff_merge_L5$maternal_age, 
                                         levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_L5)
anova(m.interaction)

# Visualization:
LT50_L5 <- ggplot(tradeoff_merge_L5, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_L5

# BmanRUS:

# Models for each maternal age:
tradeoff_merge_RUS_Y <- filter(tradeoff_merge_RUS, maternal_age == "young")
tradeoff_merge_RUS_O <- filter(tradeoff_merge_RUS, maternal_age == "old")

lm_50_RUS_Y <- lm(time_fifty ~ surv_time, data = tradeoff_merge_RUS_Y)
summary(lm_50_RUS_Y)
lm_50_RUS_O <- lm(time_fifty ~ surv_time, data = tradeoff_merge_RUS_O)
summary(lm_50_RUS_O)

# Statistics (ANCOVA):
tradeoff_merge_RUS$maternal_age <- factor(tradeoff_merge_RUS$maternal_age, 
                                          levels = c("young", "old"), labels = c("Y","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_RUS)
anova(m.interaction)

# Visualization:
LT50_RUS <- ggplot(tradeoff_merge_RUS, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_RUS

# BmanRUS-RE:

# Models for each maternal age:
tradeoff_merge_RUS_RE_Y <- filter(tradeoff_merge_RUS_RE, maternal_age == "young")
tradeoff_merge_RUS_RE_M <- filter(tradeoff_merge_RUS_RE, maternal_age == "mid")
tradeoff_merge_RUS_RE_O <- filter(tradeoff_merge_RUS_RE, maternal_age == "old")

lm_50_RUS_RE_Y <- lm(time_fifty ~ surv_time, data = tradeoff_merge_RUS_RE_Y)
summary(lm_50_RUS_RE_Y)
lm_50_RUS_RE_M <- lm(time_fifty ~ surv_time, data = tradeoff_merge_RUS_RE_M)
summary(lm_50_RUS_RE_M)
lm_50_RUS_RE_O <- lm(time_fifty ~ surv_time, data = tradeoff_merge_RUS_RE_O)
summary(lm_50_RUS_RE_O)

# Statistics (ANCOVA):
tradeoff_merge_RUS_RE$maternal_age <- factor(tradeoff_merge_RUS_RE$maternal_age, 
                                             levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_RUS_RE)
anova(m.interaction)

# Visualization:
LT50_RUS_RE <- ggplot(tradeoff_merge_RUS_RE, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_RUS_RE

# BpL1:

# Models for each maternal age:
tradeoff_merge_BpL1_Y <- filter(tradeoff_merge_BpL1, maternal_age == "young")
tradeoff_merge_BpL1_M <- filter(tradeoff_merge_BpL1, maternal_age == "mid")
tradeoff_merge_BpL1_O <- filter(tradeoff_merge_BpL1, maternal_age == "old")

lm_50_BpL1_Y <- lm(time_fifty ~ surv_time, data = tradeoff_merge_BpL1_Y)
summary(lm_50_BpL1_Y)
lm_50_BpL1_M <- lm(time_fifty ~ surv_time, data = tradeoff_merge_BpL1_M)
summary(lm_50_BpL1_M)
lm_50_BpL1_O <- lm(time_fifty ~ surv_time, data = tradeoff_merge_BpL1_O)
summary(lm_50_BpL1_O)

# Statistics (ANCOVA):
tradeoff_merge_BpL1$maternal_age <- factor(tradeoff_merge_BpL1$maternal_age, 
                                           levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_BpL1)
anova(m.interaction)

# Visualization:
LT50_BpL1 <- ggplot(tradeoff_merge_BpL1, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_BpL1

# Exploring tradeoffs: MDR versus lifespan: ---------------------------------------------------

# Making a dataframe with lifespan and MDR data:
tradeoff_MDR <- cbind(survival_data_per_ind, LRO_per_ind$max_RO)
tradeoff_MDR <- rename(tradeoff_MDR, "MDR" = "...7")
tradeoff_MDR_F1 <- filter(tradeoff_MDR, generation != "F0")

tradeoff_MDR_F1$strain <- factor(tradeoff_MDR_F1$strain, 
                                   levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

# Exploring relationships among maternal age cohorts, within each strain:
tradeoff_MDR_L5 <- filter(tradeoff_MDR_F1, strain == "BmanL5")
tradeoff_MDR_RUS <- filter(tradeoff_MDR_F1, strain == "BmanRUS")
tradeoff_MDR_RUS_RE <- filter(tradeoff_MDR_F1, strain == "BmanRUS-RE")
tradeoff_MDR_BpL1 <- filter(tradeoff_MDR_F1, strain == "BpL1")

# BmanL5:

# Models for each maternal age:
tradeoff_MDR_L5_Y <- filter(tradeoff_MDR_L5, maternal_age == "young")
tradeoff_MDR_L5_M <- filter(tradeoff_MDR_L5, maternal_age == "mid")
tradeoff_MDR_L5_O <- filter(tradeoff_MDR_L5, maternal_age == "old")

lm_MDR_L5_Y <- lm(MDR ~ surv_time, data = tradeoff_MDR_L5_Y)
summary(lm_MDR_L5_Y)
lm_MDR_L5_M <- lm(MDR ~ surv_time, data = tradeoff_MDR_L5_M)
summary(lm_MDR_L5_M)
lm_MDR_L5_O <- lm(MDR ~ surv_time, data = tradeoff_MDR_L5_O)
summary(lm_MDR_L5_O)

# Statistics (ANCOVA):
tradeoff_MDR_L5$maternal_age <- factor(tradeoff_MDR_L5$maternal_age, 
                                         levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(MDR ~ surv_time*maternal_age, data = tradeoff_MDR_L5)
anova(m.interaction)

# Visualization:
MDR_L5 <- ggplot(tradeoff_MDR_L5, aes(x=surv_time, y = MDR, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
MDR_L5

# BmanRUS:

# Models for each maternal age:
tradeoff_MDR_RUS_Y <- filter(tradeoff_MDR_RUS, maternal_age == "young")
tradeoff_MDR_RUS_O <- filter(tradeoff_MDR_RUS, maternal_age == "old")

lm_MDR_RUS_Y <- lm(MDR ~ surv_time, data = tradeoff_MDR_RUS_Y)
summary(lm_MDR_RUS_Y)
lm_MDR_RUS_O <- lm(MDR ~ surv_time, data = tradeoff_MDR_RUS_O)
summary(lm_MDR_RUS_O)

# Statistics (ANCOVA):
tradeoff_MDR_RUS$maternal_age <- factor(tradeoff_MDR_RUS$maternal_age, 
                                       levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(MDR ~ surv_time*maternal_age, data = tradeoff_MDR_RUS)
anova(m.interaction)

# Visualization:
MDR_RUS <- ggplot(tradeoff_MDR_RUS, aes(x=surv_time, y = MDR, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
MDR_RUS

# BmanRUS-RE:

# Models for each maternal age:
tradeoff_MDR_RUS_RE_Y <- filter(tradeoff_MDR_RUS_RE, maternal_age == "young")
tradeoff_MDR_RUS_RE_M <- filter(tradeoff_MDR_RUS_RE, maternal_age == "mid")
tradeoff_MDR_RUS_RE_O <- filter(tradeoff_MDR_RUS_RE, maternal_age == "old")

lm_MDR_RUS_RE_Y <- lm(MDR ~ surv_time, data = tradeoff_MDR_RUS_RE_Y)
summary(lm_MDR_RUS_RE_Y)
lm_MDR_RUS_RE_M <- lm(MDR ~ surv_time, data = tradeoff_MDR_RUS_RE_M)
summary(lm_MDR_RUS_RE_M)
lm_MDR_RUS_RE_O <- lm(MDR ~ surv_time, data = tradeoff_MDR_RUS_RE_O)
summary(lm_MDR_RUS_RE_O)

# Statistics (ANCOVA):
tradeoff_MDR_RUS_RE$maternal_age <- factor(tradeoff_MDR_RUS_RE$maternal_age, 
                                       levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(MDR ~ surv_time*maternal_age, data = tradeoff_MDR_RUS_RE)
anova(m.interaction)

# Visualization:
MDR_RUS_RE <- ggplot(tradeoff_MDR_RUS_RE, aes(x=surv_time, y = MDR, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
MDR_RUS_RE

# BpL1:

# Models for each maternal age:
tradeoff_MDR_BpL1_Y <- filter(tradeoff_MDR_BpL1, maternal_age == "young")
tradeoff_MDR_BpL1_M <- filter(tradeoff_MDR_BpL1, maternal_age == "mid")
tradeoff_MDR_BpL1_O <- filter(tradeoff_MDR_BpL1, maternal_age == "old")

lm_MDR_BpL1_Y <- lm(MDR ~ surv_time, data = tradeoff_MDR_BpL1_Y)
summary(lm_MDR_BpL1_Y)
lm_MDR_BpL1_M <- lm(MDR ~ surv_time, data = tradeoff_MDR_BpL1_M)
summary(lm_MDR_BpL1_M)
lm_MDR_BpL1_O <- lm(MDR ~ surv_time, data = tradeoff_MDR_BpL1_O)
summary(lm_MDR_BpL1_O)

# Statistics (ANCOVA):
tradeoff_MDR_BpL1$maternal_age <- factor(tradeoff_MDR_BpL1$maternal_age, 
                                       levels = c("young", "mid", "old"), labels = c("Y", "M","O"))

m.interaction <- lm(MDR ~ surv_time*maternal_age, data = tradeoff_MDR_BpL1)
anova(m.interaction)

# Visualization:
MDR_BpL1 <- ggplot(tradeoff_MDR_BpL1, aes(x=surv_time, y = MDR, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=24), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
MDR_BpL1

# Reproducing main text trade-off figure:--------------------------------------------------

LRO <- ggarrange(LRO_BpL1, LRO_RUS,
                 LRO_RUS_RE, LRO_L5,
                 ncol = 1, nrow = 4,
                 labels = c("A","E","I","M"),
                 font.label = list(size = 20),
                 common.legend = TRUE, legend = "top")

LRO <- annotate_figure(LRO, left = textGrob(expression(paste("Lifetime reproductive output (offspring ind"^"-1"*")")),
                                            rot = 90, gp = gpar(fontsize = 24)))
LRO

RPP <- ggarrange(RPP_BpL1, RPP_RUS,
                 RPP_RUS_RE, RPP_L5,
                 ncol = 1, nrow = 4,
                 labels = c("B","F","J","N"),
                 font.label = list(size = 20),
                 common.legend = TRUE, legend = "top")

RPP <- annotate_figure(RPP, left = textGrob("Reproductive period (% of lifespan)",
                                            rot = 90, gp = gpar(fontsize = 24)))
RPP

MDR <- ggarrange(MDR_BpL1, MDR_RUS,
                 MDR_RUS_RE, MDR_L5,
                 ncol = 1, nrow = 4,
                 labels = c("C","G","K","O"),
                 font.label = list(size = 20),
                 common.legend = TRUE, legend = "top")

MDR <- annotate_figure(MDR, left = textGrob("Maximum daily reproduction", 
                                            rot = 90, gp = gpar(fontsize = 24)))
MDR

LT50 <- ggarrange(LT50_BpL1, LT50_RUS,
                  LT50_RUS_RE, LT50_L5,
                  ncol = 1, nrow = 4,
                  labels = c("D","H","L","P"),
                  font.label = list(size = 20),
                  common.legend = TRUE, legend = "top")

LT50 <- annotate_figure(LT50, left = textGrob("Age of 50% LRO (d)", 
                                              rot = 90, gp = gpar(fontsize = 24)))
LT50

tiff("tradeoff.tiff", units="in", width=16, height=12, res=300)

BIG_BOI <- ggarrange(LRO, NULL, RPP, NULL, MDR, NULL, LT50,
                     ncol = 7, nrow = 1, widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1))
annotate_figure(BIG_BOI, bottom = textGrob("Lifespan (d)", gp = gpar(fontsize = 28)))

dev.off()
