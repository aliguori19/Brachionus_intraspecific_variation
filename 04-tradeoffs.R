# Code for analyzing and visualizing relationships between lifespan and reproduction (trade-offs present?)

# Required packages:
library(tidyverse)
library(reshape2)
library(survival)
library(survminer)
library(car)
library(lsmeans)

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

# Calculating LRO per individual:
LRO_per_ind <- neo_long %>%
  group_by(generation, maternal_age, strain, plate, well) %>%
  summarize(LRO = sum(num_neonates, na.rm = TRUE),
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

# Exploring differences among strains, lumping across maternal age cohorts:

# Statistics:
m.interaction <- lm(LRO ~ surv_time*strain, data = tradeoff_df_F1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "strain", var="surv_time")
pairs(m.lst)

tradeoff_df_F1$strain <- factor(tradeoff_df_F1$strain, levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

# Visualization:
TO_1 <- ggplot(tradeoff_df_F1, aes(x=surv_time, y=LRO, color = strain)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2", "goldenrod2"), name = "") +
  xlab(element_blank()) +
  ylab("Lifetime reproductive output") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text=element_text(size=18),
        legend.text=element_text(size=16))
TO_1

# Exploring relationships among maternal age cohorts, within each strain:

tradeoff_df_F1_L5 <- filter(tradeoff_df_F1, strain == "BmanL5")
tradeoff_df_F1_RUS <- filter(tradeoff_df_F1, strain == "BmanRUS")
tradeoff_df_F1_RUS_RE <- filter(tradeoff_df_F1, strain == "BmanRUS-RE")
tradeoff_df_F1_BpL1 <- filter(tradeoff_df_F1, strain == "BpL1")

# BmanL5:

# Pooling maternal ages:
lm_L5 <- lm(LRO ~ surv_time, data = tradeoff_df_F1_L5)
summary(lm_L5)

# Statistics:
tradeoff_df_F1_L5$maternal_age <- factor(tradeoff_df_F1_L5$maternal_age, levels = c("3", "6", "10"), labels = c(
  "Y", "M","O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_L5)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
LRO_L5 <- ggplot(tradeoff_df_F1_L5, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  geom_label(
    label="BmanL5", 
    size = 6,
    x=7,
    y=30,
    label.size = 0,
    color = "black") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_L5

# BmanRUS:

# Pooling maternal ages:
lm_RUS <- lm(LRO ~ surv_time, data = tradeoff_df_F1_RUS)
summary(lm_RUS)

# Statistics:
tradeoff_df_F1_RUS$maternal_age <- factor(tradeoff_df_F1_RUS$maternal_age, levels = c("3", "11"), labels = c(
  "Y", "O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_RUS)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")

# Visualization:
LRO_RUS <- ggplot(tradeoff_df_F1_RUS, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  geom_label(
    label="BmanRUS", 
    size = 6,
    x=8,
    y=32,
    label.size = 0,
    color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_RUS

# BmanRUS-RE:

# Pooling maternal ages:
lm_RUS_RE <- lm(LRO ~ surv_time, data = tradeoff_df_F1_RUS_RE)
summary(lm_RUS_RE)

# Statistics:
tradeoff_df_F1_RUS_RE$maternal_age <- factor(tradeoff_df_F1_RUS_RE$maternal_age, 
                                             levels = c("3", "6", "10"), labels = c("Y", "M","O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_RUS_RE)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
LRO_RUS_RE <- ggplot(tradeoff_df_F1_RUS_RE, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  geom_label(
    label="BmanRUS-RE", 
    size = 6,
    x=10,
    y=39,
    label.size = 0,
    color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_RUS_RE

# BpL1:

# Pooling maternal ages:
lm_BpL1 <- lm(LRO ~ surv_time, data = tradeoff_df_F1_BpL1)
summary(lm_BpL1)

# Statistics:
tradeoff_df_F1_BpL1$maternal_age <- factor(tradeoff_df_F1_BpL1$maternal_age, 
                                           levels = c("3", "6", "9"), labels = c("Y", "M","O"))

m.interaction <- lm(LRO ~ surv_time*maternal_age, data = tradeoff_df_F1_BpL1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
LRO_BpL1 <- ggplot(tradeoff_df_F1_BpL1, aes(x=surv_time, y=LRO, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  geom_label(
    label="BpL1", 
    size = 6,
    x=5,
    y=25,
    label.size = 0,
    color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LRO_BpL1

# Exploring tradeoffs: Reproductive period (percentage of lifespan) versus lifespan: ---------------------------------------------------

# Making a dataframe with lifespan and reproductive period data:
tradeoff_df_rep <- cbind(survival_data_per_ind, rep_per_ind$percent_rep)
tradeoff_df_rep <- rename(tradeoff_df_rep, "rep_percent" = "...7")
tradeoff_df_rep_F1 <- filter(tradeoff_df_rep, generation != "F0")

# Exploring differences among strains, lumping across maternal age cohorts:

# Statistics:
m.interaction <- lm(rep_percent ~ surv_time*strain, data = tradeoff_df_rep_F1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "strain", var="surv_time")
pairs(m.lst)

tradeoff_df_rep_F1$strain <- factor(tradeoff_df_rep_F1$strain, levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

# Visualization:
TO_2 <- ggplot(tradeoff_df_rep_F1, aes(x=surv_time, y=rep_percent, color = strain)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2", "goldenrod2"), name = "") +
  xlab(element_blank()) +
  ylab("Reproductive period (% lifespan)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text=element_text(size=18),
        legend.text = element_text(size=16))
TO_2

# Exploring relationships among maternal age cohorts, within each strain:

tradeoff_rep_L5 <- filter(tradeoff_df_rep_F1, strain == "BmanL5")
tradeoff_rep_RUS <- filter(tradeoff_df_rep_F1, strain == "BmanRUS")
tradeoff_rep_RUS_RE <- filter(tradeoff_df_rep_F1, strain == "BmanRUS-RE")
tradeoff_rep_BpL1 <- filter(tradeoff_df_rep_F1, strain == "BpL1")

# BmanL5:

# Statistics:
tradeoff_rep_L5$maternal_age <- factor(tradeoff_rep_L5$maternal_age, 
                                       levels = c("3", "6", "10"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_L5)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
RPP_L5 <- ggplot(tradeoff_rep_L5, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_L5

# BmanRUS:

# Statistics:
tradeoff_rep_RUS$maternal_age <- factor(tradeoff_rep_RUS$maternal_age, 
                                        levels = c("3", "11"), labels = c("Y","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_RUS)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")

# Visualization:
RPP_RUS <- ggplot(tradeoff_rep_RUS, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_RUS

# BmanRUS-RE:

# Statistics:
tradeoff_rep_RUS_RE$maternal_age <- factor(tradeoff_rep_RUS_RE$maternal_age, 
                                           levels = c("3", "6", "10"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_RUS_RE)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
RPP_RUS_RE <- ggplot(tradeoff_rep_RUS_RE, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_RUS_RE

# BpL1:

# Statistics:
tradeoff_rep_BpL1$maternal_age <- factor(tradeoff_rep_BpL1$maternal_age, 
                                         levels = c("3", "6", "9"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_percent ~ surv_time*maternal_age, data = tradeoff_rep_BpL1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
RPP_BpL1 <- ggplot(tradeoff_rep_BpL1, aes(x=surv_time, y=rep_percent, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPP_BpL1

# Exploring tradeoffs: Reproductive period (days) versus lifespan: ---------------------------------------------------

# Making a dataframe with lifespan and reproductive period data:
tradeoff_df_repd <- cbind(survival_data_per_ind, rep_per_ind$rep_time)
tradeoff_df_repd <-  rename(tradeoff_df_repd, "rep_time" = "...7")
tradeoff_df_repd_F1 <- filter(tradeoff_df_repd, generation != "F0")

# Exploring differences among strains, lumping across maternal age cohorts:

# Statistics:
m.interaction <- lm(rep_time ~ surv_time*strain, data = tradeoff_df_repd_F1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "strain", var="surv_time")
pairs(m.lst)

tradeoff_df_repd_F1$strain <- factor(tradeoff_df_repd_F1$strain, levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

TO_3 <- ggplot(tradeoff_df_repd_F1, aes(x=surv_time, y=rep_time, color = strain)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2", "goldenrod2"), name = "") +
  xlab(element_blank()) +
  ylab("Reproductive period (d)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text=element_text(size=18),
        legend.text = element_text(size=16))
TO_3

# Exploring relationships among maternal age cohorts, within each strain:

tradeoff_repd_L5 <- filter(tradeoff_df_repd_F1, strain == "BmanL5")
tradeoff_repd_RUS <- filter(tradeoff_df_repd_F1, strain == "BmanRUS")
tradeoff_repd_RUS_RE <- filter(tradeoff_df_repd_F1, strain == "BmanRUS-RE")
tradeoff_repd_BpL1 <- filter(tradeoff_df_repd_F1, strain == "BpL1")

# BmanL5:

# Pooling maternal ages:
lm_L5 <- lm(rep_time ~ surv_time, data = tradeoff_repd_L5)
summary(lm_L5)

# Statistics:
tradeoff_repd_L5$maternal_age <- factor(tradeoff_repd_L5$maternal_age, 
                                        levels = c("3", "6", "10"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_time ~ surv_time*maternal_age, data = tradeoff_repd_L5)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
RPD_L5 <- ggplot(tradeoff_repd_L5, aes(x=surv_time, y=rep_time, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPD_L5

# BmanRUS:

# Pooling maternal ages:
lm_RUS <- lm(rep_time ~ surv_time, data = tradeoff_repd_RUS)
summary(lm_RUS)

# Statistics:
tradeoff_repd_RUS$maternal_age <- factor(tradeoff_repd_RUS$maternal_age, 
                                         levels = c("3","11"), labels = c("Y", "O"))

m.interaction <- lm(rep_time ~ surv_time*maternal_age, data = tradeoff_repd_RUS)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")

# Visualization:
RPD_RUS <- ggplot(tradeoff_repd_RUS, aes(x=surv_time, y=rep_time, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPD_RUS

# BmanRUS-RE:

# Pooling maternal ages:
lm_RUS_RE <- lm(rep_time ~ surv_time, data = tradeoff_repd_RUS_RE)
summary(lm_RUS_RE)

# Statistics:
tradeoff_repd_RUS_RE$maternal_age <- factor(tradeoff_repd_RUS_RE$maternal_age, 
                                            levels = c("3", "6", "10"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_time ~ surv_time*maternal_age, data = tradeoff_repd_RUS_RE)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
RPD_RUS_RE <- ggplot(tradeoff_repd_RUS_RE, aes(x=surv_time, y=rep_time, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPD_RUS_RE

# BpL1:

# Pooling maternal ages:
lm_BpL1 <- lm(rep_time ~ surv_time, data = tradeoff_repd_BpL1)
summary(lm_BpL1)

# Statistics:
tradeoff_repd_BpL1$maternal_age <- factor(tradeoff_repd_BpL1$maternal_age, 
                                          levels = c("3", "6", "9"), labels = c("Y", "M","O"))

m.interaction <- lm(rep_time ~ surv_time*maternal_age, data = tradeoff_repd_BpL1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
RPD_BpL1 <- ggplot(tradeoff_repd_BpL1, aes(x=surv_time, y=rep_time, color = maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
RPD_BpL1

# Exploring tradeoffs: Timing of 50% LRO production versus lifespan: ---------------------------------------------------

# Making a dataframe with lifespan and 50% LRO data:
tradeoff_merge <- merge(survival_data_per_ind, sub_summ, by = c("generation","maternal_age","strain","plate","well"))
tradeoff_merge_F1 <- filter(tradeoff_merge, generation != "F0")

# Exploring differences among strains, lumping across maternal age cohorts:

# Statistics:
m.interaction <- lm(time_fifty ~ surv_time*strain, data = tradeoff_merge_F1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "strain", var="surv_time")
pairs(m.lst)

tradeoff_merge_F1$strain <- factor(tradeoff_merge_F1$strain, 
                                   levels = c("BpL1", "L5", "RUS", "RUS-RE"), labels = c("BpL1", "BmanL5", "BmanRUS", "BmanRUS-RE"))

# Visualization:
TO_4 <- ggplot(tradeoff_merge_F1, aes(x=surv_time, y=time_fifty, color = strain)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2", "goldenrod2"), name = "") +
  xlab(element_blank()) +
  ylab("Time of 50% LRO") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text=element_text(size=18),
        legend.text = element_text(size = 16))
TO_4

# Exploring relationships among maternal age cohorts, within each strain:

tradeoff_merge_L5 <- filter(tradeoff_merge_F1, strain == "BmanL5")
tradeoff_merge_RUS <- filter(tradeoff_merge_F1, strain == "BmanRUS")
tradeoff_merge_RUS_RE <- filter(tradeoff_merge_F1, strain == "BmanRUS-RE")
tradeoff_merge_BpL1 <- filter(tradeoff_merge_F1, strain == "BpL1")

# BmanL5:

# Pooling maternal ages:
lm_L5 <- lm(time_fifty ~ surv_time, data = tradeoff_merge_L5)
summary(lm_L5)

# Statistics:
tradeoff_merge_L5$maternal_age <- factor(tradeoff_merge_L5$maternal_age, 
                                         levels = c("3", "6", "10"), labels = c("Y", "M","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_L5)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
LT50_L5 <- ggplot(tradeoff_merge_L5, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_L5

# BmanRUS:

# Pooling maternal ages:
lm_RUS <- lm(time_fifty ~ surv_time, data = tradeoff_merge_RUS)
summary(lm_RUS)

# Statistics:
tradeoff_merge_RUS$maternal_age <- factor(tradeoff_merge_RUS$maternal_age, 
                                          levels = c("3", "11"), labels = c("Y","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_RUS)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")

# Visualization:
LT50_RUS <- ggplot(tradeoff_merge_RUS, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#0072B2"), labels = c("Y","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_RUS

# BmanRUS-RE:

# Pooling maternal ages:
lm_RUS_RE <- lm(time_fifty ~ surv_time, data = tradeoff_merge_RUS_RE)
summary(lm_RUS_RE)

# Statistics:
tradeoff_merge_RUS_RE$maternal_age <- factor(tradeoff_merge_RUS_RE$maternal_age, 
                                             levels = c("3", "6", "10"), labels = c("Y", "M","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_RUS_RE)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
LT50_RUS_RE <- ggplot(tradeoff_merge_RUS_RE, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_RUS_RE

# BpL1:

# Pooling maternal ages:
lm_BpL1 <- lm(time_fifty ~ surv_time, data = tradeoff_merge_BpL1)
summary(lm_BpL1)

# Statistics:
tradeoff_merge_BpL1$maternal_age <- factor(tradeoff_merge_BpL1$maternal_age, 
                                           levels = c("3", "6", "9"), labels = c("Y", "M","O"))

m.interaction <- lm(time_fifty ~ surv_time*maternal_age, data = tradeoff_merge_BpL1)
anova(m.interaction)

m.interaction$coefficients
m.lst <- lstrends(m.interaction, "maternal_age", var="surv_time")
pairs(m.lst)

# Visualization:
LT50_BpL1 <- ggplot(tradeoff_merge_BpL1, aes(x=surv_time, y = time_fifty, color=maternal_age)) +
  geom_jitter(size=2) +
  geom_smooth(method="lm", se= FALSE) +
  scale_color_manual(values = c("black", "#999999", "#0072B2"), labels = c("Y","M","O"), name="") +
  xlab(element_blank()) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        axis.title=element_text(size=22, face="bold"), 
        axis.text=element_text(size=20),
        legend.text=element_text(size=16))
LT50_BpL1

# Reproducing supplemental trade-off figure:---------------------------------------------

tiff("tradeoff_supp.tiff", units="in", width=14, height=12, res=1200)
tradeoff <- ggarrange(TO_1, TO_2,
                      TO_3, TO_4,
                      ncol = 2, nrow = 2,
                      labels = c("A","B","C","D"),
                      font.label = list(size = 20),
                      common.legend = TRUE, legend = "top")
annotate_figure(tradeoff, bottom = textGrob("Lifespan (d)", gp = gpar(fontsize = 20, fontface = 'bold')))
dev.off()

# Reproducing main text trade-off figure:--------------------------------------------------

LRO <- ggarrange(LRO_BpL1, LRO_RUS,
                 LRO_RUS_RE, LRO_L5,
                 ncol = 1, nrow = 4,
                 labels = c("A","E","I","M"),
                 font.label = list(size = 20),
                 common.legend = TRUE, legend = "top")

LRO <- annotate_figure(LRO, left = textGrob("Lifetime reproductive output", rot = 90, gp = gpar(fontsize = 24, fontface = 'bold')))
LRO

RPP <- ggarrange(RPP_BpL1, RPP_RUS,
                 RPP_RUS_RE, RPP_L5,
                 ncol = 1, nrow = 4,
                 labels = c("B","F","J","N"),
                 font.label = list(size = 20),
                 common.legend = TRUE, legend = "top")

RPP <- annotate_figure(RPP, left = textGrob("Reproductive period (% of lifespan)", rot = 90, gp = gpar(fontsize = 24, fontface = 'bold')))
RPP

RPD <- ggarrange(RPD_BpL1, RPD_RUS,
                 RPD_RUS_RE, RPD_L5,
                 ncol = 1, nrow = 4,
                 labels = c("C","G","K","O"),
                 font.label = list(size = 20),
                 common.legend = TRUE, legend = "top")

RPD <- annotate_figure(RPD, left = textGrob("Reproductive period (d)", rot = 90, gp = gpar(fontsize = 24, fontface = 'bold')))
RPD

LT50 <- ggarrange(LT50_BpL1, LT50_RUS,
                  LT50_RUS_RE, LT50_L5,
                  ncol = 1, nrow = 4,
                  labels = c("D","H","L","P"),
                  font.label = list(size = 20),
                  common.legend = TRUE, legend = "top")

LT50 <- annotate_figure(LT50, left = textGrob("Time of 50% LRO", rot = 90, gp = gpar(fontsize = 24, fontface = 'bold')))
LT50

tiff("tradeoff.tiff", units="in", width=16, height=12, res=1200)

BIG_BOI <- ggarrange(LRO, NULL, RPP, NULL, RPD, NULL, LT50,
                     ncol = 7, nrow = 1, widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1))
annotate_figure(BIG_BOI, bottom = textGrob("Lifespan (d)", gp = gpar(fontsize = 28, fontface = 'bold')))

dev.off()