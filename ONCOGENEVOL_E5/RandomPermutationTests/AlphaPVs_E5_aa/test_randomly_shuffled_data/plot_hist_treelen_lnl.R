setwd("/DATA/work/ONCOGENEVOL/E5/E5_Original_files/PermutationTests/AlphaPVs_E5_aa_new/test_randomly_shuffled_data");

library(reshape2); #to melt data
library(ggplot2); #for plotting
library(pgirmess); #for kruskalmc
library(car); #for LeveneTest

#data <- read.table("table_treelen_lnl_baliphy_100iter.txt", header=T);
data <- read.table("table_treelen_lnl_phyml_100iter.txt", header=T);

tlen_H0 <- data$treelen_ALL;
tlen_H1 <- data$treelen_part_H4 + data$treelen_beta;
tlen_H2 <- data$treelen_part_H3 + data$treelen_beta + data$treelen_epsilon_dseta;
tlen_H3 <- data$treelen_part_H5 + data$treelen_part_H2 + data$treelen_beta;
tlen_H4 <- data$treelen_part_H5 + data$treelen_beta + data$treelen_gammadelta + data$treelen_delta;
tlen_H5 <- data$treelen_part_H1 + data$treelen_part_H2 + data$treelen_beta + data$treelen_epsilon_dseta;
tlen_H6 <- data$treelen_alpha1 + data$treelen_alpha2 + data$treelen_gammadelta + data$treelen_delta + data$treelen_beta + data$treelen_epsilon_dseta;
treelen <- cbind(tlen_H0, tlen_H1, tlen_H2, tlen_H3, tlen_H4, tlen_H5, tlen_H6);
treelen_melt <- melt(treelen);


lnl_H0 <- data$lnl_ALL;
lnl_H1 <- data$lnl_part_H4 + data$lnl_beta;
lnl_H2 <- data$lnl_part_H3 + data$lnl_beta + data$lnl_epsilon_dseta;
lnl_H3 <- data$lnl_part_H5 + data$lnl_part_H2 + data$lnl_beta;
lnl_H4 <- data$lnl_part_H5 + data$lnl_beta + data$lnl_gammadelta + data$lnl_delta;
lnl_H5 <- data$lnl_part_H1 + data$lnl_part_H2 + data$lnl_beta + data$lnl_epsilon_dseta;
lnl_H6 <- data$lnl_alpha1 + data$lnl_alpha2 + data$lnl_gammadelta + data$lnl_delta + data$lnl_beta + data$lnl_epsilon_dseta;
loglik <- cbind(lnl_H0, lnl_H1, lnl_H2, lnl_H3, lnl_H4, lnl_H5, lnl_H6);
loglik_melt <- melt(loglik);

bf_H0 <- lnl_H0 - lnl_H0;
bf_H1 <- lnl_H0 - lnl_H1;
bf_H2 <- lnl_H0 - lnl_H2;
bf_H3 <- lnl_H0 - lnl_H3; 
bf_H4 <- lnl_H0 - lnl_H4; 
bf_H5 <- lnl_H0 - lnl_H5; 
bf_H6 <- lnl_H0 - lnl_H6; 
bayes <- cbind(bf_H0, bf_H1, bf_H2, bf_H3, bf_H4, bf_H5, bf_H6);
bayes_melt <- melt(bayes);
bayes2 <- cbind(bf_H1, bf_H2, bf_H3, bf_H4, bf_H5, bf_H6);
bayes_melt2 <- melt(bayes2);


svg(file="phyml_hist_MLtreelen_100iter.svg");
ggplot(treelen_melt, aes(x=value, fill=Var2)) +
  geom_histogram(bins=60) +
  xlab("ML tree length under CA") +
  ggtitle("AlphaPVs E5 aa");
dev.off();

svg(file="phyml_density_MLtreelen_100iter.svg");
ggplot(treelen_melt, aes(x=value, fill=Var2)) +
  geom_density(alpha=0.5) +
  xlab("ML tree length under CA") +
  ggtitle("AlphaPVs E5 aa");
dev.off();

svg(file="phyml_hist_LnL_100iter.svg");
ggplot(loglik_melt, aes(x=value, fill=Var2)) +
  geom_histogram(bins=30) +
  xlab("Log likelihood") +
  ggtitle("AlphaPVs E5 aa");
dev.off();

svg(file="phyml_density_LnL_100iter.svg");
ggplot(loglik_melt, aes(x=value, fill=Var2)) +
  geom_density(alpha=0.5) +
  xlab("Log likelihood") +
  ggtitle("AlphaPVs E5 aa");
dev.off();

svg(file="phyml_hist_BF_100iter.svg");
ggplot(bayes_melt, aes(x=value, fill=Var2)) +
  geom_histogram(bins=30) +
  xlab("Bayes Factor") +
  ggtitle("AlphaPVs E5 aa");
dev.off();

svg(file="phyml_density_BF_100iter.svg");
ggplot(bayes_melt2, aes(x=value, fill=Var2)) +
  geom_density(alpha=0.5) +
  xlab("Bayes Factor") +
  ggtitle("AlphaPVs E5 aa");
dev.off();

##################### STATS #####################
### ANOVA ###
maov <- aov(treelen_melt$value ~ treelen_melt$Var2);
summary(maov);
TukeyHSD(maov);

###################################################
# 1. Homogeneity of variances
plot(maov, 1);

leveneTest(treelen_melt$value ~ treelen_melt$Var2);
### This results suggest to reject the hypothesis of Homogeneity of Variance

### Perform ANOVA test with no assumption of equal variances
oneway.test(treelen_melt$value ~ treelen_melt$Var2);

###################################################
# 2. Normality
# The normal probability plot of residuals is used to check the assumption that
# the residuals are normally distributed. It should approximately follow a 
# straight line.
plot(maov, 2);
# Extract the residuals
maov_residuals <- residuals(object = maov);
# Run Shapiro-Wilk test
shapiro.test(x = maov_residuals);
### Residuals are not normally distributed

### Non-parametric alternative tests to one-way ANOVA test
kruskal.test(treelen_melt$value ~ treelen_melt$Var2);
kruskalmc(treelen_melt$value ~ treelen_melt$Var2);
wilcox.test(tlen_H0, tlen_H1);
wilcox.test(tlen_H0, tlen_H2);
wilcox.test(tlen_H0, tlen_H3);
wilcox.test(tlen_H0, tlen_H4);
wilcox.test(tlen_H0, tlen_H5);
wilcox.test(tlen_H0, tlen_H6);


kruskal.test(loglik_melt$value ~ loglik_melt$Var2);
kruskalmc(loglik_melt$value ~ loglik_melt$Var2);
wilcox.test(lnl_H0, lnl_H1);
wilcox.test(lnl_H0, lnl_H2);
wilcox.test(lnl_H0, lnl_H3);
wilcox.test(lnl_H0, lnl_H4);
wilcox.test(lnl_H0, lnl_H5);
wilcox.test(lnl_H0, lnl_H6);

kruskal.test(bayes_melt2$value ~ bayes_melt2$Var2);
kruskalmc(bayes_melt2$value ~ bayes_melt2$Var2);
