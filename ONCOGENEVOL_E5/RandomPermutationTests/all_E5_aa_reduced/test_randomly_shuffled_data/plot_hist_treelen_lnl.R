setwd("/DATA/work/ONCOGENEVOL/E5/E5_Original_files/PermutationTests/all_E5_aa_new_reduced/test_randomly_shuffled_data");

library(reshape2); #to melt data
library(ggplot2); #for plotting
library(pgirmess); #for kruskalmc
library(car); #for LeveneTest

#data <- read.table("table_treelen_lnl_baliphy_100iter.txt", header=T);
data <- read.table("table_treelen_lnl_phyml_100iter.txt", header=T);

tlen_H0 <- data$treelen_ALL;
tlen_H1 <- data$treelen_red + data$treelen_blue;
treelen <- cbind(tlen_H0, tlen_H1);
treelen_melt <- melt(treelen);

lnl_H0 <- data$lnl_ALL;
lnl_H1 <- data$lnl_red + data$lnl_blue;
loglik <- cbind(lnl_H0, lnl_H1);
loglik_melt <- melt(loglik);

bf_H0 <- lnl_H0 - lnl_H0;
bf_H1 <- lnl_H0 - lnl_H1;
bayes <- cbind(bf_H0, bf_H1);
bayes_melt <- melt(bayes);
bayes2 <- bf_H1;
bayes_melt2 <- melt(bayes2);


svg(file="phyml_hist_MLtreelen_100iter.svg");
ggplot(treelen_melt, aes(x=value, fill=Var2)) +
  geom_histogram(bins=60) +
  xlab("ML tree length under CA") +
  ggtitle("red-blue E5 aa reduced");
dev.off();

svg(file="phyml_density_MLtreelen_100iter.svg");
ggplot(treelen_melt, aes(x=value, fill=Var2)) +
  geom_density(alpha=0.5) +
  xlab("ML tree length under CA") +
  ggtitle("red-blue E5 aa reduced");
dev.off();

svg(file="phyml_hist_LnL_100iter.svg");
ggplot(loglik_melt, aes(x=value, fill=Var2)) +
  geom_histogram(bins=30) +
  xlab("Log likelihood") +
  ggtitle("red-blue E5 aa reduced");
dev.off();

svg(file="phyml_density_LnL_100iter.svg");
ggplot(loglik_melt, aes(x=value, fill=Var2)) +
  geom_density(alpha=0.5) +
  xlab("Log likelihood") +
  ggtitle("red-blue E5 aa reduced");
dev.off();

svg(file="phyml_hist_BF_100iter.svg");
ggplot(bayes_melt, aes(x=value, fill=Var2)) +
  geom_histogram(bins=30) +
  xlab("Bayes Factor") +
  ggtitle("red-blue E5 aa reduced");
dev.off();

svg(file="phyml_density_BF_100iter.svg");
ggplot(bayes_melt2, aes(x=value)) +
  geom_density(alpha=0.5) +
  xlab("Bayes Factor") +
  ggtitle("red-blue E5 aa reduced");
dev.off();

##################### STATS ##################### 
##Non-parametric alternative to one-way ANOVA test
kruskal.test(treelen_melt$value ~ treelen_melt$Var2);
kruskalmc(treelen_melt$value ~ treelen_melt$Var2);
wilcox.test(tlen_H0, tlen_H1);


kruskal.test(loglik_melt$value ~ loglik_melt$Var2);
kruskalmc(loglik_melt$value ~ loglik_melt$Var2);
wilcox.test(lnl_H0, lnl_H1);


quit();
