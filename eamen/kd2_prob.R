library(tidyverse)
library(plotly)
library(ggplot2)
library(colorspace)

setwd("/Users/eamenho/Qiu_Lab/sars_epistatis/")

final_variant_scores <- read.csv("final_variant_scores.csv")
final_variant_scores.wide <- final_variant_scores %>% select(target, position, mutant, bind) %>% pivot_wider(names_from = "mutant", values_from = "bind")

#####################################
# Qiu code using group_by & join
x <- final_variant_scores %>% select(target, position, mutant, bind) %>% mutate(bind2 = exp(-log(bind)))
x.sum <- x %>% group_by(target, position) %>% summarise(sum.bind = sum(bind2))
x1 <- x %>% left_join(x.sum, c("target", "position"))
x1 <- x1 %>% mutate(prob = bind2/sum.bind)
x1 %>% ggplot(aes(x= position, y= prob, color = target)) + 
  geom_point() + 
  facet_wrap(~mutant) + 
  labs(title = "Probability Analog For Mutation Affiny At Each Site", caption = "Scatterplot of the probability of affinity (dissociation constant KD) of each mutant at each site") +
  theme_bw()

x1 %>% ggplot(aes(x= mutant, y= prob, color = target)) + 
  geom_point() + 
  facet_wrap(~position) + 
  labs(title = "Probability Analog For Mutation Affiny At Each Site", caption = "Scatterplot of the probability of affinity (dissociation constant KD) of each mutant at each site") +
  theme_bw()

x1 %>% select(target, position, mutant, bind) %>% filter(position==498 & target %in% c("Beta", "Wuhan-Hu-1")) %>% group_by(target) %>% pivot_wider(names_from = target, values_from = bind) %>% rename(Wuhan = 'Wuhan-Hu-1') %>% ggplot(aes(x = Wuhan, y = Beta, label = mutant)) + geom_text() + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +labs(title = "site 498", x = "Binding affinity Wuhan-Hu-1 RBD", y = "Binding affinity Beta RBD") + theme_bw()

x1 %>% select(target, position, mutant, bind) %>% filter(position==498 & target %in% c("Beta", "Wuhan-Hu-1")) %>% group_by(target) %>% pivot_wider(names_from = target, values_from = bind) %>% rename(Wuhan = 'Wuhan-Hu-1') %>% ggplot(aes(x = Wuhan, y = Beta, label = mutant)) + geom_text() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +labs(title = "site 498", x = "Binding affinity Wuhan-Hu-1 RBD", y = "Binding affinity Beta RBD") + xlim(4,12) + ylim(4, 12) +theme_bw() +theme(aspect.ratio=1)



###########################################

#find the e^-ln for each mutation at each position
log_final_variant_scores <- exp(-log(final_variant_scores.wide[ , c(3:18)]))

#find the total sum for each position
log_final_variant_scores <- log_final_variant_scores %>% mutate(total = rowSums(across(where(is.numeric)), na.rm=TRUE))

#find the probability
prob_final_variant_scores <- log_final_variant_scores / log_final_variant_scores[,17]
prob_final_variant_scores <- cbind(position= final_variant_scores.wide$position, prob_final_variant_scores)

#pivot longer
prob_final_variant_scores.long <- prob_final_variant_scores %>% select(!total) %>% pivot_longer(!position, names_to = "mutation", values_to = "kd2prob")

# ggplot(prob_final_variant_scores.long, aes(x= position, y= kd2prob, color = mutation)) + geom_line()
ggplot(prob_final_variant_scores.long, aes(x= position, fill = mutation)) + geom_bar(position="stack") + labs(title = "Probability Analog For Mutation Affiny At Each Site")

ggplot(prob_final_variant_scores.long, aes(x= position, y= kd2prob)) + geom_line() + facet_wrap(~mutation) + labs(title = "Probability Analog For Mutation Affiny At Each Site", caption = "Scatterplot of the probability of affinity (dissociation constant KD) of each mutant at each site")


kd2prob <- function(x) {
  log_final_variant_scores <- exp(-log(x[ , c(3:18)]))
  log_final_variant_scores <- log_final_variant_scores %>% mutate(total = rowSums(across(where(is.numeric)), na.rm=TRUE))
  prob_final_variant_scores <- log_final_variant_scores / log_final_variant_scores[,17]
  prob_final_variant_scores <- cbind(position= x$position, prob_final_variant_scores)
  return(prob_final_variant_scores)
}

kd2prob(final_variant_scores.wide)
