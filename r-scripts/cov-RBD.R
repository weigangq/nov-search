library(tidyverse)
#library(plotly)
#library(ggplot2)
#library(colorspace)

final_variant_scores <- read.csv('https://media.githubusercontent.com/media/jbloomlab/SARS-CoV-2-RBD_DMS_variants/main/results/final_variant_scores/final_variant_scores.csv')
final_variant_scores.wide <- final_variant_scores %>% 
  select(target, position, mutant, bind) %>% 
  pivot_wider(names_from = "mutant", values_from = "bind")

#####################################
# Qiu code using group_by & join
x <- final_variant_scores %>% 
  select(target, position, mutant, bind) %>% 
  mutate(bind2 = exp(-log(bind)))

x.sum <- x %>% 
  group_by(target, position) %>% 
  summarise(sum.bind = sum(bind2))

x1 <- x %>% left_join(x.sum, c("target", "position"))

x1 <- x1 %>% mutate(prob = bind2/sum.bind)

x1 %>% ggplot(aes(x= position, y= prob, color = target)) + 
  geom_point() + 
  facet_wrap(~mutant) + 
  labs(title = "Probability Analog For Mutation Affiny At Each Site", caption = "Scatterplot of the probability of affinity (dissociation constant KD) of each mutant at each site") +
  theme_bw()

#Plot the KD probaility for each mutation at each position
x1 %>% ggplot(aes(x= mutant, y= prob, color = target)) + 
  geom_point(shape = 1) + 
  facet_wrap(~position) + 
#  labs(title = "Probability Analog For Mutation Affiny At Each Site", caption = "Scatterplot of the probability of affinity (dissociation constant KD) of each mutant at each site") +
  theme_bw()

#Get top 30 most epistatic sites
x1 %>% 
  group_by(position, mutant) %>% 
  summarise(sd = sd(prob)) %>% 
  filter(!is.na(sd))

x1 %>% 
  group_by(position, mutant) %>% 
  summarise(sd = sd(prob)) %>% 
  filter(!is.na(sd)) %>% 
  filter(position == 338) %>% 
  arrange(sd)

#Plot only position 338
x1 %>% 
  filter(position == 338) %>% 
  ggplot(aes(x= mutant, y= prob, color = target)) + 
  geom_point() + 
  #  facet_wrap(~position) + 
  labs(title = "Probability Analog For Mutation Affiny At Each Site", caption = "Scatterplot of the probability of affinity (dissociation constant KD) of each mutant at each site") +
  theme_bw()

x1 %>% 
  group_by(position, mutant) %>% 
  summarise(sd = sd(prob)) %>% 
  filter(!is.na(sd)) %>% 
  filter(position == 338) %>% 
  summarise(max=max(sd))

#print max sd for each position
x1 %>% 
  group_by(position, mutant) %>% 
  summarise(sd = sd(prob)) %>% 
  filter(!is.na(sd)) %>% 
  group_by(position) %>% 
  summarise(max=max(sd))

x1 %>% 
  group_by(position, mutant) %>% 
  summarise(sd = sd(prob)) %>% 
  filter(!is.na(sd)) %>% 
  group_by(position) %>% 
  mutate(diff = sd - max(sd)) %>%  
  head()

#print only max sd for each position
x1 %>% 
  group_by(position, mutant) %>% 
  summarise(sd = sd(prob)) %>% 
  filter(!is.na(sd)) %>% 
  group_by(position) %>% 
  filter(sd == max(sd)) %>%  
  head()

# get top shifts
x.sd <- x1 %>% 
  group_by(position, mutant) %>% 
  summarise(sd = sd(prob)) %>% #  variability among 5 targets
  filter(!is.na(sd)) %>% 
  group_by(position) %>% 
  filter(sd == max(sd)) # pick the top variable AA for each position

x.sd %>% 
  arrange(sd) %>% 
  ggplot(aes(x=position, y= sd)) + geom_bar(stat = "identity")
x.sd %>% 
  arrange(sd) %>% 
  pull(position)

pos <- x.sd %>% 
  arrange(sd) %>% 
  pull(position)

x.sd$position <- factor(x.sd$position, levels = pos)
x.sd <- x.sd %>% mutate(top = if_else(sd > 0.005, TRUE, FALSE))
x.sd %>% 
  filter(top == TRUE) %>% 
  arrange(sd) %>% 
  ggplot(aes(x=position, y= sd, label = paste(position, mutant))) + geom_bar(stat = "identity") +
  theme_bw() + 
  coord_flip() + 
  geom_text(nudge_y = 0.001)

x.sd %>% 
  arrange(sd) %>% 
  ggplot(aes(sd)) + 
  geom_boxplot()

x.sd2 <- x.sd %>% 
  arrange(sd)

x.sd2 %>% filter(sd > 5e-3)
x.mad <- x1 %>% 
  group_by(position, mutant) %>% 
  summarise(mad = mad(prob)) %>% 
  filter(!is.na(mad)) %>% 
  group_by(position) %>% 
  filter(mad == max(mad))

write_tsv(x.sd, "epi-shift.tsv")
#replicate fig 2c
x1 %>% 
  select(target, position, mutant, bind) %>% 
  filter(position==498 & target %in% c("Beta", "Wuhan-Hu-1")) %>%
  group_by(target) %>% 
  pivot_wider(names_from = target, values_from = bind) %>% 
  rename(Wuhan = 'Wuhan-Hu-1') %>% 
  ggplot(aes(x = Wuhan, y = Beta, label = mutant)) + 
  geom_text() + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "site 498", x = "Binding affinity Wuhan-Hu-1 RBD", y = "Binding affinity Beta RBD") + theme_bw()

x1 %>% 
  select(target, position, mutant, bind) %>% 
  filter(position==498 & target %in% c("Beta", "Wuhan-Hu-1")) %>%
  group_by(target) %>% pivot_wider(names_from = target, values_from = bind) %>% 
  rename(Wuhan = 'Wuhan-Hu-1') %>% 
  ggplot(aes(x = Wuhan, y = Beta, label = mutant)) + 
  geom_text() + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  labs(title = "site 498", x = "Binding affinity Wuhan-Hu-1 RBD", y = "Binding affinity Beta RBD") + 
  xlim(4,12) + 
  ylim(4, 12) +
  theme_bw() +
  theme(aspect.ratio=1)




