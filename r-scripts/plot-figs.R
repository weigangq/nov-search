library(tidyverse)
setwd("C:/Users/weiga/Dropbox/nov-search")
#setwd("Dropbox/nov-search")
library(ggrepel)

##############
# GB1 LON analysis
# Basins of Attractions

# tree doesn't work. Delete
# library(ggtree)
# tree <- read.tree("../data/gb1-mid.dnd")
# basin <- read_tsv("../data/gb1-basins.tsv", col_names = F)
# top25<- basin %>% arrange(X2) %>% tail(25) 

# tips <- tibble(id = tree$tip.label,
#               basin = c(rep("sample", 10), "ref", rep("sample",16# )))

# p <- ggtree(tree) +
#  xlim(0, 0.3) +
#  geom_treescale(x=0.2, y=25) +
#  geom_tiplab(align = T)

# p %<+% tips +
#  geom_tiplab(aes(color = ref), align = T) +
#  geom_tippoint(alpha=.5) +
#  scale_color_manual(values = c("sample" = 1,  "ref" = 2)) +
#  theme(legend.position = "none") +
#  ggtitle(label = "E. faecalis ML tree by IQTREE (ref: OG1RF); 76,752 SNVs", subtitle = "version: 11/6/2022 (W.Qiu, Hunter)")

# basins
gb1 <- read_tsv("data/gb1-basins.tsv", col_names = c("pep", "basin", "fit"))

# main
gb1 %>% 
  ggplot(aes(fit, basin)) + 
  geom_point(shape = 1) +
  theme_bw() + 
  geom_text_repel(data = gb1 %>% filter(basin > 1500), aes(fit, basin, label = pep)) +
  geom_vline(xintercept = 1, linetype = 2) 
  #scale_y_log10() + 
  #scale_x_log10() +
  #geom_smooth(method = "lm")

summary(lm(data = gb1, log(basin) ~ fit))

# insert (log10 scaled)
gb1 %>% 
  ggplot(aes(fit, basin)) + 
  geom_point(shape = 1) +
  theme_bw() + 
  geom_text_repel(data = gb1 %>% filter(basin > 1500), aes(fit, basin, label = pep)) +
  geom_vline(xintercept = 1, linetype = 2) 
  scale_y_log10() + 
#scale_x_log10() +
geom_smooth(method = "lm")

# Vasins of Slippage
# gb1 <- read_tsv("Brandon/gb1-vasins.tsv", col_names = c("pep", "vasin", "fit"))

gb1 %>% 
  ggplot(aes(fit, vasin)) + 
  geom_point(shape = 1) +
  theme_bw() + 
  #geom_text_repel(data = gb1 %>% filter(basin > 1500), aes(fit, basin, label = pep)) +
  geom_vline(xintercept = 1, linetype = 2)  
  #scale_y_log10() + 
  #scale_x_log10() 
  #geom_smooth(method = "lm")

gb1 %>% ggplot(aes(vasin)) +
  geom_histogram(bins=100) + 
  theme_bw() +
  scale_y_log10()

#summary(lm(data = gb1, log(vasin) ~ fit))
###################
# NA1 search
#######################
df.na <- read_tsv("data/na-obj-nov-comb.out") 
df.na <- df.na %>% 
  mutate(alg = if_else(algorithm == 1, '1-obj', if_else(algorithm == 2, "2-nov", "3-combo")))
df.na <- df.na %>% 
  mutate(trial = str_remove(tag, "alg.-"))

counts <- df.na %>% 
  group_by(tag, algorithm) %>% count()

counts <- counts %>% 
  mutate(alg = if_else(algorithm == 1, '1-obj', if_else(algorithm == 2, "2-nov", "3-combo")))

# plot num of generations before global peak
counts %>% 
  ggplot(aes(x = alg, y = n)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(width = 0.25, shape=1, alpha = 0.5) + 
  theme_bw() +
  xlab("Search algorithms") + 
  ylab("Generations (til global peak or max run)")

# fitness changes
fit.max = max(df.na$elite_fitness)

# plot num of generations before global peak
counts %>% 
  ggplot(aes(x = alg, y = n)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(width = 0.25, shape=1, alpha = 0.5) + 
  theme_bw() +
  xlab("Search algorithms") + 
  ylab("Generations (til global peak or max run)")

# fitness changes
#fit.max = max(df.gb$elite_fitness)

df.na %>% 
  filter(trial %in% c("trial1", "trial10", "trial20")) %>% 
  filter(gen <21) %>% 
  ggplot(aes(x = gen, y = elite_fitness, group = tag, label = elite_hap)) + 
  geom_line() + 
  geom_point(shape = 21, size=2, fill="white") + 
  theme_bw() + 
  labs(title = "Influenza NA protein landscape") +
  ylab("affinity") +
  geom_hline(yintercept = fit.max, linetype = 2) +
  #facet_wrap(~alg, nrow = 3) + 
  facet_grid(rows = vars(alg), cols = vars(trial)) +
  theme(legend.position = "bottom") +
  geom_text_repel(color = "gray", size = 3) 


###################
# GB1 search
#######################
df.gb <- read_tsv("data/gb1-obj-nov-comb.out") 
df.gb <- df.gb %>% 
  mutate(alg = if_else(algorithm == 1, '1-obj', if_else(algorithm == 2, "2-nov", "3-combo")))
df.gb <- df.gb %>% 
  mutate(trial = str_remove(tag, "alg.-"))

counts <- df.gb %>% 
  group_by(tag, algorithm) %>% count()

counts <- counts %>% 
  mutate(alg = if_else(algorithm == 1, '1-obj', if_else(algorithm == 2, "2-nov", "3-combo")))

# plot num of generations before global peak
counts %>% 
  ggplot(aes(x = alg, y = n)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(width = 0.25, shape=1, alpha = 0.5) + 
  theme_bw() +
  xlab("Search algorithms") + 
  ylab("Generations (til global peak or max run)")

# fitness changes
fit.max = max(df.gb$elite_fitness)

df.gb %>% 
  filter(trial %in% c("trial1", "trial10", "trial20")) %>% 
  filter(gen <25) %>% 
  ggplot(aes(x = gen, y = elite_fitness, group = tag, label = elite_hap)) + 
  geom_line() + 
  geom_point(shape = 21, size=2, fill="white") + 
  theme_bw() + 
  #labs(title = "G protein B1 domain landscape") +
  ylab("affinity") +
  geom_hline(yintercept = fit.max, linetype = 2) +
  #facet_wrap(~alg, nrow = 3) + 
  facet_grid(rows = vars(alg), cols = vars(trial)) +
  theme(legend.position = "bottom") +
  geom_text_repel(color = "gray") 

df.gb %>% 
  filter(trial %in% c("trial1", "trial10", "trial20", "trial30")) %>% 
  ggplot(aes(x = gen, y = elite_fitness, group = tag)) + 
  geom_line(alpha = 0.75) + 
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  labs(title = "G protein B1 domain landscape") +
  ylab("affinity") +
  geom_hline(yintercept = fit.max, linetype = 2) +
  facet_grid(rows = vars(alg), cols = vars(trial)) + 
  theme(legend.position = "bottom")

# objective search
pick_obj <- df.gb %>% filter(algorithm == 1 & trial == 'trial1' & gen < 20)
pick_obj %>% 
  ggplot(aes(x = gen, y = elite_fitness, label = elite_hap)) +
  geom_line() +
  geom_point(shape = 21, size = 4, fill = "white") + 
  geom_text_repel() + 
  geom_hline(yintercept = fit.max, linetype = 2) +
  theme_bw()

# novelty search
pick_nov <- df.gb %>% filter(algorithm == 2 & trial == 'trial20' & gen < 30)
pick_nov %>% 
  ggplot(aes(x = gen, y = elite_fitness, label = elite_hap)) +
  geom_line() +
  geom_point(shape = 21, size = 3, fill = "white") + 
  geom_text_repel() + 
  geom_hline(yintercept = fit.max, linetype = 2) +
  theme_bw()

# combo search
pick_combo <- df.gb %>% filter(algorithm == 3 & trial == 'trial10' & gen < 40)
pick_combo %>% 
  ggplot(aes(x = gen, y = elite_fitness, label = elite_hap)) +
  geom_line() +
  geom_point(shape = 21, size = 3, fill = "white") + 
  geom_text_repel() + 
  geom_hline(yintercept = fit.max, linetype = 2) +
  theme_bw()
#################################
# binary search
df <- read_tsv("data/nk10/search-obj-nov.combo.out") # binary string as behavior
#df <- read_tsv("data/nk10/nk10.out2") # normal fitness
#dfN <- read_tsv("data/nk15.out2")  # normal fitness
#dfE <- read_tsv("data/nk15-exp-fitness/nk15-exp.out2") # exponential fitness
#dfN <- dfN %>% mutate(fit_distri = 'normal')
#dfE <- dfE %>% mutate(fit_distri = 'exponential')

df <- df %>% mutate(alg = if_else(search_algo == 1, '1-obj', if_else(search_algo == 2, "2-nov", "3-combo")))

counts <- df %>% 
  group_by(Tag, search_algo, K) %>% count()

counts <- counts %>% 
  mutate(alg = if_else(search_algo == 1, '1-obj', if_else(search_algo == 2, "2-nov", "3-combo")))

# plot num of generations before global peak
counts %>% 
#  filter(n < 100) %>% # do not count failures
  ggplot(aes(x = alg, y = n)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(width = 0.2, aes(color = as.factor(K))) + 
  theme_bw() +
#  scale_y_log10() + 
#  facet_grid(rows = vars(fit_distri), cols = vars(alg)) + 
#  facet_wrap(~alg, nrow= 3) + 
#  theme() + 
  xlab("K") + 
  ylab("Generations (til global peak or max run)")

# fitness changes (Fig 3B)
df %>% 
  filter(K %in% c(0,2,5,9)) %>% 
  filter(str_detect(Tag, "-trial3$")) %>% 
  filter(Gen < 16) %>% 
  ggplot(aes(x = Gen, y = fit, color=as.factor(K))) +
  geom_line() + 
  geom_point(shape = 21, size=2, fill="white") +
  theme_bw() + 
  facet_wrap(~alg, nrow = 3) + 
  geom_hline(yintercept = 1, linetype = 2) 
#  theme(legend.position = "bottom") 
  #ylim(0,1.1)

#############
counts %>% 
  #  filter(n < 100) %>% # do not count failures
  ggplot(aes(x = as.factor(K), y = n)) +
  geom_boxplot(outlier.shape = NA ) +
  #geom_violin() +
  geom_jitter(width = 0.2, shape=1, color = "gray") + 
  theme_bw() +
  #scale_y_log10() + 
  #facet_grid(rows = vars(fit_distri), cols = vars(alg)) + 
  facet_wrap(~alg, nrow = 3) +
  #theme(legend.position = "bottom") + 
  xlab("K") + 
  ylab("Generations (til global peak or max run)")

counts %>% 
  ggplot(aes(x = as.factor(K), y = n, color = alg)) +
  geom_boxplot(outlier.shape = NA ) +
  #geom_violin() +
  geom_point(shape=1, position = position_jitterdodge()) + 
  theme_bw() +
  #scale_y_log10() + 
  #facet_grid(rows = vars(fit_distri), cols = vars(alg)) + 
  #facet_wrap(~alg, nrow = 3) +
  #theme(legend.position = "bottom") + 
  xlab("K") + 
  ylab("Generations (til global peak or max run)") + 
  theme(legend.position = "bottom")

# fig 3A
counts %>% 
  #  filter(n < 100) %>% # do not count failures
  ggplot(aes(x = K, y = n)) +
  #geom_point(outlier.shape = NA) +
  #geom_violin() +
  geom_smooth(method = "lm") + 
  geom_jitter(width = 0.2, shape=1, color = "gray") + 
  theme_bw() +
  #scale_y_log10() + 
  #facet_grid(rows = vars(fit_distri), cols = vars(alg)) + 
  facet_wrap(~alg, nrow = 3) +
  #theme(legend.position = "bottom") + 
  xlab("K") + 
  ylab("Generations (til global peak or max run)") + 
  scale_x_continuous(breaks = 0:9, labels = 0:9)


counts %>% 
  ggplot(aes(x = alg, y = n, colour = fit_distri)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_point(position=position_jitterdodge(), shape = 1, alpha = 0.4) +  
  theme_bw() +
  #  scale_y_log10() + 
  #facet_wrap(~alg) + 
  theme(legend.position = "bottom") + 
  xlab("search algorithms") + 
  ylab("Generations (til global peak or max run)")

#####################
# plot avg fitness

df2 <- df %>% group_by(alg, fit_distri, K, Gen) %>% summarise(fit_avg = mean(fit), fit_sd = sd(fit), fit_avg_norm = mean(fit)/max(fit))

# facet by algorithm
df2 %>% filter(fit_distri == 'normal' & K %in% c(0, 2, 5, 14)) %>% 
  ggplot(aes(x = Gen, y = fit_avg, group = K, color = as.factor(K))) + 
  geom_line(alpha = 0.75) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~alg) + 
  geom_hline(yintercept = 1, linetype = 2) + 
  theme(legend.position = "bottom")

df2 %>% filter(fit_distri != 'normal' & K %in% c(0, 2, 5, 14)) %>% 
  ggplot(aes(x = Gen, y = fit_avg_norm, group = K, color = as.factor(K))) + 
  geom_line(alpha = 0.75) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~alg) + 
  geom_hline(yintercept = 1, linetype = 2) + 
  theme(legend.position = "bottom")

# facet by algorithm
df2 %>% filter(fit_distri == 'normal' & K %in% c(0, 2, 5, 14)) %>% 
  ggplot(aes(x = Gen, y = fit_avg, group = K, color = as.factor(K))) + 
  geom_line(alpha = 0.75) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~alg) + 
  geom_hline(yintercept = 1, linetype = 2) + 
  theme(legend.position = "bottom")

# facet by K (not good)
df2 %>% filter(fit_distri == 'normal' & K %in% c(0, 2, 5, 14)) %>% 
  ggplot(aes(x = Gen, y = fit_avg_norm, group = K, color = alg)) + 
  geom_line(alpha = 0.75) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~as.factor(K)) + 
  geom_hline(yintercept = 1, linetype = 2) + 
  theme(legend.position = "bottom")

###################
# cov search
#######################
df.cov <- read_tsv("data/wuhan.out2") # wuhan 
df.cov <- df.cov %>% 
  mutate(alg = if_else(algorithm == 1, '1-obj', if_else(algorithm == 2, "2-nov", "3-combo")))
df.cov <- df.cov %>% 
  mutate(trial = str_remove(tag, "alg.-"))

counts <- df.cov %>% 
  group_by(tag, algorithm) %>% count()

counts <- counts %>% 
  mutate(alg = if_else(algorithm == 1, '1-obj', if_else(algorithm == 2, "2-nov", "3-combo")))

# plot num of generations before global peak
counts %>% 
  ggplot(aes(x = alg, y = n)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(width = 0.25, shape=1, alpha = 0.5) + 
  theme_bw() +
  xlab("Search algorithms") + 
  ylab("Generations (til global peak or max run)")

# fitness changes
bind.max = max(df.cov$elite_fitness)

df.cov %>% 
  filter(trial %in% c("trial1", "trial10", "trial20", "trial30")) %>% 
  ggplot(aes(x = gen, y = elite_fitness, group = tag, color = trial)) + 
  geom_line(alpha = 0.75) + 
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  labs(title = "SARS-Cov-2 RBD landscape by deep mutation", subtitle = "Genetic background: Wuhan-Hu-1; Length: 201 AA") +
  ylab("delta_bind") +
  geom_hline(yintercept = bind.max, linetype = 2) +
  facet_wrap(~alg, nrow = 3)

########################
# skip all below
# plot time
log <- read_tsv("data/nk15-exp-fitness/nk15-exp.logs")
colnames(log) <- c("run", "alg", "rep", "time")
log <- log %>% mutate(alg2 = if_else(alg == 1, '1-obj', if_else(alg == 2, "2-nov", "3-combo")), interact_site = as.integer(str_replace(run, "nk15-(\\d+).tsv", "\\1")))

log %>% filter(time < 100) %>% ggplot(aes(x = interact_site, y = time)) + 
  geom_jitter(shape = 1, width = 0.1) + 
  facet_wrap(~alg2, nrow = 3) +
#  scale_y_log10() + 
  theme_bw()

#################################
counts %>% 
  ggplot(aes(x = as.factor(search_algo), y = n)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(width = 0.2, shape=1, color = "gray") + 
  theme_bw() +
#  scale_y_log10() + 
#  facet_wrap(~as.factor(search_algo), nrow=3) + 
#  theme(legend.position = "bottom") + 
  xlab("Algorithm") + 
  ylab("Generations (til global peak or max run)")

counts <- counts %>% 
  mutate(interact = if_else(K < 3, "0-low", "1-high"))

counts %>%
#  filter(search_algo != 3) %>% 
  #  filter(n < 100) %>% # do not count failures
  ggplot(aes(color = as.factor(search_algo), y = n, x = interact)) +
  geom_boxplot(outlier.shape = NA ) +
#  geom_point(shape=1, alpha = 0.5) + 
  geom_point(position=position_jitterdodge(), shape = 1) +
  theme_bw() +
  scale_y_log10() + 
#  facet_wrap(~interact, nrow=3) + 
  theme(legend.position = "bottom") + 
  xlab("Num interactions") + 
  ylab("Generations (til global peak or max run)")



#plot success rate on y and k value on x 
#expect a downward trend

nk_df <- read_tsv("~/nov-search/nk_out2.tsv")

nk_counts <- nk_df %>% group_by(Tag) %>% count()
nk_counts <- nk_counts %>% mutate(alg = str_split(Tag,"-")[[1]][2]) %>% mutate(k = str_split(Tag,"-")[[1]][3])

nk_success_com <- nk_counts %>% filter(alg == "com") %>% group_by(k) %>% summarise(count = sum(n<100)) %>% mutate(alg = "com")

nk_success_nov <- nk_counts %>% filter(alg == "nov") %>% group_by(k) %>% summarise(count = sum(n<100)) %>% mutate(alg = "nov")

nk_success_obj <- nk_counts %>% filter(alg == "obj") %>% group_by(k) %>% summarise(count = sum(n<100)) %>% mutate(alg = "obj")

nk_success <- rbind(nk_success_com, nk_success_nov, nk_success_obj)

ggplot(nk_counts, aes(x = k, y = n)) + geom_jitter(shape=1) + facet_wrap(~alg, nrow=3) +scale_y_log10()

#plot number of success for each alg for each k
ggplot(nk_success, aes(x = k, y = count, color = alg)) + geom_boxplot() + geom_jitter()
ggplot(nk_success, aes(x = k, y = count, color = alg, group = alg)) + geom_line() + geom_jitter()