library(tidyverse)

#BD1
bd_df <- read_tsv("~/nov-search/bd1_out2.tsv")
bd_df1 <- bd_df %>% mutate(alg = str_remove(tag, "\\-[^-]*$"))

bd_counts <- bd_df1 %>% group_by(tag) %>% count()
bd_counts <- bd_counts %>% mutate(alg = str_remove(tag, "\\-[^-]*$"))

ggplot(bd_counts, aes(x = alg, y = n))+ geom_boxplot() +geom_jitter()
###aa search for BD1 shows that combo and obj search performed better than nov

#influenza
na_df <- read_tsv("~/nov-search/na_out2.tsv")
na_df1 <- na_df %>% mutate(alg = str_remove(tag, "\\-[^-]*$"))

na_counts <- na_df1 %>% group_by(tag) %>% count()
na_counts <- na_counts %>% mutate(alg = str_remove(tag, "\\-[^-]*$"))

ggplot(na_counts, aes(x = alg, y = n))+ geom_boxplot() +geom_jitter()
###aa search for influenza shows that combo performed the best out of the three