library(tidyverse)

df <- read_tsv("~/nov-search/out2.tsv")
df1 <- df %>% mutate(alg = str_remove(Tag, "-[0-9]+"))

counts <- df1 %>% group_by(Tag) %>% count()
counts <- counts %>% mutate(alg = str_remove(Tag, "-[0-9]+"))

ggplot(counts, aes(x = alg, y = n))+ geom_boxplot() +geom_jitter()

df <- read_tsv("~/Downloads/output_2.tsv")
df1 <- df %>% mutate(alg = str_remove(tag, "\\-[^-]*$"))

counts <- df1 %>% group_by(tag) %>% count()
counts <- counts %>% mutate(alg = str_remove(tag, "\\-[^-]*$"))

ggplot(counts, aes(x = alg, y = n))+ geom_boxplot() +geom_jitter()

counts %>% group_by(alg) %>% summarise(n<100) %>% table() %>% as_tibble() %>% ggplot(aes(x=alg, y=TRUE)) +geom_bar(stat="identity")

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