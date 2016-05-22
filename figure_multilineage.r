library(ggplot2)
library(dplyr)
theme_set(theme_bw(10))

load('figure_multilineage.rdata')
qplot(data=output %>% filter(name!="virus"), x=t, y=frequency, geom="area" , fill=name)
ggsave("figures/multilineage.pdf", height=4, width=7)
