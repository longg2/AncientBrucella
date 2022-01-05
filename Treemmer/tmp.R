library(ggplot2)
library(dplyr)

dat <- read.delim("ST11RTL.tab", header = F, col.names = c("RTL", "Leaves"))

dat %>% ggplot(aes(x = Leaves, y = RTL)) + geom_point() + scale_x_reverse() + theme_bw() + geom_hline(yintercept = 0.95, col = "red", lty = 2)
ggsave("TreemerTrimmingST11.pdf", width = 6, height = 4)
