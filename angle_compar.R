library(tidyverse)

angle_compar <- read.csv("NK1_Middle_DMR_Angles.csv", sep = ";")

angle_compar <- angle_compar %>% 
  rowwise() %>% 
  mutate(mean = mean(c(X2, X3, X4), na.rm = T))


ggplot(angle_compar, aes(x=Program, y=mean)) + 
  geom_boxplot()


group_by(angle_compar, Program) %>%
  summarise(
    count = n(),
    median = median(mean, na.rm = TRUE),
    IQR = IQR(mean, na.rm = TRUE)
  )



ME78.1049g

library("ggpubr")


ggboxplot(angle_compar, x = "Program", y = "mean", 
          color = "Program", palette = c("#00AFBB", "#E7B800"),
          order = c("R", "Meshlab"),
          ylab = "mean", xlab = "Program")

ggline(angle_compar, x = "Program", y = "mean", 
       add = c("mean_se", "jitter"), 
       ylab = "Mean angle", xlab = "Program")


wilcox.test(mean ~ Program, data = angle_compar, paired = TRUE,)



R <- subset(angle_compar,  Program == "R", X3,
                 drop = TRUE)
# subset weight data after treatment
Meshlab <- subset(angle_compar,  Program == "Meshlab", X3,
                drop = TRUE)


wilcox.test(R, Meshlab, paired = TRUE)


library(PairedData)
pd <- paired(R, Meshlab)

plot(pd, type = "profile") + theme_bw()
