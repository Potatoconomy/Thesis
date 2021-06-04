library(ggplot2)
library(ggpubr)

path <- '/home/patrick/Desktop/Masters Immuno Figs/Counting/'
f_ctrl <- paste(path, '360_ctrl_Image.csv', sep='')
f_heat <- paste(path, '360_heat_Image.csv', sep='')

ctrl <- read.csv(f_ctrl)
heat <- read.csv(f_heat)

mg_ctrl <- unlist(ctrl["Count_IbaPositiveNuclei"])
mg_heat <- unlist(heat["Count_IbaPositiveNuclei"])



t.test(mg_ctrl, mg_heat)


mean_ctrl <- 46.167
mean_heat <- 62.370
std_ctrl <- sd(mg_ctrl)
std_heat <- sd(mg_heat)

df <- data.frame(group=c('Ctrl', 'Heat'), mean=c(mean_ctrl, mean_heat), std=c(std_ctrl, std_heat))

#barplot(height = c(mean_ctrl, mean_heat), names.arg = c('Ctrl', 'Heat'))

p <- ggplot(df) +
  geom_bar(aes(x=group, y=mean), stat="identity") +
  geom_errorbar( aes(x=group, ymin=mean-std, ymax=mean+std), width=0.4, colour="orange", alpha=1, size=1)


svg(filename = paste(path, 'CountPlot.svg'))
print(p)
dev.off()



ggdensity(ctrl$)

ggdensity(ctrl$Classify_IbaPos_NumObjectsPerBin, 
          main = 'Ctrl Cell Counts')

ggdensity(heat$Classify_IbaPos_NumObjectsPerBin,
          main = 'Heat Cell Counts')

