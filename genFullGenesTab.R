library(tidyverse)
library(ggpubr)
library(Maaslin2)
library(ComplexHeatmap)
library(vegan)

fullLocusFitness<-read_tsv('full/fit_t.tab')

metadata = read_tsv('fullbarseqMeta.txt')
metadata$numericDay = 0
metadata$numericDay[metadata$day== 'day1'] = 1
metadata$numericDay[metadata$day== 'day3'] = 3
metadata$numericDay[metadata$day== 'day 7'] = 7
metadata$numericDay[metadata$day== 'day 14'] = 14
metadata$numericDay = as.integer(metadata$numericDay)
metadata=metadata%>%
  mutate(mouseDay = paste(mouse, day, sep = '-'))
metadata=metadata%>%
  mutate(mouseDayTissue = paste(mouse, day, tissue,sep = '-'))

genes=read_tsv('genesWithAnvioAnnotations.tsv')
fullLocusFitnessTransposed<-fullLocusFitness %>%
  pivot_longer(cols = 4:45, values_to = 'fitnessScore', names_to = 'sample') %>%
  left_join(metadata, by ='sample')

fullLocusFitnessTransposed=fullLocusFitnessTransposed %>%
  left_join(genes, by = c('locusId', "sysName", "desc" ))%>%
  as.data.frame()

write_tsv(fullLocusFitnessTransposed, file = 'fullGenesTab.tsv')