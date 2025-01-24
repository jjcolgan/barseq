library(tidyverse)
library(vegan)
library(ggpubr)
fullGenesTabAndMeta=read_tsv('fullGenesTab.tsv')
tscores = read_tsv('full/fit_t.tab')

#these include experiments that did not 'pass'
fitnessScores= read_tsv('full/fit_logratios.tab')
metadata= read_tsv('fullbarseqMeta.txt')
metadata= metadata%>%
  mutate('mouseTissueDay'= paste(mouse, tissueDay, sep = '-'))

metadata$day = factor(metadata$day, levels = c('day1','day3', 'day7', 'day14'))
geneFitnessscores=fitnessScores %>%
  select(-desc,
         -sysName)%>%
  column_to_rownames('locusId')%>%
  t() %>%
  prcomp(center = T)
geneFitnessscores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = day,
             col = tissue))+
  geom_point()

mouseSamples = metadata %>%
  filter(tissue == 'dj'|
         tissue == 'colon')%>%
  .$sample

geneFitnessscores=fitnessScores %>%
  select(-desc,
         -sysName)%>%
  column_to_rownames('locusId')%>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% mouseSamples)%>%
  column_to_rownames('sample')%>%
  prcomp()

geneFitnessscores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = day,
             col = tissue))+
  geom_point()+
  labs(x = "PC1 - 25.82%",
       y = 'PC2 - 20.12%')
summary(geneFitnessscores)

adonisIn=fitnessScores %>%
  select(-desc,
         -sysName)%>%
  column_to_rownames('locusId')%>%
  t() %>%
  as.data.frame()
adonisMeta=metadata %>%
  filter(sample %in% rownames(adonisIn))
geneFitnessscores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             shape = day,
             col = tissue))+
  geom_point()

geneFitnessscores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = tissueDay))+
  geom_point()+
  stat_ellipse()

geneFitnessscores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissueDay))+
  geom_point()+
  stat_ellipse()

fitnessScoresLong=fitnessScores%>%
  select(-desc,
         -sysName)%>%
  pivot_longer(cols = c(2:43), values_to = 'log2ratio', names_to = 'sample')

tscoresLong=tscores %>%
  select(-desc,
         -sysName)%>%
  pivot_longer(cols = c(2:43), values_to = 'tScore', names_to = 'sample')

mergedFitnessScoresLong=fitnessScoresLong %>%
  left_join(tscoresLong, by = c('sample', 'locusId'))

signficantFitnessEffectsLong=mergedFitnessScoresLong %>%
  filter(abs(tScore) >= 5,
         abs(log2ratio)>=2)

significantFitnessEffectsWide=signficantFitnessEffectsLong %>%
  pivot_wider(names_from = 'locusId',
              id_cols = 'sample',
              values_from = 'log2ratio')
significantFitnessEffectsWide[is.na(significantFitnessEffectsWide)] = 0

signficantFitnessEffectsScores=significantFitnessEffectsWide %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(signficantFitnessEffectsScores)
signficantFitnessEffectsScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
  y = PC2,
  shape = day,
  col = tissue))+
  geom_point()+
  labs(x = 'PC1 = 33.41%',
       y = 'PC2 = 21.51%')

heatmapIn=signficantFitnessEffectsLong %>%
  left_join(metadata, by = 'sample')%>%
  pivot_wider(names_from = 'locusId',
              id_cols = 'mouseTissueDay',
              values_from = 'log2ratio')%>%
  column_to_rownames('mouseTissueDay')
heatmapIn[is.na(heatmapIn)] = 0
heatmapIn = t(heatmapIn)
ComplexHeatmap::Heatmap(heatmapIn, show_row_names = F)

significantFitnessPassingFilter=signficantFitnessEffectsLong%>%
  group_by(locusId)%>%
  summarise('numberOfObservations'= n()) %>%
  filter(numberOfObservations >= 20)%>%
  .$locusId

significantFitnessEffectsWideFiltered=signficantFitnessEffectsLong %>%
  filter(locusId %in% significantFitnessPassingFilter)%>%
  pivot_wider(names_from = 'locusId',
              id_cols = 'sample',
              values_from = 'log2ratio')

significantFitnessEffectsWideFiltered[is.na(significantFitnessEffectsWideFiltered)] = 0

significantFitnessEffectsWideFilteredScores=significantFitnessEffectsWideFiltered %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(significantFitnessEffectsWideFilteredScores)
significantFitnessEffectsWideFilteredScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = day,
             col = tissue))+
  geom_point()+
  labs(x = 'PC1 = 42.21%',
       y = 'PC2 = 24.09%')

heatmapInFiltered=signficantFitnessEffectsLong %>%
  filter(locusId %in% significantFitnessPassingFilter)%>%
  left_join(metadata, by = 'sample')%>%
  pivot_wider(names_from = 'locusId',
              id_cols = 'mouseTissueDay',
              values_from = 'log2ratio')%>%
  column_to_rownames('mouseTissueDay')
heatmapInFiltered[is.na(heatmapInFiltered)] = 0
heatmapInFiltered = t(heatmapInFiltered)
ComplexHeatmap::Heatmap(heatmapInFiltered, show_row_names = F)

highConfidenceFitnessEffectsLong=mergedFitnessScoresLong %>%
  filter(abs(tScore) >= 5)

highConfidenceFitnessEffectsWide=highConfidenceFitnessEffectsLong %>%
  pivot_wider(names_from = 'locusId',
              id_cols = 'sample',
              values_from = 'log2ratio')
highConfidenceFitnessEffectsWide[is.na(highConfidenceFitnessEffectsWide)] = 0

highConfidenceFitnessEffectsScores=highConfidenceFitnessEffectsWide %>%
  column_to_rownames('sample')%>%
  prcomp()
summary(highConfidenceFitnessEffectsScores)
highConfidenceFitnessEffectsScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = day,
             col = tissue))+
  geom_point()+
  labs(x = 'PC1 = 33.41%',
       y = 'PC2 = 21.51%')