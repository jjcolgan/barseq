library(tidyverse)
library(ggpubr)
library()
metadata = read_tsv('fullbarseqMeta.txt')
counts=read_tsv('lane1Counts/combine.poolcount')
counts=counts[c(1,6:25)]
counts =counts %>%
  select(! 'EC-CG-2s-pl1-Filler68')
countsTransposed <- counts %>%
  column_to_rownames('barcode')%>%
  t()
countsTransposed = as.data.frame(countsTransposed)
countsTransposedFiltered <- countsTransposed[, sapply(countsTransposed, function(x) var(x, na.rm = TRUE) > 0)]
barseqScores<-countsTransposedFiltered %>%
  prcomp(center =TRUE, scale = TRUE)
barseqScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  left_join(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1,
             col = tissue,
             shape = day,
             y = PC2))+
  geom_point()
'There are 3 samples missing from the colonic groups, and 1 from the small intestine. '
barseqScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  left_join(metadata, by = 'sample')%>%
  view()

strainFitness = read_tsv('lane1Counts/barseqResults/strain_fit.tab')
strainFitnessFiltered=strainFitness %>%
  filter(used ==T)
strainFitnessFiltered=strainFitnessFiltered[c(1,10:29)]
strainFitnessFiltered = strainFitnessFiltered %>%
  select(!'EC_CG_2s_pl1_Filler68')
strainFitnessFilteredTransposed <- strainFitnessFiltered %>%
  column_to_rownames('barcode')%>%
  t()

strainFitnessFilteredTransposed = as.data.frame(strainFitnessFilteredTransposed)
strainFitnessFilteredTransposed <- strainFitnessFilteredTransposed[, sapply(strainFitnessFilteredTransposed, function(x) var(x, na.rm = TRUE) > 0)]


strainFitnessScores<-strainFitnessFilteredTransposed %>%
  prcomp(center =TRUE, scale = TRUE)
summary(strainFitnessScores)
pcaData<-strainFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')
pcaData$sample <- gsub("_", "-", pcaData$sample)

pcaData %>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = tissue,
             shape = day,
             label = sample,
             y = PC2))+
  geom_point()+
  geom_label()

pcaData %>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             col = tissue,
             shape = day,
             label = sample,
             y = PC3))+
  geom_point()

pcaData %>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = tissue,
             label = sample,
             y = PC2))+
  geom_point()

fitnessScores = read_tsv('lane1Counts/barseqResults/fit_t.tab')
fitnessScoresLoci= fitnessScores[c(1, 4:23)]
fitnessScoresLoci=fitnessScoresLoci%>%
  select(!'EC-CG-2s-pl1-Filler68')
fitnessScoresLociTransposed<-fitnessScoresLoci %>%
  column_to_rownames('locusId') %>%
  t()

locusFintessScores<-fitnessScoresLociTransposed %>%
  prcomp(center =TRUE, scale = TRUE)
summary(locusFintessScores)
locusFintessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = tissue,
             label = sample,
             shape = day,
             y = PC2))+
  geom_point()+
  labs(title = 'locus fitness PCA',
       y = 'PC2 -17.11% of variance',
       x = 'PC1 - 22.39% of variance')


locusFintessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             col = tissue,
             label = sample,
             shape = day,
             y = PC3))+
  geom_point()


smallIntestinalSamples<-metadata %>%
  filter(tissue =='dj') %>%
  .$sample

fitnessScoresLoci%>%
  select(all_of(smallIntestinalSamples))