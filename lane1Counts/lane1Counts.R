library(tidyverse)
library(ggpubr)
library(Maaslin2)
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


largeIntestinalSamples<-metadata %>%
  filter(tissue =='colon') %>%
  .$sample
#there is probably an easier way to do this
d1D3FitnessScores<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% largeIntestinalSamples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()%>%
  prcomp(center = TRUE,
         scale = TRUE)
d1D3FitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = day,
             y = PC2))+
  geom_point()

smallIntestinalSamples<-metadata %>%
  filter(tissue =='dj') %>%
  .$sample
#there is probably an easier way to do this
d1D3FitnessScores<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% smallIntestinalSamples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()%>%
  prcomp(center = TRUE,
         scale = TRUE)
d1D3FitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = day,
             y = PC2))+
  geom_point()

maaslin2SiFitness<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% smallIntestinalSamples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()
maaslin2SiMetadata<-metadata%>%
  filter(tissue == 'dj')%>%
  column_to_rownames('sample')

fit_data = Maaslin2(
  input_data = maaslin2SiFitness,
  input_metadata = maaslin2SiMetadata,
  output = "day1Vday2Maaslin",
  normalization = 'none',
  transform = 'none',
  fixed_effects = c("day"))

maaslin2CoFitness<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% largeIntestinalSamples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()
maaslin2CoMetadata<-metadata%>%
  filter(tissue == 'colon')%>%
  column_to_rownames('sample')
fit_data = Maaslin2(
  input_data = maaslin2CoFitness,
  input_metadata = maaslin2CoMetadata,
  output = "day1Vday2MaaslinColon",
  normalization = 'none',
  transform = 'none',
  fixed_effects = c("day"))


day1Samples<-metadata %>%
  filter(day =='day1') %>%
  .$sample
#there is probably an easier way to do this
d1FitnessScores<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% day1Samples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()%>%
  prcomp(center = TRUE,
         scale = TRUE)
summary(d1FitnessScores)
d1FitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             label = sample,
             col = tissue,
             y = PC2))+
  geom_point()

maaslin2Day1Fitness<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% day1Samples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()
maaslin2CoMetadata<-metadata%>%
  filter(day == 'day1')%>%
  column_to_rownames('sample')
fit_data = Maaslin2(
  input_data = maaslin2Day1Fitness,
  input_metadata = maaslin2CoMetadata,
  output = "day1DjVColon",
  normalization = 'none',
  transform = 'none',
  fixed_effects = c("tissue"))

day3Samples<-metadata %>%
  filter(day =='day3') %>%
  .$sample
#there is probably an easier way to do this
d3FitnessScores<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% day3Samples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()%>%
  prcomp(center = TRUE,
         scale = TRUE)
summary(d3FitnessScores)
d3FitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             label = sample,
             col = tissue,
             y = PC2))+
  geom_point()

maaslin2Day1Fitness<-fitnessScoresLoci %>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% day1Samples)%>%
  column_to_rownames('sample')%>%
  as.data.frame()
maaslin2CoMetadata<-metadata%>%
  filter(day == 'day1')%>%
  column_to_rownames('sample')
fit_data = Maaslin2(
  input_data = maaslin2Day1Fitness,
  input_metadata = maaslin2CoMetadata,
  output = "day1DjVColon",
  normalization = 'none',
  transform = 'none',
  fixed_effects = c("tissue"))