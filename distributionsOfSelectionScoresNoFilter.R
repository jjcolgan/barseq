library(tidyverse)

fitness = read_tsv('barseqAdjustedParams/fit_logratios.tab')
colnames(fitness) <- sub("setA", "", colnames(fitness))
colnames(fitness) <- sub("_.*", "", colnames(fitness))
colnames(fitness) <- sub("CO$", "Co", colnames(fitness))
colnames(fitness) <- sub("DJ$", "Dj", colnames(fitness))

metadata = read_tsv('fullbarseqMeta.txt')
metadata$cage = as.factor(metadata$cage)
metadata$lane = as.factor(metadata$lane)
quality=read_tsv('barseqAdjustedParams/fit_quality.tab')
quality$name <- sub("setA", "", quality$name)
quality$name <- sub("_.*", "", quality$name)
quality$name <- sub("CO$", "Co", quality$name)
quality$name <- sub("DJ$", "Dj", quality$name)
quality = quality %>%
  rename('sample'= name)

metadata = metadata %>%
  left_join(quality, by = 'sample')

metadata$day = factor(metadata$day, levels = c('t0', 'day1','day3', 'day7', 'day14'))


fitnessFullLong = fitness%>%
  pivot_longer(cols = c(4:45), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')
'Impact of most mutations are 0. Data is left skewed, which makes sense'
fitnessFullLong%>%
  ggplot(aes(x = fitnessScore))+
  geom_histogram()

fitnessFullLong%>%
  filter(day != 't0')%>%
  ggplot(aes(x = fitnessScore))+
  geom_histogram()+
  facet_wrap(~day, nrow = 1)

fitnessFullLong%>%
  filter(day != 't0')%>%
  ggplot(aes(x = fitnessScore))+
  geom_density()+
  facet_wrap(~day, ncol = 1)

fitnessFullLong%>%
  ggplot(aes(y = fitnessScore,
             x = day))+
  geom_boxplot(outliers = F)




fitnessFullLong%>%
  ggplot(aes(x = fitnessScore))+
  geom_density()

hist(fitnessFullLong$fitnessScore)

'Mean fitness score is weakly negative, -.38. Median is also weak, but more so
at - .05058. '
'From density and histogram, there are more observations of strongly negative
phenotypes than strongly positve'
summary(fitnessFullLong$fitnessScore)

fitnessFullLong%>%
  filter(tissue != 'T0')%>%
  ggplot(aes(x = fitnessScore))+
  geom_histogram()+
  facet_wrap(~tissue)

fitnessFullLong%>%
  filter(tissue != 'T0')%>%
  ggplot(aes(x = fitnessScore))+
  geom_density()+
  facet_wrap(~tissue, ncol =1)

fitnessFullLong%>%
  filter(tissue != 'T0') %>%
  ggplot(aes(x = fitnessScore))+
  geom_histogram()+
  facet_wrap(~tissue, ncol = 1)+
  geom_vline(xintercept = 1)+
  geom_vline(xintercept = -1)

fitnessFullLong%>%
  filter(tissue != 'T0') %>%
  ggplot(aes(x = fitnessScore))+
  geom_density()+
  facet_wrap(~tissue, ncol = 1)

ks.test(fitnessFullLong$fitnessScore[fitnessFullLong$tissue == "dj"],
          fitnessFullLong$fitnessScore[fitnessFullLong$tissue == "colon"])

ks.test(fitnessFullLong$fitnessScore[fitnessFullLong$tissue == "dj"],
        fitnessFullLong$fitnessScore[fitnessFullLong$tissue == "colon"],alternative = 'less')


fitnessFullLong%>%
  filter(tissue != 'T0') %>%
  ggplot(aes(x = tissue,y = fitnessScore))+
  geom_violin()

meanFitnessScoresTissue=fitnessFullLong %>%
  group_by(tissue, locusId) %>%
  summarise(meanFitnessScore = mean(fitnessScore)) %>%
  arrange(meanFitnessScore)

hundredGenesWithMeanLowestScoresDj=meanFitnessScoresTissue %>%
  filter(tissue == 'dj')%>%
  head(100)

hundredGenesWithMeanLowestScoresColon=meanFitnessScoresTissue %>%
  filter(tissue == 'colon')%>%
  head(100)
'Of the 100 genes with the lowest mean values across time points, only
69 are shared between tissues, 62 are not shared.'
sharedLowFitnessMutations=intersect(hundredGenesWithMeanLowestScoresColon$locusId,
hundredGenesWithMeanLowestScoresDj$locusId)

notSharedLowMutationGenes=rbind(hundredGenesWithMeanLowestScoresDj,hundredGenesWithMeanLowestScoresColon)%>%
  filter(!locusId %in% sharedLowFitnessMutations)%>%
  ungroup()%>%
  select(locusId)%>%
  distinct()%>%
  .$locusId


fitnessFullLong$color <- NA
fitnessFullLong$color[fitnessFullLong$locusId %in% sharedLowFitnessMutations] <- 'shared'
fitnessFullLong$color[fitnessFullLong$locusId %in% notSharedLowMutationGenes &
                        fitnessFullLong$locusId %in% hundredGenesWithMeanLowestScoresColon$locusId] <- 'not shared, low colon'
fitnessFullLong$color[fitnessFullLong$locusId %in% notSharedLowMutationGenes &
                        fitnessFullLong$locusId %in% hundredGenesWithMeanLowestScoresDj$locusId]<- 'not shared, low dj'

'Shared mutations seem like they have more deleterious effects in the colon, and earlier on in colonization'
'Day 14 DJ is the only one in which these genes cross to positive effects'
'Shared lowest mutations look to be have more deleterious consquences regardless of tissue
with the ones unique to a tissue being secont lowest'
fitnessFullLong %>%
  filter(tissue != 'T0') %>%
  ggplot(aes(x = fitnessScore, fill = factor(color))) +  # Ensure color is treated as a factor
  geom_histogram()+
  facet_wrap(~tissue,
           nrow = 2,
           ncol = 4)+
  geom_vline(xintercept = -2)+
  geom_vline(xintercept = -5)
fitnessFullLong %>%
  filter(tissue != 'T0') %>%
  ggplot(aes(x = fitnessScore)) +  # Ensure color is treated as a factor
  geom_histogram()+
  facet_wrap(~tissue+day,
             nrow = 2,
             ncol = 4)

fitnessFullLong %>%
  filter(tissue != 'T0') %>%
  ggplot(aes(x = fitnessScore)) +  # Ensure color is treated as a factor
  geom_density()+
  facet_wrap(~tissue+day,
             nrow = 2,
             ncol = 4)

ggsave(filename = 'histogramOfFitnessScores.pdf', width = 10, height = 10)
annotatedGenes = read_tsv('genesWithAnvioAnnotations.tsv')

annotatedGenes%>%
  filter(locusId %in% sharedLowFitnessMutations)

annotatedGenes%>%
  filter(locusId %in% notSharedLowMutationGenes & locusId%in% hundredGenesWithMeanLowestScoresDj$locusId)
annotatedGenes%>%
  filter(locusId %in% notSharedLowMutationGenes & locusId%in% hundredGenesWithMeanLowestScoresColon$locusId)

hundredGenesWithMeanHighestScoresDj=meanFitnessScoresTissue %>%
  filter(tissue == 'dj')%>%
  tail(100)

hundredGenesWithMeanHighestScoresColon=meanFitnessScoresTissue %>%
  filter(tissue == 'colon')%>%
  tail(100)


#56 shared high fitness genes
sharedHighFitnessMutations=intersect(hundredGenesWithMeanHighestScoresColon$locusId,
                                     hundredGenesWithMeanHighestScoresDj$locusId)

length(sharedHighFitnessMutations)

# 144 unique high fitness genes
notSharedHighMutationGenes=rbind(hundredGenesWithMeanHighestScoresDj,hundredGenesWithMeanHighestScoresColon)%>%
  filter(!locusId %in% sharedHighFitnessMutations)%>%
  ungroup()%>%
  select(locusId)%>%
  distinct()%>%
  .$locusId

notSharedHighMutationGenes%>% length()

fitnessFullLong$color <- NA
fitnessFullLong$color[fitnessFullLong$locusId %in% sharedHighFitnessMutations] <- 'shared high fitness'
fitnessFullLong$color[fitnessFullLong$locusId %in% notSharedHighMutationGenes &
                        fitnessFullLong$locusId %in% hundredGenesWithMeanHighestScoresColon$locusId] <- 'not shared, high colon'
fitnessFullLong$color[fitnessFullLong$locusId %in% notSharedHighMutationGenes &
                        fitnessFullLong$locusId %in% hundredGenesWithMeanHighestScoresDj$locusId]<- 'not shared, high dj'

fitnessFullLong %>%
  filter(tissue != 'T0') %>%
  ggplot(aes(x = fitnessScore, fill = factor(color))) +  # Ensure color is treated as a factor
  geom_histogram()+
  facet_wrap(~tissue)

fitnessFullLong %>%
  group_by(locusId)%>%
  summarise('geneFitnessVariance' = var(fitnessScore))%>%
  ggplot(aes(x = geneFitnessVariance))+
  geom_histogram()

lowestVarGenes=fitnessFullLong %>%
  group_by(locusId)%>%
  summarise('geneFitnessVariance' = var(fitnessScore))%>%
  arrange(geneFitnessVariance)%>%
  head(20)%>%
  .$locusId

highesetVarGenes=fitnessFullLong %>%
  group_by(locusId)%>%
  summarise('geneFitnessVariance' = var(fitnessScore))%>%
  arrange(geneFitnessVariance)%>%
  tail(20)%>%
  .$locusId

annotatedGenes %>%
  filter(locusId %in% lowestVarGenes)%>%
  view()

# a lot of these genes do not seem to have a huge impact on fitness
for (g in lowestVarGenes){
  p=fitnessFullLong %>%
    filter(locusId ==g)%>%
    ggplot(aes(x = fitnessScore))+
    geom_histogram()+
    facet_wrap(~tissue)
  plot(p)

}
# these plots are interesting and the genes themselves warrant
# more investigation, but generally the trends seem to be
#the same in both tissues over time
for (g in highesetVarGenes){
  p=fitnessFullLong %>%
    filter(locusId ==g)%>%
    ggplot(aes(x = fitnessScore,
               fill = day))+
    geom_histogram()+
    facet_wrap(~tissue)+
    labs(title = g)
  plot(p)
}

for (g in highesetVarGenes){
  p=fitnessFullLong %>%
    filter(locusId ==g)%>%
    ggplot(aes(y = fitnessScore,
               x = day,
               group = tissue,
               col = day))+
    geom_point()+
    geom_smooth(method = 'lm')+
    facet_wrap(~tissue)+
    labs(title = g)
  plot(p)
}

annotatedGenes %>%
  filter(locusId %in% highesetVarGenes)%>%
  view()


ftsXLoci=annotatedGenes %>%
  filter(grepl(name, pattern = 'ftsX') == T)%>%
  .$locusId
# really wierd to see fts proteins among those with highest variability given their role in cell division.
#FtsX is in that list, so want to look at just the dist of those proteins from the entire
#set

fitnessFullLong%>%
  filter(locusId %in% ftsXLoci)%>%
  ggplot(aes(x = fitnessScore,
             fill = day))+
  geom_histogram()+
  facet_wrap(~tissue)

fitnessFullLong%>%
  filter(tissue != 'T0')%>%
  filter(locusId %in% ftsXLoci)%>%
  ggplot(aes(y = fitnessScore,
             x = day,
             group = tissue,
             col = day))+
  geom_point()+
  facet_wrap(~tissue+locusId)+
  geom_smooth(method = 'lm')

ftsLoci=annotatedGenes %>%
  filter(grepl(name, pattern = 'fts') == T)%>%
  .$locusId

annotatedGenes %>%
  filter(locusId %in% ftsLoci)%>%
  view()

fitnessFullLong%>%
  filter(tissue != 'T0')%>%
  filter(locusId %in% ftsLoci)%>%
  ggplot(aes(y = fitnessScore,
             x = day,
             group = tissue,
             col = day))+
  geom_point()+
  facet_wrap(~tissue+locusId, nrow = 2)+
  geom_smooth()

perTissueFitnessVariances=fitnessFullLong %>%
  filter(tissue != 'T0')%>%
  group_by(tissue, locusId)%>%
  summarise('geneTissueVariance'= var(fitnessScore))

top10LowestVarDjGenes=perTissueFitnessVariances %>%
  filter(tissue == 'dj')%>%
  arrange(geneTissueVariance)%>%
  head(10)


top10LowestVarColonicGenes=perTissueFitnessVariances %>%
  filter(tissue == 'colon')%>%
  arrange(geneTissueVariance)%>%
  head(10)

'This is interesting, the gene variances in the colon are more or less
poisson distributed, while those in the dj have a heavy right skewed tail'
perTissueFitnessVariances %>%
  ggplot(aes(x = geneTissueVariance))+
  geom_histogram()+
  facet_wrap(~tissue)

'Only 4 of these are shared'
intersect(top10LowestVarColonicGenes$locusId,
          top10LowestVarDjGenes$locusId)

annotatedGenes %>%
  filter(locusId %in% top10LowestVarColonicGenes$locusId |
           locusId %in% top10LowestVarDjGenes$locusId)%>%
  view()
'Generally the trends are similar for each gene in each population
though the SI samples seem to have a greater spread than the colonic
samples'
for (g in top10LowestVarDjGenes$locusId){
  p=fitnessFullLong %>%
    filter(locusId == g)%>%
    ggplot(aes(x = day,
               group = tissue,
               y = fitnessScore))+
    geom_point()+
    facet_wrap(~tissue)+
    geom_smooth()+
    labs(title = g)
  plot(p)
}

for (g in top10LowestVarColonicGenes$locusId){
  p=fitnessFullLong %>%
    filter(locusId == g)%>%
    ggplot(aes(x = day,
               group = tissue,
               y = fitnessScore))+
    geom_point()+
    facet_wrap(~tissue)+
    geom_smooth()+
    labs(title = g)
  plot(p)
}

top10HighestVarDjGenes=perTissueFitnessVariances %>%
  filter(tissue == 'dj')%>%
  arrange(geneTissueVariance)%>%
  tail(10)

top10HighestVarColonicGenes=perTissueFitnessVariances %>%
  filter(tissue == 'colon')%>%
  arrange(geneTissueVariance)%>%
  tail(10)
'only 3 are shared'
intersect(top10HighestVarColonicGenes$locusId,
          top10HighestVarDjGenes$locusId)

annotatedGenes %>%
  filter(locusId %in% top10HighestVarDjGenes$locusId |
           locusId %in% top10HighestVarColonicGenes$locusId)%>%
  view()
'These are kinda of interesting'
for (g in top10HighestVarDjGenes$locusId){
  p=fitnessFullLong %>%
    filter(locusId == g)%>%
    ggplot(aes(x = day,
               group = tissue,
               y = fitnessScore))+
    geom_point()+
    facet_wrap(~tissue)+
    geom_smooth()+
    labs(title = g)
  plot(p)
}


