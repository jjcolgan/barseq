metadata = read_tsv('fullbarseqMeta.txt')

fitness = read_tsv('barseqAdjustedParams/fit_logratios.tab')
colnames(fitness) <- sub("setA", "", colnames(fitness))
colnames(fitness) <- sub("_.*", "", colnames(fitness))
colnames(fitness) <- sub("CO$", "Co", colnames(fitness))
colnames(fitness) <- sub("DJ$", "Dj", colnames(fitness))

fitnessLong = fitness %>%
  pivot_longer(cols = c(4:ncol(fitness)),
               names_to = 'sample',
               values_to = 'fitnessScore')

tscores = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/barseqAdjustedParams/fit_t.tab')
colnames(tscores) <- sub("setA", "", colnames(tscores))
colnames(tscores) <- sub("_.*", "", colnames(tscores))
colnames(tscores) <- sub("CO$", "Co", colnames(tscores))
colnames(tscores) <- sub("DJ$", "Dj", colnames(tscores))

tscoresLong = tscores %>%
  pivot_longer(cols = c(4:ncol(tscores)),
               names_to = 'sample',
               values_to = 'tscore')

mergedTandFitnessScoresLong=tscoresLong %>%
  left_join(fitnessLong,
            by = c('locusId',
                   'desc',
                   'sysName',
                   'sample'))%>%
  left_join(metadata, by = 'sample')

mergedTandFitnessScoresLong$weightedFitness = abs(mergedTandFitnessScoresLong$tscore)*mergedTandFitnessScoresLong$fitnessScore

mergedTandFitnessScoresLong%>%
  filter(tissue != 'T0')%>%
  group_by(tissue, locusId)%>%
  summarise('meanWeightedScore' = (mean(weightedFitness)))%>%
  pivot_wider(id_cols = 'locusId', names_from = 'tissue', values_from = 'meanWeightedScore')%>%
  ggplot(aes(x = colon, y =dj ))+
  geom_point(alpha = .05)

 pcaScaledOut=mergedTandFitnessScoresLong%>%
  filter(tissue != 'T0')%>%
  pivot_wider(names_from = 'locusId',
              id_cols = 'sample',
              values_from = 'weightedFitness')%>%
  column_to_rownames('sample')%>%
  prcomp(center = T,
         scale = T)

pcaScaledOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissue))+
  geom_point()

pcaOut=mergedTandFitnessScoresLong%>%
  filter(tissue != 'T0')%>%
  pivot_wider(names_from = 'locusId',
              id_cols = 'sample',
              values_from = 'weightedFitness')%>%
  column_to_rownames('sample')%>%
  prcomp(center = T)

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissue))+
  geom_point()


pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = day))+
  geom_point()