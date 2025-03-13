'For analysis of mutants with change in mean rank betweeen
colonic and dj tissues of less than 100, and ranked within
the top 100 lowest ranked fitness mutants in both tissues '

library(tidyverse)

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
                   'sample'))

lowRankMutants=read_tsv('conservedLowRankAcrossTissues.tsv')
annotatedMutants = read_tsv('genesWithAnvioAnnotations.tsv')

keggOrthology = read_tsv('/Users/johnjamescolgan/Downloads/KO_Orthology_ko00001 (3).txt', col_names = F)
keggOrthology$KO<-str_split_i(keggOrthology$X4, pattern = ' ', i = 1)
keggAnnotations=annotatedMutants%>%
  select(locusId, desc, kofamFunction, kofamAccession)%>%
  distinct()

keggAnnotations$
pfamAnnotations= annotatedMutants%>%
select(locusId, desc, pfamFunction, pfamAcession)%>%
distinct()

'Can use the t-scores and binomial statistics to determine if this is
significant(?)'
lowRankMutants%>%
  left_join(keggAnnotations, by = 'locusId')%>%
  view()

lowRankMutants%>%
  mutate(geneFunction = '')%>%
  left_join(keggAnnotations, by = 'locusId')%>%
  rename(KO = kofamAccession)%>%
  left_join(keggOrthology,by = 'KO' )%>%
  ggplot(aes(x = X1, fill = X2))+
  geom_bar()+
  labs(fill = 'Kegg pathway',
       x = 'Kegg superpathway',
       title = 'Low rank genes across tissues')+
  facet_wrap(~X1, scales = 'free_x')+
  theme_bw()


'Carbohydrate metabolism seems to be the largest component,
need to figure out how to do some sort of overrepresentation analysis'
lowRankMutants%>%
  left_join(pfamAnnotations, by = 'locusId')%>%
  view()

significantObservationsPerLocus=mergedTandFitnessScoresLong %>%
  filter(locusId %in% lowRankMutants$locusId)%>%
  filter(abs(tscore)>= 2)%>%
  group_by(locusId)%>%
  summarise('significantObservationsPerGene' = n())

significantObservations =mergedTandFitnessScoresLong%>%
  filter(locusId %in% lowRankMutants$locusId)%>%
  filter(abs(tscore)>= 2)%>%
  left_join(metadata, by = 'sample')


'Can take this and use binomial statistics to determine whether this was due to chance per group, per tissue,
and overall '
significantObservationsPerGroup = significantObservations%>%
  filter(tissue != 'T0')%>%
  group_by(tissueDay, locusId)%>%
  summarise('perLocusSignificantGroupObservations' = n())%>%
  ungroup()
'I dont think this makes sense. Trying to establish whether the mutatation is over all deletrious for the bacteria,
need to test sample wide, not just at each time point right?'
nMutants = fitness%>%
  select(locusId)%>%
  distinct()%>%
  nrow()

loci = lowRankMutants%>%
  select(locusId)%>%
  distinct()%>%
  .$locusId

prob = 1 / nMutants

groups = significantObservationsPerGroup%>%
  select(tissueDay)%>%
  distinct()%>%
  .$tissueDay

groupCounts=metadata %>%
  filter(sample %in% colnames(fitness))%>%
  group_by(tissueDay)%>%
  summarise(groupCounts = n())
output = data.frame('locus'=character(),
                    'group' = character(),
                    'pvalue' = double())

for (g in groups) {
  for (l in loci) {
    testDf = significantObservationsPerGroup %>%
      filter(tissueDay == g, locusId == l)

    if (nrow(testDf) > 0) {
      nSuccess = testDf$perLocusSignificantGroupObservations

      nTrials = groupCounts %>%
        filter(tissueDay == g) %>%
        pull(groupCounts)

      binomOut = binom.test(x = nSuccess,
                            n = nTrials,
                            p = prob)

      output = rbind(output,
                     data.frame('locus' = l,
                                'group' = g,
                                'pvalue' = binomOut$p.value))
    }
  }
}

output$padjust = p.adjust(output$pvalue, method = 'BH')

output%>%
  filter(padjust < .05)
output$group = factor(output$group, levels = c('colonday1',
                                               'colonday3',
                                               'colonday7',
                                               'colonday14',
                                               'djday1',
                                               'djday3',
                                               'djday7',
                                               'djday14'))

nTrials = metadata %>%
  filter(sample %in% colnames(fitness))%>%
  nrow()
output %>%
  rename(locusId = locus)%>%
  left_join(keggAnnotations,by = 'locusId')%>%
  left_join(keggOrthology, by = )
  ggplot(aes(y = kofamAccession,
             x = group,
             fill = log10(padjust)))+
  geom_tile()

'testing overall significance'
output = data.frame('locus'=character(),
                   'pvalue' = double())

for (l in loci) {
  testDf = significantObservations %>%
    filter(locusId == l)

  if (nrow(testDf) > 0) {
    nSuccess = nrow(testDf)
    binomOut = binom.test(x = nSuccess,
                          n = nTrials,
                          p = prob)

    output = rbind(output,
                   data.frame('locus' = l,
                              'pvalue' = binomOut$p.value))
  }
}

output$padjust = p.adjust(output$pvalue, method = 'BH')

output %>%
  filter(padjust < 0.05) %>%
  rename(locusId = locus) %>%
  left_join(keggAnnotations, by = "locusId") %>%
  rename(KO = kofamAccession) %>%
  left_join(keggOrthology, by = "KO") %>%
  mutate(X1 = if_else(is.na(X1), "Unknown or unintegrated", X1)) %>%
  mutate(X1 = fct_infreq(X1)) %>%
  mutate(lowFitnessPhenotypes = '') %>%
  ggplot(aes(y= lowFitnessPhenotypes, fill = X1)) +
  geom_bar(position = "fill") +
  labs(fill = 'KEGG superpathway', y = "Proportion") +
  theme_bw()