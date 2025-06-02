
library(clusterProfiler)
library(KEGGREST)
library(ggplot2)
library(purrr)
library(progressr)
library(tidyverse)
library(ComplexUpset)

get_module <- function(accession) {
  keggout <- keggGet(accession)
  if (!is.null(keggout[[1]]$BRITE)) {
    brite <- keggout[[1]]$BRITE
    level1 <- if (length(brite) >= 2) as.character(brite[2]) else NA_character_
    level2 <- if (length(brite) >= 3) as.character(brite[3]) else NA_character_
    return(c(level1, level2))
  } else {
    return(c(NA_character_, NA_character_))
  }
}

annotate_modules <- function(df) {
  df_no_na <- df %>% filter(!is.na(kofamAccession)) %>%
    mutate(moduleLevel1 = NA_character_, moduleLevel2 = NA_character_)

  with_progress({
    p <- progressor(steps = nrow(df_no_na))
    for (i in seq_len(nrow(df_no_na))) {
      p()
      keggOut <- get_module(df_no_na$kofamAccession[i])
      df_no_na$moduleLevel1[i] <- keggOut[1]
      df_no_na$moduleLevel2[i] <- keggOut[2]
    }
  })

  df_annotated <- df %>%
    left_join(df_no_na %>% select(kofamAccession, moduleLevel1, moduleLevel2),
              by = "kofamAccession")

  return(df_annotated)
}

spearmanRes = read_tsv('djClusterIntegration/spearmanRes/sigDjBarseqGenesStrongMutantSpearmanCorrelationEstimatedPvalues.tsv')
clusters = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqJejunumOutputs/manualClusteringClusters.tsv')
expression = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqJejunumOutputs/lrtResults/significantGenes.tsv')
fitnessMeta = read_tsv('fullbarseqMeta.txt')
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction)%>%
  distinct()
metadataRNA = read.csv('metadataRnaSeq.csv')%>%
  as.data.frame()

sampleLibraryRna = metadataRNA %>%
  filter(Tissue == 'Jejunum'& Treatment == 'Bar-seq')%>%
  select(sample, library)
'Clean up column names'
expression=expression %>%
  rename_with(~ sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", .))


fitnessScores = read_tsv('barseqAdjustedParams/fit_logratios.tab')
colnames(fitnessScores) <- sub("setA", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("_.*", "", colnames(fitnessScores))
colnames(fitnessScores) <- sub("CO$", "Co", colnames(fitnessScores))
colnames(fitnessScores) <- sub("DJ$", "Dj", colnames(fitnessScores))
fitnessScores=fitnessScores %>%
  select(-c(sysName, desc))%>%
  as.data.frame()%>%
  column_to_rownames('locusId')%>%
  t()

fitnessMeta=fitnessMeta%>%
  filter(tissue == 'dj')
fitnessMeta=fitnessMeta %>%
  mutate(mergeCol = paste0(mouse, 'Je'))

fitnessMeta = fitnessMeta %>%
  distinct(mergeCol, .keep_all = TRUE)


expression= expression[rownames(expression) %in% fitnessMeta$mergeCol, , drop = FALSE]

fitnessMeta=fitnessMeta %>%
  filter(mergeCol %in% row.names(expression))
fitnessScores=fitnessScores %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(fitnessMeta[c('sample', 'mergeCol')], by = 'sample')

expression= expression[rownames(expression) %in% fitnessScores$mergeCol, , drop = FALSE]
fitnessScores = fitnessScores %>%
  select(-sample)%>%
  column_to_rownames('mergeCol')

spearmanRes$padjust%>%
  summary()

sigCor = spearmanRes%>%
  filter(padjust < .05)

sigCor%>%
  ggplot(aes(x = abs(Spearman_Correlation),
                     fill = padjust))+
  geom_histogram()

mergedClusterSpearmanRes=clusters %>%
  rename_with(~"Gene", .cols = 1)%>%
  left_join(sigCor, by = 'Gene', relationship = 'many-to-many')

mergedClusterSpearmanRes%>%
  select(Gene, Mutant, clusterId)%>%
  mutate(clusterId = as.factor(clusterId))%>%
  distinct()%>%
  group_by(clusterId)%>%
  summarise(numberOfAssociations = n())%>%
  ggplot(aes(x = fct_reorder(clusterId, numberOfAssociations),
             y = log10(numberOfAssociations)))+
  geom_col()

'There are like 40K genexbug assoications and like idk 90% of them are in cluster 2'
mergedClusterSpearmanRes%>%
  select(Gene, Mutant, clusterId)%>%
  mutate(clusterId = as.factor(clusterId))%>%
  distinct()%>%
  group_by(clusterId)%>%
  summarise(numberOfAssociations = n())%>%
  ggplot(aes(x = fct_reorder(clusterId, numberOfAssociations),
             y = numberOfAssociations))+
  geom_col()+
  labs(x = 'Cluster ID',
       y = 'Number of gene-mutant assoications',
       title = 'Number of gene-mutant assoications with each cluster')


mergedClusterSpearmanRes%>%
  select(Mutant, clusterId)%>%
  mutate(clusterId = as.factor(clusterId))%>%
  distinct()%>%
  group_by(clusterId)%>%
  summarise(numberOfAssociations = n())%>%
  ggplot(aes(x = fct_reorder(clusterId, numberOfAssociations),
             y = numberOfAssociations))+
  geom_col()+
  labs(x = 'Cluster ID',
       y = 'Number of distinct mutants assoicated with cluster',
       title = 'Number of mutants assoicated with each cluster')

mergedClusterSpearmanRes$relationship = 'Negative'
mergedClusterSpearmanRes$relationship[mergedClusterSpearmanRes$Spearman_Correlation > 0 ] = 'Positive'

annotatedDistinctMutants=mergedClusterSpearmanRes%>%
  rename(locusId = Mutant)%>%
  left_join(keggs, by = 'locusId')%>%
  select(locusId,kofamAccession) %>%
  distinct()%>%
  annotate_modules()

annotatedDistinctMutants%>%
  rename(Mutant = locusId)%>%
  left_join(mergedClusterSpearmanRes, by = 'Mutant')%>%
  select(Mutant,
         moduleLevel1,
         moduleLevel2,
         clusterId,
         relationship)%>%
  distinct()%>%
  ggplot(aes(y =moduleLevel1,
             fill = moduleLevel2 ))+
  geom_bar()+
  facet_wrap(~clusterId+relationship)+
  labs(title = 'Mutant functions and relationship with each cluster jejunum padj < .05',
       y = 'BRITE level 1',
       fill = 'BRITE level 2')

