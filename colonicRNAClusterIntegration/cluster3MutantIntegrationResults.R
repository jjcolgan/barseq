library(clusterProfiler)
library(KEGGREST)
library(ggplot2)
library(purrr)
library(progressr)
library(tidyverse)

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

spearmanRes = read_tsv('colonicRNAClusterIntegration/spearmanRes/spearmanRes.tsv')
cluster3Expresssion = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/colonicKmeansResults/cluster3SigGenes.tsv')
expression = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/sigGeneTab.tsv')
fitnessMeta = read_tsv('fullbarseqMeta.txt')
sigMutants = read_tsv('colonLogRatiosMaaslin2/significant_results.tsv')
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction)%>%
  distinct()



metadataRNA = read.csv('metadataRnaSeq.csv')%>%
  as.data.frame()

sampleLibraryRna = metadataRNA %>%
  filter(Tissue == 'Colon'& Treatment == 'Bar-seq')%>%
  select(sample, library)
'Clean up column names'
expression=expression %>%
  rename_with(~ sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", .))

expression= expression %>%
  rename('2490Co'='2490Co-d')%>%
  column_to_rownames('SYMBOL')%>%
  t()%>%
  as.data.frame()

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
  filter(tissue == 'colon')
fitnessMeta=fitnessMeta %>%
  mutate(mergeCol = paste0(mouse, 'Co'))

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

spearmanRes=spearmanRes%>%
  filter(Gene %in% cluster3Expresssion$SYMBOL)

sigRes=spearmanRes%>%
  filter(padjust < .1)

'In this case a positive assoication should also be decreasing in fitness as expression decreases. This
cluster trends down over time'

'10,031 mutant-gene pairs'
positiveRelationship = sigRes %>%
  filter(Spearman_Correlation > 0)
'22,471 negative mutant-gene pairs'
negativeRelationship = sigRes %>%
  filter(Spearman_Correlation < 0)

sigRes%>%
  pivot_wider(names_from = 'Gene',
              id_cols = "Mutant",
              values_from = 'Spearman_Correlation',
              values_fill = 0)%>%
  column_to_rownames('Mutant')%>%
  ComplexHeatmap::Heatmap()

plotIn=positiveRelationship%>%
  arrange(desc(-log10(padjust)), Spearman_Correlation)%>%
  head(20)

for (i in 1:nrow(plotIn)){
  gene = plotIn$Gene[i]
  mutant = plotIn$Mutant[i]
  expressionTemp=expression%>%
    select(gene)%>%
    rownames_to_column('mergeCol')
  fitnesstemp = fitnessScores%>%
    select(mutant)%>%
    rownames_to_column('mergeCol')

  padj = positiveRelationship%>%
    filter(Gene == gene & Mutant == mutant)%>%
    .$padjust

  rho = positiveRelationship%>%
    filter(Gene == gene & Mutant == mutant)%>%
    .$Spearman_Correlation

  p=merge(expressionTemp,
          by = 'mergeCol',
          fitnesstemp) %>%
    left_join(fitnessMeta, by = 'mergeCol')%>%
    ggplot(aes(y = (.data[[gene]]),
               x = .data[[mutant]],
               group = tissue,
               col = day))+
    geom_point() +
    geom_smooth(method = 'lm')+
    labs(y = paste0('log2', gene),
         x= mutant,
         caption = paste0('Padj: ', padj, '\nRho: ', rho))

  plot(p)
}

plotIn=negativeRelationship%>%
  arrange(desc(-log10(padjust)), Spearman_Correlation)%>%
  head(20)

for (i in 1:nrow(plotIn)){
  gene = plotIn$Gene[i]
  mutant = plotIn$Mutant[i]
  expressionTemp=expression%>%
    select(gene)%>%
    rownames_to_column('mergeCol')
  fitnesstemp = fitnessScores%>%
    select(mutant)%>%
    rownames_to_column('mergeCol')

  padj = negativeRelationship%>%
    filter(Gene == gene & Mutant == mutant)%>%
    .$padjust

  rho = negativeRelationship%>%
    filter(Gene == gene & Mutant == mutant)%>%
    .$Spearman_Correlation

  p=merge(expressionTemp,
          by = 'mergeCol',
          fitnesstemp) %>%
    left_join(fitnessMeta, by = 'mergeCol')%>%
    ggplot(aes(y = (.data[[gene]]),
               x = .data[[mutant]],
               group = tissue,
               col = day))+
    geom_point() +
    geom_smooth(method = 'lm')+
    labs(y = paste0('log2', gene),
         x= mutant,
         caption = paste0('Padj: ', padj, '\nRho: ', rho))

  plot(p)
}


annotatedNegativeRelationship=negativeRelationship%>%
  rename(locusId = Mutant)%>%
  left_join(keggs, by = 'locusId')%>%
  select(locusId,kofamAccession) %>%
  distinct()%>%
  annotate_modules()

p=annotatedNegativeRelationship%>%
  select(locusId, moduleLevel2, moduleLevel1) %>%
  distinct()%>%
  ggplot(aes(y =moduleLevel1,
             fill = moduleLevel2 ))+
  geom_bar()+
  xlim(0,100)+
  labs(title = 'Cluster 3 negative mutant function')

ggsave(p,
       file = 'colonicRNAClusterIntegration/spearmanRes/cluster3NegativeAssociationMutantFunctions.pdf',
       height = 20,
       width = 20)

annotatedPositiveRelationship=positiveRelationship%>%
  rename(locusId = Mutant)%>%
  left_join(keggs, by = 'locusId')%>%
  select(locusId,kofamAccession) %>%
  distinct()%>%
  annotate_modules()

p=annotatedPositiveRelationship%>%
  ggplot(aes(y =moduleLevel1,
             fill = moduleLevel2 ))+
  geom_bar()+
  xlim(0,100)+
  labs(title = 'Cluster 3 positive mutant function')

ggsave(p,
       file = 'colonicRNAClusterIntegration/spearmanRes/cluster3PositiveAssociationMutantFunctions.pdf',
       height = 20,
       width = 20)
