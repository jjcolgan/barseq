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

spearmanRes = read_tsv('colonicRNAClusterIntegration/spearmanRes/spearmanRes.tsv')
cluster1Expresssion = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/colonicKmeansResults/cluster1SigGenes.tsv')
cluster2Expresssion = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/colonicKmeansResults/cluster2SigGenes.tsv')
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

cluster1SpearmanRes = spearmanRes%>%
  filter(Gene %in% cluster1Expresssion$SYMBOL & padjust < .1)%>%
  mutate(cluster = 1)

cluster2SpearmanRes = spearmanRes%>%
  filter(Gene %in% cluster2Expresssion$SYMBOL & padjust < .1)%>%
  mutate(cluster = 2)

cluster3SpearmanRes = spearmanRes%>%
  filter(Gene %in% cluster3Expresssion$SYMBOL & padjust < .1)%>%
  mutate(cluster = 3)

allClustersSigSpearmanres=rbind(cluster1SpearmanRes,
      cluster2SpearmanRes,
      cluster3SpearmanRes)


allClustersSigSpearmanres$cluster = as.factor(allClustersSigSpearmanres$cluster)

allClustersSigSpearmanres$relationship = 'negative'
allClustersSigSpearmanres$relationship[allClustersSigSpearmanres$Spearman_Correlation > 0] = 'positive'

allClustersSigSpearmanres%>%
  group_by(cluster,
           relationship)%>%
  summarise(significantAssociations=n())%>%
  ggplot(aes(x = cluster,
             fill = cluster,
             y = significantAssociations))+
  geom_col()+
  facet_wrap(~relationship)+
  labs(x = 'Gene cluster',
       y ='Number of significant gene-microbial associations')

allClustersSigSpearmanres%>%
  group_by(cluster)%>%
  summarise(significantAssociations=n())%>%
  ggplot(aes(x = fct_reorder(cluster, significantAssociations),
             fill = cluster,
             y = significantAssociations))+
  geom_col()+
  labs(x = 'Gene cluster',
       y ='Number of significant gene-microbial associations')

allClustersSigSpearmanres%>%
  select(Mutant,cluster, relationship)%>%
  distinct()%>%
  group_by(cluster,
           relationship)%>%
  summarise(significantAssociations=n())%>%
  ggplot(aes(x = cluster,
             fill = cluster,
             y = significantAssociations))+
  geom_col()+
  facet_wrap(~relationship)+
  labs(x = 'Gene cluster',
       y ='Number of significant distinct microbe associations')

allClustersSigSpearmanres%>%
  select(Mutant,cluster, relationship)%>%
  distinct()%>%
  group_by(cluster)%>%
  summarise(significantAssociations=n())%>%
  ggplot(aes(x = fct_reorder(cluster, significantAssociations),
             fill = cluster,
             y = significantAssociations))+
  geom_col()+
  labs(x = 'Gene cluster',
       y ='Number of significant gene-microbial associations')

allClustersSigSpearmanresWide <- allClustersSigSpearmanres %>%
  distinct(Mutant, cluster) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = cluster, values_from = value, values_fill = 0)

upset(allClustersSigSpearmanresWide, intersect = colnames(allClustersSigSpearmanresWide)[2:4])+
  labs(title = 'Shared mutants by cluster')

allClustersSigSpearmanresWide <- allClustersSigSpearmanres %>%
  mutate(clusterRelationship = paste(cluster, relationship, sep = '-'))%>%
  distinct(Mutant, clusterRelationship) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = clusterRelationship, values_from = value, values_fill = 0)

upset(allClustersSigSpearmanresWide, intersect = colnames(allClustersSigSpearmanresWide)[2:7])+
  labs(title = 'Shared mutants by cluster and relationship')

mutantsIn=allClustersSigSpearmanresWide%>%
  filter(`1-negative` == 1, `1-positive` == 1)%>%
  .$Mutant

plotIn=allClustersSigSpearmanres %>%
  filter(cluster == 1,
         Mutant %in% mutantsIn)

# for (i in 1:50){
#   gene = plotIn$Gene[i]
#   mutant = plotIn$Mutant[i]
#   expressionTemp=expression%>%
#     select(gene)%>%
#     rownames_to_column('mergeCol')
#   fitnesstemp = fitnessScores%>%
#     select(mutant)%>%
#     rownames_to_column('mergeCol')
#
#   padj = plotIn%>%
#     filter(Gene == gene & Mutant == mutant)%>%
#     .$padjust
#
#   rho = plotIn%>%
#     filter(Gene == gene & Mutant == mutant)%>%
#     .$Spearman_Correlation
#
#   p=merge(expressionTemp,
#           by = 'mergeCol',
#           fitnesstemp) %>%
#     left_join(fitnessMeta, by = 'mergeCol')%>%
#     ggplot(aes(y = (.data[[gene]]),
#                x = .data[[mutant]],
#                group = tissue,
#                col = day))+
#     geom_point() +
#     geom_smooth(method = 'lm')+
#     labs(y = paste0('log2', gene),
#          x= mutant,
#          caption = paste0('Padj: ', padj, '\nRho: ', rho))
#
#   plot(p)
# }
#

annotatedDistinctMutants=allClustersSigSpearmanres%>%
  rename(locusId = Mutant)%>%
  left_join(keggs, by = 'locusId')%>%
  select(locusId,kofamAccession, cluster, relationship) %>%
  distinct()%>%
  annotate_modules()

annotatedDistinctMutants%>%
  ggplot(aes(y =moduleLevel1,
             fill = moduleLevel2 ))+
  geom_bar()+
  facet_wrap(~relationship+cluster)+
  labs(title = 'Cluster 2 negative assocication mutant function')