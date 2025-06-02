
library(parallel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)
library(progressr)
library(tidyverse)
library(KEGGREST)
library(scales)


compute_correlation <- function(pair, method) {
  gene <- pair[[1]]
  mutant <- pair[[2]]
  test <- cor.test(scale(expression[, gene]), scale(fitnessScores[, mutant]), method = method, exact = F)
  return(data.frame(Gene = gene,
                    Mutant = mutant,
                    Spearman_Correlation = test$estimate,
                    P_Value = test$p.value))
}

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

expression = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqJejunumOutputs/lrtResults/significantGenes.tsv')
fitnessMeta = read_tsv('fullbarseqMeta.txt')
sigMutants = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/linear models/mixedEffectsJejunumOnly/mixedEffectsRes.tsv')
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction)%>%
  distinct()

mutantsToTest=sigMutants %>%
  as.data.frame()%>%
  select(-sample)%>%
  colnames()

metadataRNA = read.csv('metadataRnaSeq.csv')%>%
  as.data.frame()

sampleLibraryRna = metadataRNA %>%
  filter(Tissue == 'Jejunum'& Treatment == 'Bar-seq')%>%
  select(sample, library)
'Clean up column names'
expression=expression %>%
  rename_with(~ sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", .))

expression= expression %>%
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
  filter(locusId %in% mutantsToTest)%>%
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

fitnessScores <- fitnessScores[order(rownames(fitnessScores)), , drop = FALSE]
expression <- expression[order(rownames(expression)), , drop = FALSE]



# Generate all gene-mutant pairs
gene_mutant_pairs <- expand.grid(colnames(expression), colnames(fitnessScores), stringsAsFactors = FALSE)

# Use parallel processing
num_cores <- 4
cor_results <- do.call(rbind, mclapply(1:nrow(gene_mutant_pairs), function(i) {
  compute_correlation(gene_mutant_pairs[i, ], method = 'spearman')
}, mc.cores = num_cores))

# Adjust p-values
cor_results$padjust <- p.adjust(cor_results$P_Value, method = 'fdr')

# Filter significant results
significant_results <- cor_results %>%
  filter(padjust < 0.1)


significant_results %>%
  filter(Spearman_Correlation > 0)%>%
  nrow()

significant_results%>%
  filter(Spearman_Correlation < 0)%>%
  nrow()

cor_results%>%
  ggplot(aes(x = P_Value))+
  geom_histogram(bins = 100)


cor_results%>%
  arrange(desc(-log10(P_Value)))%>%
  view()
cor_results%>%
  write_tsv(file = 'djClusterIntegration/spearmanRes/sigDjBarseqGenesStrongMutantSpearmanCorrelationEstimatedPvalues.tsv')

mutants = unique(significant_results$Mutant)

# for (m in mutants){
#   assoications=significant_results%>%
#     filter(Mutant==m)%>%
#     .$Gene
#   print(m)
#   goOut=enrichGO(gene = assoications,
#            OrgDb=org.Mm.eg.db,
#            ont = 'bp',
#            keyType = 'ENSEMBL',
#            pvalueCutoff = .1 ,
#            pAdjustMethod = 'fdr'
#   )
#   goOut=goOut%>%
#     as.data.frame()%>%
#     filter(Count >= 3)
# print(goOut)
# }

cor_results %>%
  filter(Gene %in% significant_results$Gene)%>%
  filter(Mutant %in% significant_results$Mutant)%>%
  ggplot(aes(x = Mutant,
             y = Gene,
             fill = Spearman_Correlation))+
  geom_tile()

plotIn=significant_results%>%
  arrange(desc(-log10(padjust)), Spearman_Correlation)%>%
  head(20)

for (i in 1:nrow(plotIn)){
  gene = plotIn$Gene[i]
  mutant = plotIn$Mutant[i]
  expressionTemp=expression%>%
    select(gene)
  fitnesstemp = fitnessScores%>%
    select(mutant)
  # p=cbind(expressionTemp, fitnesstemp) %>%
  #   rownames_to_column('mergeCol')%>%
  #   left_join(fitnessMeta, by = 'mergeCol')%>%
  #   ggplot(aes(x = (.data[[gene]]),
  #              y = .data[[mutant]],
  #              group = tissue,
  #              col = day))+
  #   geom_point() +
  #   geom_smooth(method = 'lm')+
  #   labs(x = paste0('log2', gene), y = mutant)
  # plot(p)

  padj = significant_results%>%
    filter(Gene == gene & Mutant == mutant)%>%
    .$padjust

  rho = significant_results%>%
    filter(Gene == gene & Mutant == mutant)%>%
    .$Spearman_Correlation

  p=cbind(expressionTemp, fitnesstemp) %>%
    rownames_to_column('mergeCol')%>%
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
  ggsave(filename = paste0('djGeneExpressionFitnessScoresSpearmanIntegration/sigRes/',
                           mutant,'_',
                           gene,'.pdf'),
         units = 'in',
         height = 4.5,
         width = 4.5,
         p )
}
#write_tsv(plotIn, file = 'top10VstRnaMutantSpearmanCorrelations.tsv')

significant_results%>%
  group_by(Gene)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))%>%
  view()

genesWithSignificantAssociations=significant_results%>%
  group_by(Gene)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))%>%
  .$Gene

significant_results%>%
  group_by(Mutant)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))

mutants = significant_results%>%
  group_by(Mutant)%>%
  summarise(nAssociations = n())%>%
  filter(nAssociations > 9)%>%
  arrange(desc(nAssociations))%>%
  select(Mutant)%>%
  distinct()%>%
  .$Mutant

# for (m in mutants){
#   print('-----------------')
#   print(m)
#   g=significant_results%>%
#     filter(Mutant == m)%>%
#     .$Gene
#
#   enrichGO(gene = g,
#            universe = colnames(expression),
#            keyType = 'SYMBOL',
#            org.Mm.eg.db)%>%
#     as.data.frame()%>%
#     print()
#   print('-----------------')
#
# }

significant_results%>%
  pivot_wider(names_from = 'Gene',
              values_from = 'Spearman_Correlation',
              id_cols = 'Mutant', values_fill = 0)%>%
  column_to_rownames('Mutant')%>%
  ComplexHeatmap::Heatmap()




interestingGenes=significant_results%>%
  group_by(Gene)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))%>%
  filter(`number of associations`>=10)

significant_results_annotated=significant_results%>%
  filter(gene %in% interestingGenes$Gene)%>%
  rename(locusId=Mutant)%>%
  left_join(keggs, by = 'locusId')

'Need to redo this such that it takes into account the directionality of the
correlation'
for (g in interestingGenes$Gene){
  temp = significant_results_annotated%>%
    filter(Gene == g)

  temp = annotate_modules(temp)
  p=temp %>%
    mutate(moduleLevel1 = str_remove(moduleLevel1, "^\\s*[[:digit:]]+\\s+"), # Remove leading spaces & number
           moduleLevel2 = str_remove(moduleLevel2, "^\\s*[[:digit:]]+\\s+"),
           moduleLevel1 = wrap_format(20)(moduleLevel1),  # Wrap text for readability
           moduleLevel1 = fct_reorder(moduleLevel1, .x = moduleLevel1, .fun = function(x) length(x), .desc = TRUE)) %>%  # Order by count
    ggplot(aes(x = moduleLevel1, fill = moduleLevel2)) +
    geom_bar() +
    coord_flip() +  # Flip bars for better readability
    theme_minimal()+
    labs(x = 'Count',
         y = 'Brite level 1',
         fill = 'Brite level 2')+
    labs(title = paste('Mutants assoicated with', g, sep = ' '))
  plot(p)
}

plotIn = significant_results%>%
  filter(Gene == 'Reg3g')
for (i in 1:nrow(plotIn)){
  gene = plotIn$Gene[i]
  mutant = plotIn$Mutant[i]
  expressionTemp=expression%>%
    select(gene)
  fitnesstemp = fitnessScores%>%
    select(mutant)
  # p=cbind(expressionTemp, fitnesstemp) %>%
  #   rownames_to_column('mergeCol')%>%
  #   left_join(fitnessMeta, by = 'mergeCol')%>%
  #   ggplot(aes(x = (.data[[gene]]),
  #              y = .data[[mutant]],
  #              group = tissue,
  #              col = day))+
  #   geom_point() +
  #   geom_smooth(method = 'lm')+
  #   labs(x = paste0('log2', gene), y = mutant)
  # plot(p)

  p=cbind(expressionTemp, fitnesstemp) %>%
    rownames_to_column('mergeCol')%>%
    left_join(fitnessMeta, by = 'mergeCol')%>%
    ggplot(aes(y = (.data[[gene]]),
               x = .data[[mutant]],
               group = tissue,
               col = day))+
    geom_point() +
    geom_smooth(method = 'lm')+
    labs(y = paste0('log2', gene), x= mutant)
  plot(p)

}
