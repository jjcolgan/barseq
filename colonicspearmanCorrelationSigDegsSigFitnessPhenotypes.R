library(parallel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)
library(progressr)
library(tidyverse)
library(KEGGREST)
library(scales)

'Using the DEGs with the GF mice included produces very strange results, this set is not worth using'

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

expression = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve RNA seq/2025RnaAnalysis/DeSeq2/barseqColonicOutputs/lrtResults/sigGeneTab.tsv')
fitnessMeta = read_tsv('fullbarseqMeta.txt')
sigMutants=read_tsv('linear models/mixedEffectsColonOnly/sigRes.tsv')
'maaslin2res'
# numericDayRes = read_tsv('colonLogRatiosMaaslin2NumericDay/significant_results.tsv')
# categoricalDayRes = read_tsv('colonLogRatiosMaaslin2/significant_results.tsv')
#
# numericDayRes=numericDayRes%>%
#   filter(metadata ==  'dayNumeric',
#          qval < .1)%>%
#   .$feature
#
# categoricalDayRes=categoricalDayRes%>%
#   filter(metadata ==  'day',
#          qval < .1)%>%
#   .$feature%>%
#   unique()

sigMutants = unique(sigMutants$locusId)
annotations = read_tsv('genesWithAnvioAnnotations.tsv')

keggs = annotations %>%
  select(locusId, kofamAccession, kofamFunction)%>%
  distinct()

mutantsToTest=sigMutants
print(length(mutantsToTest))

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
  filter(locusId %in% mutantsToTest)%>%
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

write_tsv(cor_results, 'colonicRNAFitnesssIntegretion/spearmanRes/spearmanRes.tsv')

# Filter significant results
significant_results <- cor_results %>%
  filter(padjust < 0.05)
cor_results%>%
  ggplot(aes(x = P_Value))+
  geom_histogram(bins = 100)


cor_results%>%
  arrange(desc(-log10(P_Value)))%>%
  view()


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
  head(100)

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

  ggsave(filename = paste0('colonicGeneExpressionFitnessScoresSpearmanIntegration/sigRes/',
                           mutant,'_',
                           gene,'.pdf'),
         units = 'in',
         height = 4.5,
         width = 4.5,
         p )
}

plotIn=significant_results%>%
  arrange(desc(-log10(padjust)), Spearman_Correlation)%>%
  tail(100)

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

  ggsave(filename = paste0('colonicGeneExpressionFitnessScoresSpearmanIntegration/sigRes/',
                           mutant,'_',
                           gene,'.pdf'),
         units = 'in',
         height = 4.5,
         width = 4.5,
         p )
}

significant_results%>%
  group_by(Gene)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))%>%
  view()

significant_results%>%
  filter(Spearman_Correlation > 0)%>%
  nrow()

significant_results%>%
  filter(Spearman_Correlation < 0)%>%
  nrow()

genesWithSignificantAssociations=significant_results%>%
  group_by(Gene)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))%>%
  .$Gene

significant_results%>%
  group_by(Mutant)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))

significant_results%>%
  group_by(Mutant)%>%
  summarise('number of associations' = n())%>%
  arrange(desc(`number of associations`))%>%
  ggplot(aes(x = `number of associations`))+
  geom_histogram()

mutants = significant_results%>%
  group_by(Mutant)%>%
  summarise(nAssoc = n())%>%
  filter(nAssoc > 9)%>%
  arrange(desc(nAssoc))%>%
  distinct()%>%
  .$Mutant

for (m in mutants[1:20]){
  print('-----------------')
  print(m)
  print('inverse relationship genes')
  g=significant_results%>%
    filter(Mutant == m & Spearman_Correlation < 0 )%>%
    .$Gene

  enrichGO(gene = g,
           keyType = 'SYMBOL',
           org.Mm.eg.db)%>%
    as.data.frame()%>%
    filter(Count > 2)%>%
    print()

  print('positve relationship genes')
  g=significant_results%>%
    filter(Mutant == m & Spearman_Correlation > 0 )%>%
    .$Gene

  enrichGO(gene = g,
           keyType = 'SYMBOL',
           org.Mm.eg.db)%>%
    as.data.frame()%>%
    filter(Count > 2)%>%
    print()

  print('-----------------')

}

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

significant_results%>%
  group_by(Gene)%>%
  summarise('number of associations' = n())%>%
  ggplot(aes(x = `number of associations`))+
  geom_histogram(bins = 100)

significant_results_annotated=significant_results%>%
  filter(gene %in% interestingGenes$Gene)%>%
  rename(locusId=Mutant)%>%
  left_join(keggs, by = 'locusId')

'Need to redo this such that it takes into account the directionality of the
correlation'
for (g in interestingGenes$Gene[1:5]){
  temp = significant_results_annotated%>%
    filter(Gene == g & Spearman_Correlation > 0)

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
    labs(title = paste('Mutants with positive assoication with', g, sep = ' '))
  plot(p)

  temp = significant_results_annotated%>%
    filter(Gene == g & Spearman_Correlation < 0)

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
    labs(title = paste('Mutants with positive negative with', g, sep = ' '))
  plot(p)
}

plotIn = significant_results%>%
    filter(Gene == 'Gp2')
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