library(tidyverse)
library(parallel)

compute_correlation <- function(pair, method) {
  gene <- pair[[1]]
  mutant <- pair[[2]]
  test <- cor.test(expression[, gene], fitnessScores[, mutant], method = method)
  return(data.frame(Gene = gene,
                    Mutant = mutant,
                    Spearman_Correlation = test$estimate,
                    P_Value = test$p.value))
}

expression = read_tsv('sigGenesCounts.tsv')
fitnessMeta = read_tsv('fullbarseqMeta.txt')
sigMutants = read_tsv('barseqAdjustedParams/strong.tab')

sigMutants$name <- sub("setA", "", sigMutants$name)
sigMutants$name<- sub("_.*", "", sigMutants$name)
sigMutants$name <- sub("CO$", "Co", sigMutants$name)
sigMutants$name <- sub("DJ$", "Dj", sigMutants$name)

mutantsToTest=sigMutants %>%
  rename(sample = name)%>%
  left_join(fitnessMeta, by = 'sample')%>%
  filter(lrn >2)%>%
  group_by(locusId)%>%
  summarise('observed'=n())%>%
  filter(observed > 1)%>%
  .$locusId


metadataRNA = read.csv('metadataRnaSeq.csv')%>%
  as.data.frame()

sampleLibraryRna = metadataRNA %>%
  filter(Tissue == 'Colon'& Treatment == 'Bar-seq')%>%
  filter(sample != '2493Co')%>%
  select(sample, library)
'Clean up column names'
expression=expression %>%
  rename_with(~ sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", .))


expression=expression %>%
  column_to_rownames('gene')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('library')%>%
  merge(sampleLibraryRna, by = 'library') %>%
  select(-library)%>%
  column_to_rownames('sample')

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
num_cores <- 2
cor_results <- do.call(rbind, mclapply(1:nrow(gene_mutant_pairs), function(i) {
  compute_correlation(gene_mutant_pairs[i, ], method = 'spearman')
}, mc.cores = num_cores))

# Adjust p-values
cor_results$padjust <- p.adjust(cor_results$P_Value, method = 'fdr')

# Filter significant results
significant_results <- cor_results %>%
  filter(padjust < 0.05)

write_tsv(file = 'sigGeneStrongMutantSpearmanCorrelation.tsv', x = cor_results)
?write_tsv()

