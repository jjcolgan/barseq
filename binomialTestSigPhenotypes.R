'The idea here is to use a binomial test to verify that the significant phenotypes are indeed assoicated with
a given region and time point. '

library(tidyverse)



strongFitnessEffects=read_tsv('barseqAdjustedParams/strong.tab')
strongFitnessEffects %>%
  select(locusId)%>%
  distinct()%>%
  nrow()
metadata = read_tsv('fullbarseqMeta.txt')
metadata = metadata %>%
  mutate('tissueDayMouse'=paste(tissueDay, mouse, sep = '-'))
quality=read_tsv('barseqAdjustedParams/fit_quality.tab')
quality$name <- sub("setA", "", quality$name)
quality$name <- sub("_.*", "", quality$name)
quality$name <- sub("CO$", "Co", quality$name)
quality$name <- sub("DJ$", "Dj", quality$name)
metadata = quality %>%
  rename('sample'= name)%>%
  left_join(metadata) %>%
  filter(cor12 >= .2)

# clean up sample names
strongFitnessEffects$name <- sub("setA", "", strongFitnessEffects$name)
strongFitnessEffects$name <- sub("_.*", "", strongFitnessEffects$name)
strongFitnessEffects$name <- sub("CO$", "Co", strongFitnessEffects$name)
strongFitnessEffects$name <- sub("DJ$", "Dj", strongFitnessEffects$name)
strongFitnessEffects%>%
  select(name)%>%
  distinct()%>%
  nrow()

strongFitnessEffects=strongFitnessEffects %>%
  rename('sample'= name)

strongFitnessEffects=strongFitnessEffects%>%
  left_join(metadata, by = 'sample')

'Count the number of genes annotated in barseq'
gene=read_tsv(file = 'barseqAdjustedParams/genes')
numberGenes=gene %>%
  nrow()

strongPhenotypes= strongFitnessEffects %>%
  select(locusId)%>%
  distinct()%>%
  .$locusId

testInput=strongFitnessEffects %>%
  filter(tissue =='colon',
         day == 'day1')
testInput%>%
  select(sample)%>%
  distinct()%>%
  nrow()

groups = strongFitnessEffects %>%
  select(tissueDay)%>%
  filter(is.na(tissueDay)!= T)%>%
  distinct()%>%
  .$tissueDay
groups
selectionRes = data.frame('gene' = character(),
                               'group' = character(),
                               'pvalue' = double(),
                               'direction' = character())
# test negative selection
for (g in 1:length(groups)){
  inputMat = strongFitnessEffects %>%
    filter(tissueDay == groups[g])%>%
    filter(lrn < 0)
  nSamples = strongFitnessEffects %>%
    filter(tissueDay == groups[g])%>%
    select(sample)%>%
    distinct()%>%
    nrow()
  for (i in 1:length(strongPhenotypes) ){
    sucesses = inputMat %>%
      filter(locusId == strongPhenotypes[i])%>%
      nrow()
    trials = nSamples
    binomRes = (binom.test(x = sucesses, n = trials, p = 1 / numberGenes))
    selectionRes = rbind(negativeSelection,
                              data.frame('gene' = strongPhenotypes[i],
                                         'group' = groups[g],
                                         'pvalue' = binomRes$p.value,
                                         'direction' = 'negative'))

  }
}

# test positive selection
for (g in 1:length(groups)){
  inputMat = strongFitnessEffects %>%
    filter(tissueDay == groups[g])%>%
    filter(lrn > 0)
  nSamples = strongFitnessEffects %>%
    filter(tissueDay == groups[g])%>%
    select(sample)%>%
    distinct()%>%
    nrow()
  for (i in 1:length(strongPhenotypes) ){
    sucesses = inputMat %>%
      filter(locusId == strongPhenotypes[i])%>%
      nrow()
    trials = nSamples
    binomRes = (binom.test(x = sucesses, n = trials, p = 1 / numberGenes))
    selectionRes = rbind(negativeSelection,
                              data.frame('gene' = strongPhenotypes[i],
                                         'group' = groups[g],
                                         'pvalue' = binomRes$p.value,
                                         'direction' = 'positive'))

  }
}

selectionRes$padj = p.adjust(selectionRes$pvalue, method = 'fdr')

sigRes = selectionRes %>%
  filter(padj < .05)

sigRes %>%
  group_by(group, direction)%>%
  summarise(n())
gene=read_tsv(file = 'barseqAdjustedParams/genes')

write_tsv(file = 'significantRepatedSelectionPhenotypes.tsv', sigRes)