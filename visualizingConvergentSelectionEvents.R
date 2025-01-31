library(tidyverse)
library(ggVennDiagram)
sigGenes=read_tsv('significantRepatedSelectionPhenotypes.tsv')

geneAnnotations = read_tsv('genesWithAnvioAnnotations.tsv')

keggs = geneAnnotations%>%
  select(locusId, kofamFunction,kofamAccession)%>%
  filter(is.na(kofamAccession)!= T)%>%
  distinct()

sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct()%>%
  view()

positiveSelectionDJDay1 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday1',
         direction == 'positive')%>%
  .$locusId

positiveSelectionColonDay1 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday1',
         direction == 'positive')%>%
  .$locusId

positiveSelectionDJDay3 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday3',
         direction == 'positive')%>%
  .$locusId

positiveSelectionColonDay3 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday3',
         direction == 'positive') %>%
  .$locusId

positiveSelectionDJDay7 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday7',
         direction == 'positive')%>%
  .$locusId

positiveSelectionColonDay7 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday7',
         direction == 'positive')%>%
  .$locusId

positiveSelectionDJDay14 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday14',
         direction == 'positive')%>%
  .$locusId

positiveSelectionColonDay14 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday14',
         direction == 'positive')%>%
  .$locusId

vennIN = list('colonDay1'=positiveSelectionColonDay1,
              'colonDay3'= positiveSelectionColonDay3,
              'colonDay7' = positiveSelectionColonDay7,
              'colonDay14' = positiveSelectionColonDay7,
              'djDay1' = positiveSelectionDJDay1,
              'djDay3' = positiveSelectionDJDay3,
              'djDay7' = positiveSelectionDJDay7,
              'djDay14' = positiveSelectionDJDay14)

ggVennDiagram(vennIN)

negativeSelectionDJDay1 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday1',
         direction == 'negative')%>%
  .$locusId

negativeSelectionColonDay1 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday1',
         direction == 'negative')%>%
  .$locusId

negativeSelectionDJDay3 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday3',
         direction == 'negative')%>%
  .$locusId

negativeSelectionColonDay3 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday3',
         direction == 'negative') %>%
  .$locusId

negativeSelectionDJDay7 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday7',
         direction == 'negative')%>%
  .$locusId

negativeSelectionColonDay7 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday7',
         direction == 'negative')%>%
  .$locusId

negativeSelectionDJDay14 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'djday14',
         direction == 'negative')%>%
  .$locusId

negativeSelectionColonDay14 = sigGenes %>%
  rename('locusId' = gene)%>%
  merge(keggs, by ='locusId')%>%
  distinct() %>%
  filter(group == 'colonday14',
         direction == 'negative')%>%
  .$locusId

vennInNegative = list('colonDay1'=negativeSelectionColonDay1,
              'colonDay3'= negativeSelectionColonDay3,
              'colonDay7' = negativeSelectionColonDay7,
              'colonDay14' = negativeSelectionColonDay7,
              'djDay1' = negativeSelectionDJDay1,
              'djDay3' = negativeSelectionDJDay3,
              'djDay7' = negativeSelectionDJDay7,
              'djDay14' = negativeSelectionDJDay14)

ggVennDiagram(vennInNegative)
'Not really sure if this makes the most sense to do'
' Day 14 dj has the most unique significant negative selection events, lets see what those are'

negativeSelectionDjDay14DF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group =='djday14')

filterDF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group !='djday14')

sigUniqueDay14DjNegative=negativeSelectionDjDay14DF %>%
  filter(!gene %in% filterDF$gene)%>%
  rename('locusId'= gene)%>%
  distinct()%>%
  merge(keggs, by = 'locusId')
'Day 3 dj has the second most'

negativeSelectionDjDay3DF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group =='djday3')

filterDF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group !='djday3')

sigUniqueDay3DjNegative=negativeSelectionDjDay3DF %>%
  filter(!gene %in% filterDF$gene)%>%
  rename('locusId'= gene)%>%
  distinct()%>%
  merge(keggs, by = 'locusId')

#day 7 dj negative unique
negativeSelectionDjDay7DF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group =='djday7')

filterDF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group !='djday7')

sigUniqueDay7DjNegative=negativeSelectionDjDay7DF %>%
  filter(!gene %in% filterDF$gene)%>%
  rename('locusId'= gene)%>%
  distinct()%>%
  merge(keggs, by = 'locusId')

#day 1 dh negative unique
negativeSelectionDjDay1DF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group =='djday1')

filterDF=sigGenes %>%
  filter(direction == 'negative')%>%
  filter(group !='djday1')

sigUniqueDay1DjNegative=negativeSelectionDjDay1DF %>%
  filter(!gene %in% filterDF$gene)%>%
  rename('locusId'= gene)%>%
  distinct()%>%
  merge(keggs, by = 'locusId')\

#unique dj
negativeSelectionDjDF=sigGenes %>%
filter(direction == 'negative')%>%
filter(grepl(pattern = 'dj',group)==T)

filterDF=sigGenes %>%
filter(direction == 'negative')%>%
filter(grepl(pattern = 'dj',group)==F)

sigUniqueDjNegative=negativeSelectionDjDF %>%
filter(!gene %in% filterDF$gene)%>%
rename('locusId'= gene)%>%
select(-c(group, pvalue, padj))%>%
distinct()%>%
merge(keggs, by = 'locusId')