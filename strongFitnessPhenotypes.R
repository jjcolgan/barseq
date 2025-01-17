library(tidyverse)
library(ggpubr)
library(Maaslin2)
library(ComplexHeatmap)

strongFitnessTab=read_tsv('full/strong.tab')
metadata = read_tsv('fullbarseqMeta.txt')
pcaIn<-strongFitnessTab %>%
  pivot_wider(names_from = locusId,
              id_cols = name,
              values_from = t)%>%
  column_to_rownames('name')
pcaIn[is.na(pcaIn)]<-0
write_tsv(rownames_to_column(pcaIn, 'name'), file = 'full/strongFitnessMatrix0Imputed.tsv')
pcaIn<-read_tsv('full/strongFitnessMatrix0Imputed.tsv')%>%
  column_to_rownames('name')
strongFitnessScores<-pcaIn%>%
  prcomp(center = TRUE)
summary(strongFitnessScores)
strongFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  merge(metadata, by = 'sample') %>%
  ggplot(aes(x = PC1,
             col = tissue,
             shape = day,
             y = PC2))+
  geom_point()

pcaIn %>%
  t()%>%
  Heatmap(show_row_names = FALSE)

tScoresWithMetaLong<-pcaIn %>%
  rownames_to_column('sample')%>%
  pivot_longer(2:345,
               names_to = 'locusId',
               values_to = 't') %>%
  merge(metadata, by = 'sample')

siStrongFitnessScores<-tScoresWithMetaLong %>%
  filter(tissue == 'dj')%>%
  pivot_wider(id_cols = sample,
              names_from = locusId,
              values_from = t) %>%
  column_to_rownames('sample') %>%
  prcomp(center = TRUE)
summary(siStrongFitnessScores)

siStrongFitnessScores$x %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  merge(metadata,
        by = "sample") %>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = day))+
  geom_point()

coStrongFitnessScores<-tScoresWithMetaLong %>%
  filter(tissue == 'colon')%>%
  pivot_wider(id_cols = sample,
              names_from = locusId,
              values_from = t) %>%
  column_to_rownames('sample') %>%
  prcomp(center = TRUE)
summary(coStrongFitnessScores)

coStrongFitnessScores$x %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  merge(metadata,
        by = "sample") %>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = day))+
  geom_point()