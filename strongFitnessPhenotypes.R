library(tidyverse)
library(ggpubr)
library(Maaslin2)
library(ComplexHeatmap)

strongFitnessTab=read_tsv('full/strong.tab')
metadata = read_tsv('fullbarseqMeta.txt')
metadata$numericDay = 0
metadata$numericDay[metadata$day== 'day1'] = 1
metadata$numericDay[metadata$day== 'day3'] = 3
metadata$numericDay[metadata$day== 'day 7'] = 7
metadata$numericDay[metadata$day== 'day 14'] = 14
metadata=metadata%>%
  mutate(mouseDayTissue = paste(mouse, day, tissue,sep = '-'))
pcaIn<-strongFitnessTab %>%
  pivot_wider(names_from = locusId,
              id_cols = name,
              values_from = t)%>%
  column_to_rownames('name')
pcaIn[is.na(pcaIn)]<-0
#write_tsv(rownames_to_column(pcaIn, 'name'), file = 'full/strongFitnessMatrix0Imputed.tsv')
pcaIn<-read_tsv('full/strongFitnessMatrix0Imputed.tsv')%>%
  column_to_rownames('name')
strongFitnessScores<-pcaIn%>%
  prcomp(center = TRUE)
summary(strongFitnessScores)
strongFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = tissue,
             shape = day,
             y = PC2))+
  geom_point()

strongFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             col = numericDay,
             y = PC2))+
  geom_point()

strongFitnessScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  view()
  ggplot(aes(x = PC1,
             col = mouse,
             group = mouse,
             y = PC2))+
  geom_point()+
  geom_line()

pcaIn %>%
  rownames_to_column('sample') %>%
  pivot_longer(cols = c(2:345), names_to = 'gene', values_to = 'fitnessScore')%>%
  left_join(metadata,by = 'sample')%>%
  pivot_wider(names_from = mouseDayTissue, id_cols = gene, values_from = fitnessScore)%>%
  column_to_rownames('gene')%>%
  Heatmap(show_row_names = FALSE)

pcaIn %>%
  rownames_to_column('sample') %>%
  pivot_longer(cols = c(2:345), names_to = 'gene', values_to = 'fitnessScore')%>%
  left_join(metadata,by = 'sample')%>%
  filter(tissue != 'colon') %>%
  pivot_wider(names_from = mouseDayTissue, id_cols = gene, values_from = fitnessScore)%>%
  column_to_rownames('gene')%>%
  Heatmap(show_row_names = FALSE)

pcaIn %>%
  rownames_to_column('sample') %>%
  pivot_longer(cols = c(2:345), names_to = 'gene', values_to = 'fitnessScore')%>%
  left_join(metadata,by = 'sample')%>%
  filter(tissue == 'colon') %>%
  pivot_wider(names_from = mouseDayTissue, id_cols = gene, values_from = fitnessScore)%>%
  column_to_rownames('gene')%>%
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