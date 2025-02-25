library(tidyverse)

fitness = read_tsv('barseqAdjustedParams/fit_t.tab')
colnames(fitness) <- sub("setA", "", colnames(fitness))
colnames(fitness) <- sub("_.*", "", colnames(fitness))
colnames(fitness) <- sub("CO$", "Co", colnames(fitness))
colnames(fitness) <- sub("DJ$", "Dj", colnames(fitness))

metadata = read_tsv('fullbarseqMeta.txt')
metadata$cage = as.factor(metadata$cage)
metadata$lane = as.factor(metadata$lane)
quality=read_tsv('barseqAdjustedParams/fit_quality.tab')
quality$name <- sub("setA", "", quality$name)
quality$name <- sub("_.*", "", quality$name)
quality$name <- sub("CO$", "Co", quality$name)
quality$name <- sub("DJ$", "Dj", quality$name)
quality = quality %>%
  rename('sample'= name)

metadata = metadata %>%
  left_join(quality, by = 'sample')

pcaOut=fitness %>%
  select(-c(desc,
            sysName))%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% metadata$sample)%>%
  column_to_rownames('sample')%>%
  prcomp(center = TRUE)

summary(pcaOut)
dir.create('qcPcas')
dir.create(path= 'qcPcas/tscoresNoFilterPCAs')
dir.create(path= 'qcPcas/tscoresNoFilterPCAs/pc2Pc3')
columns= colnames(metadata)
for (c in 1:length(columns)){
  p= pcaOut$x %>%
    as.data.frame()%>%
    rownames_to_column('sample')%>%
    left_join(metadata, by = 'sample')%>%
    ggplot(aes(x = PC2,
               y = PC3,
               col = .data[[columns[c]]]))+
    geom_point()+
    labs(x = 'PC2 23.13% var',
         y = 'PC3 7.626% var')
  path = paste0('qcPcas/tscoresNoFilterPCAs/pc2Pc3/', columns[c],'.pdf')
  ggsave(plot = p,
         filename = path,
         units = c('in'),
         width = 4,
         height = 4)

}

pcaOut=fitness %>%
  select(-c(desc,
            sysName))%>%
  column_to_rownames('locusId')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% metadata$sample)%>%
  column_to_rownames('sample')%>%
  prcomp(center = TRUE, scale. = T)

summary(pcaOut)
dir.create('qcPcas')
dir.create(path= 'qcPcas/tscoresNoFilterPCAsScaled')
dir.create(path= 'qcPcas/tscoresNoFilterPCAsScaled/pc2Pc3')
columns= colnames(metadata)
for (c in 1:length(columns)){
  p= pcaOut$x %>%
    as.data.frame()%>%
    rownames_to_column('sample')%>%
    left_join(metadata, by = 'sample')%>%
    ggplot(aes(x = PC2,
               y = PC3,
               col = .data[[columns[c]]]))+
    geom_point()+
    labs(x = 'PC2 22.13% var',
         y = 'PC3 13.21% var')
  path = paste0('qcPcas/tscoresNoFilterPCAsScaled/pc2Pc3/', columns[c],'.pdf')
  ggsave(plot = p,
         filename = path,
         units = c('in'),
         width = 4,
         height = 4)

}


pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissue))+
  geom_point()+
  labs(x = 'PC1 - 25.43%',
       y = 'PC2 - 22.75%')

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = as.factor(lane)))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             shape = tissue,
             col = day))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = day,
             group = mouse))+

  geom_point()+
  geom_line()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = as.factor(cage),
             group = mouse))+

  geom_point()+
  geom_line()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = as.factor(cage),
             shape = day))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissue,
             group = mouse))+

  geom_point()+
  geom_line()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = log10(millionBases)))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = meanQualityScore ))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = cor12 ))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = gccor ))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = `>=Q30 bases` ))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = mad12 ))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = adjcor ))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = u ))+

  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = maxFit))+
  geom_point()

pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = opcor))+
  geom_point()


dir.create(path= 'qcPcas/noFilterPCAs/pc2Pc3')
columns= colnames(metadata)
for (c in 1:length(columns)){
  p= pcaOut$x %>%
    as.data.frame()%>%
    rownames_to_column('sample')%>%
    left_join(metadata, by = 'sample')%>%
    ggplot(aes(x = PC2,
               y = PC3,
               col = .data[[columns[c]]]))+
    geom_point()+
    labs(x = 'PC2 - 22.75%',
         y = 'PC3 - 7.53%')
  path = paste0('qcPcas/noFilterPCAs/pc2Pc3/', columns[c],'.pdf')
  ggsave(plot = p,
         filename = path,
         units = c('in'),
         width = 4,
         height = 4)

}