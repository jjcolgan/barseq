library(tidyverse)
library(KEGGREST)
getKeggBrite = function(geneVec){

  out=keggGet(geneVec$kofamAccession[1])
  gene = geneVec$kofamAccession[1]
  superpathway = out[[1]]$BRITE[2]
  subpathway=out[[1]]$BRITE[3]
  output = cbind(gene, superpathway, subpathway)
  for (i in 2:nrow(geneVec)){
    out=keggGet(geneVec$kofamAccession[i])
    gene = geneVec$kofamAccession[i]
    superpathway = out[[1]]$BRITE[2]
    subpathway=out[[1]]$BRITE[3]
    output = rbind(output,cbind(gene, superpathway, subpathway))

  }
  return(output)
}
'In general, it seems like very few gene deletions have positive effects in the colon'
fullgenesFit = read_tsv('fullGenesTab.tsv')

fitnessScoresMat=fullgenesFit %>%
  filter(tissue == 'colon')%>%
  select(sample, locusId, fitnessScore) %>%
  distinct()%>%
  pivot_wider(names_from = 'sample',id_cols = 'locusId', values_from = 'fitnessScore')

functions =fullgenesFit %>%
  select(locusId, name, desc,kofamAccession,kofamAccession, pfamAcession, pfamFunction) %>%
  distinct()

fitnessScoresMat%>%
  column_to_rownames('locusId')%>%
  Heatmap(show_row_names = F)

wss <- (nrow(fitnessScoresMat) - 1) * sum(apply(column_to_rownames(fitnessScoresMat,'locusId'), 2, var))
# Loop through cluster numbers from 2 to 20 to calculate WSS for each k-means solution
for (i in 2:20) {
  # Run k-means clustering with 'i' clusters and store the WSS for the solution
  wss[i] <- sum(kmeans(column_to_rownames(fitnessScoresMat,'locusId'),
                       centers = i)$withinss)
}
plot(1:20, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")
clusters=fitnessScoresMat%>%
  column_to_rownames('locusId')%>%
  kmeans(centers = 4)

clusterMembership=clusters$cluster%>%
  as.data.frame()

colnames(clusterMembership) = c('cluster')

clusterMembership=clusterMembership %>%
  rownames_to_column('locusId')

clusterMembership %>%
  group_by(cluster)%>%
  summarise('clusterSize'= n())

cluster1=clusterMembership %>%
  filter(cluster == 1) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  filter(tissue == 'colon')%>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster1$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  labs(title = 'Cluster 1 - colon negative, worse at later timepoints')+
  theme(legend.position = 'none')

cluster2=clusterMembership %>%
  filter(cluster == 2) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  filter(tissue == 'colon')%>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster2$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  theme(legend.position = 'none')+
  labs(title = 'cluster 2 - early neutral, late strong positive')

cluster3=clusterMembership %>%
  filter(cluster == 3) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  filter(tissue == 'colon')%>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster3$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+

  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  theme(legend.position = 'none')+
  labs(title = 'Cluster 3 - colon weak, but real weird on day 3')

cluster4=clusterMembership %>%
  filter(cluster == 4) %>%
  merge(functions, by = 'locusId')

fullgenesFit %>%
  filter(tissue == 'colon')%>%
  select(sample, locusId, fitnessScore, numericDay, tissue) %>%
  filter(locusId %in% cluster4$locusId,
         tissue != 'T0')%>%
  distinct() %>%
  group_by(tissue, locusId, numericDay)%>%
  summarise('geneMeanFitScore' = mean(fitnessScore))%>%
  ggplot(aes(x = numericDay,
             col = locusId,
             group = locusId,
             y = geneMeanFitScore,))+
  geom_point()+
  geom_line()+
  facet_wrap(~tissue)+
  theme(legend.position = 'none')+
  labs(title = 'Cluster 4 - colon, decreasing fitness later on')

cluster1Keggs=cluster1 %>%
  select(kofamAccession)%>%
  filter(is.na(kofamAccession)== FALSE) %>%
  distinct()

cluster1GeneCategories=getKeggBrite(cluster1Keggs)
cluster1GeneCategories=cluster1GeneCategories %>%
  as.data.frame()
cluster1GeneCategories$superpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster1GeneCategories$superpathway)
cluster1GeneCategories %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(superpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()+
  labs(title = 'Cluster 1 KO Brite superpathway')

cluster1GeneCategories$subpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster1GeneCategories$subpathway)
cluster1GeneCategories %>%
  mutate(subpathwayCleaned = str_wrap(subpathwayCleaned, width = 10)) %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(subpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()

cluster2Keggs=cluster2 %>%
  select(kofamAccession)%>%
  filter(is.na(kofamAccession)== FALSE) %>%
  distinct()

cluster2GeneCategories=getKeggBrite(cluster2Keggs)
cluster2GeneCategories=cluster2GeneCategories %>%
  as.data.frame()
cluster2GeneCategories$superpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster2GeneCategories$superpathway)
cluster2GeneCategories %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(superpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()+
  labs(title = 'Cluster 2 KO Brite superpathway')
cluster2GeneCategories$subpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster2GeneCategories$subpathway)
cluster2GeneCategories %>%
  mutate(subpathwayCleaned = str_wrap(subpathwayCleaned, width = 10)) %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(subpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()

cluster3Keggs=cluster3 %>%
  select(kofamAccession)%>%
  filter(is.na(kofamAccession)== FALSE) %>%
  distinct()

cluster3GeneCategories=getKeggBrite(cluster3Keggs)
cluster3GeneCategories=cluster3GeneCategories %>%
  as.data.frame()
cluster3GeneCategories$superpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster3GeneCategories$superpathway)
cluster3GeneCategories %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(superpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()+
  labs(title = 'Cluster 2 KO Brite superpathway')
cluster3GeneCategories$subpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster3GeneCategories$subpathway)
cluster3GeneCategories %>%
  mutate(subpathwayCleaned = str_wrap(subpathwayCleaned, width = 10)) %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(subpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()

cluster4Keggs=cluster4 %>%
  select(kofamAccession)%>%
  filter(is.na(kofamAccession)== FALSE) %>%
  distinct()

cluster4GeneCategories=getKeggBrite(cluster4Keggs)
cluster4GeneCategories=cluster4GeneCategories %>%
  as.data.frame()
cluster4GeneCategories$superpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster4GeneCategories$superpathway)
cluster4GeneCategories %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(superpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()+
  labs(title = 'Cluster 4 KO Brite superpathway')
cluster4GeneCategories$subpathwayCleaned<- sub("^\\s+\\d+\\s+", "", cluster4GeneCategories$subpathway)
cluster4GeneCategories %>%
  mutate(subpathwayCleaned = str_wrap(subpathwayCleaned, width = 10)) %>%
  mutate(superpathwayCleaned = str_wrap(superpathwayCleaned, width = 10)) %>%
  ggplot(aes(x = fct_infreq(subpathwayCleaned),
             fill = superpathwayCleaned))+
  geom_bar()