library(tidyverse)

'I am not sure how to use binom test, but I can for sure use a fishers exact
test'

'I think just filtering for the significant and strong selective events unique to a tissue, then
ranking by the number of times seen and the score would probably be the best way to do this.
report like the top 25 or something. '
metadata = read_tsv('fullbarseqMeta.txt')

fitness = read_tsv('barseqAdjustedParams/fit_logratios.tab')
colnames(fitness) <- sub("setA", "", colnames(fitness))
colnames(fitness) <- sub("_.*", "", colnames(fitness))
colnames(fitness) <- sub("CO$", "Co", colnames(fitness))
colnames(fitness) <- sub("DJ$", "Dj", colnames(fitness))

fitnessLong = fitness %>%
  pivot_longer(cols = c(4:ncol(fitness)),
               names_to = 'sample',
               values_to = 'fitnessScore')

tscores = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/barseqAdjustedParams/fit_t.tab')
colnames(tscores) <- sub("setA", "", colnames(tscores))
colnames(tscores) <- sub("_.*", "", colnames(tscores))
colnames(tscores) <- sub("CO$", "Co", colnames(tscores))
colnames(tscores) <- sub("DJ$", "Dj", colnames(tscores))

tscoresLong = tscores %>%
  pivot_longer(cols = c(4:ncol(tscores)),
               names_to = 'sample',
               values_to = 'tscore')

mergedTandFitnessScoresLong=tscoresLong %>%
  left_join(fitnessLong,
            by = c('locusId',
                   'desc',
                   'sysName',
                   'sample'))%>%
  left_join(metadata, by = 'sample')

sigDay1=mergedTandFitnessScoresLong %>%
  filter(day == 'day1')%>%
  filter(abs(tscore) >4,
         abs(fitnessScore)> 0)

sigDay1DjLoci = sigDay1 %>%
  filter(tissue == 'dj')

sigDay1CoLoci = sigDay1 %>%
  filter(tissue == 'colon')


sigDay1 %>%
  filter(!locusId %in% sigDay1CoLoci$locusId)%>%
  view()



sigDay1 %>%
  group_by(tissue, locusId)%>%
  summarise(meanFitnessScore = mean(fitnessScore))%>%
  pivot_wider(id_cols = 'locusId', names_from = 'tissue', values_from = 'meanFitnessScore', values_fill = 0)%>%
  ggplot(aes(x = colon,
             y = dj))+
  geom_point()+
  xlim(-10,10)+
  ylim(-10,10)

sigDay3=mergedTandFitnessScoresLong %>%
  filter(day == 'day3')%>%
  filter(abs(tscore) >4)%>%
  filter(abs(fitnessScore) > .5)

sigDay3%>%
  group_by(tissue, locusId)%>%
  summarise(meanFitnessScore = mean(fitnessScore))%>%
  pivot_wider(id_cols = 'locusId', names_from = 'tissue', values_from = 'meanFitnessScore', values_fill = 0)%>%
  ggplot(aes(x = colon,
             y = dj))+
  geom_point()+
  xlim(-10,10)+
  ylim(-10,10)


sigDay7=mergedTandFitnessScoresLong %>%
  filter(day == 'day7')%>%
  filter(abs(tscore) >4)%>%
  filter(abs(fitnessScore)>.5)

sigDay7%>%
  filter(tissue == 'dj')%>%
  filter(fitnessScore > 0)


sigDay7%>%
  group_by(tissue, locusId)%>%
  summarise(meanFitnessScore = mean(fitnessScore))%>%
  pivot_wider(id_cols = 'locusId', names_from = 'tissue', values_from = 'meanFitnessScore', values_fill = 0)%>%
  ggplot(aes(x = colon,
             y = dj))+
  geom_point()+
  xlim(-10,10)+
  ylim(-10,10)

sigDay14=mergedTandFitnessScoresLong %>%
  filter(day == 'day14')%>%
  filter(abs(tscore) >4)

sigDay14%>%
  group_by(tissue, locusId)%>%
  summarise(meanFitnessScore = mean(fitnessScore))%>%
  pivot_wider(id_cols = 'locusId', names_from = 'tissue', values_from = 'meanFitnessScore', values_fill = 0)%>%
  ggplot(aes(x = colon,
             y = dj))+
  geom_point()+
  xlim(-10,10)+
  ylim(-10,10)

mergedTandFitnessScoresLong %>%
  filter(day != 'T0')%>%
  filter(abs(tscore) >4)%>%
  filter(abs(fitnessScore) > .5)%>%
  group_by(tissue, locusId)%>%
  summarise(meanFitnessScore = mean(fitnessScore))%>%
  pivot_wider(id_cols = 'locusId', names_from = 'tissue', values_from = 'meanFitnessScore', values_fill = 0)%>%
  ggplot(aes(x = colon,
             y = dj))+
  geom_point()+
  xlim(-10,10)+
  ylim(-10,10)
