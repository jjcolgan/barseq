library(tidyverse)
library(ggpubr)
qualityInfo=read_tsv('full/fit_quality.tab')
view(qualityInfo)
qualityInfo$tissue <- 'Co'
qualityInfo$tissue[grepl(qualityInfo$name, pattern = 'dj', ignore.case =T)] = 'Dj'
qualityInfo$tissue[grepl(qualityInfo$name, pattern = 't0', ignore.case =T)] = 't0'

qualityInfo$short = factor(qualityInfo$short, levels = c('day1',
                                                            'day3',
                                                            'day 7',
                                                            'day 14',
                                                            'Time0'))

qualityInfo%>%
  ggplot(aes(x = fUsedStrainsSeen,
             fill = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_histogram(binwidth = .1)+
  labs(title = 'Histogram % of genic strains in each tissue')

qualityInfo%>%
  filter(tissue != 't0')%>%
  ggplot(aes(y = fUsedStrainsSeen,
             x = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_point()+
  geom_smooth()


qualityInfo%>%
  ggplot(aes(x = fUsedStrainsSeen,
             fill = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_histogram(binwidth = .1)+
  labs(title = 'Histogram % of genic strains')

qualityInfo%>%
  filter(tissue != 't0') %>%
  group_by(tissue, short)%>%
  summarise('meanPercentageGenicStrains'=mean(fUsedStrainsSeen))%>%
  ggplot(aes(y = meanPercentageGenicStrains,
             x = short,
             fill = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_point()+
  geom_line()+
  labs(title = 'Mean percentage of genic strains at each time point')

qualityInfo%>%
  ggplot(aes(y = fUsedStrainsSeen,
             x = tissue))+
  geom_boxplot()+stat_compare_means(comparisons = list(c('Co','Dj'),
                                                       c('Dj', 't0'),
                                                       c('Co','t0')))+
  labs(title = 'Low percentage of genic strains in Dj indicate a population bottleneck')

qualityInfo%>%
  ggplot(aes(x = fNGSeen,
             fill = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_histogram(binwidth = .1)+
  labs(title = 'Histogram of the percentage of intergenic strains')

qualityInfo%>%
  ggplot(aes(y = fNGSeen,
             x = tissue))+
  geom_boxplot()+stat_compare_means(comparisons = list(c('Co','Dj'),
                                                       c('Dj', 't0'),
                                                       c('Co','t0')))+
  labs(title = 'Low percentage of intergenic strains in Dj indicate a population bottleneck')


qualityInfo %>%
  ggplot(aes(x = maxFit,
             fill = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_histogram(binwidth = 1)

qualityInfo %>%
  ggplot(aes(x = cor12,
             fill = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_histogram(binwidth = .1)

qualityInfo %>%
  filter(tissue != 't0')%>%
  ggplot(aes(x = mad12,
             fill = short))+
  facet_wrap(~tissue, ncol = 1)+
  geom_histogram(binwidth = .1)+
  geom_vline(xintercept= .5)

qualityInfo %>%
  filter(tissue != 't0')%>%
  ggplot(aes(y = cor12,
             x = tissue))+
  geom_boxplot()+
  stat_compare_means()


qualityInfo %>%
  filter(short == 'day 7',
         tissue =='Dj')%>%
  view()
