library (tidyverse)
day1PositiveModules = read_tsv('positiveDay1ModulesRanked.tsv')
day1NegativeModules = read_tsv('negativeDay1ModulesRanked.tsv')
day1DjPositiveModules = read_tsv('day1DjPositiveModulesRanked.tsv')
day1DjNegativeModules = read_tsv('day1DjNegativeModulesRanked.tsv')
day1ColonPositiveModules = read_tsv('day1ColonPositiveModulesRanked.tsv')
day1ColonNegativeModules = read_tsv('day1ColonNegativeModulesRanked.tsv')

day1PositiveModules$tissue = 'both'
day1NegativeModules$tissue = 'both'
day1PositiveModules$selection = 'positive'
day1NegativeModules$selection = 'negative'

day1DjPositiveModules$tissue = 'dj'
day1DjNegativeModules$tissue = 'dj'
day1DjPositiveModules$selection = 'positive'
day1DjNegativeModules$selection = 'negative'

day1ColonPositiveModules$tissue = 'colon'
day1ColonNegativeModules$tissue = 'colon'
day1ColonPositiveModules$selection = 'positive'
day1ColonNegativeModules$selection = 'negative'

day1PositiveModules=day1PositiveModules %>%
  filter(!kofamAccession%in% day1DjPositiveModules$kofamAccession)%>%
  filter(!kofamAccession%in% day1ColonPositiveModules$kofamAccession)

day1NegativeModules=day1NegativeModules %>%
  filter(!kofamAccession%in% day1DjNegativeModules$kofamAccession)%>%
  filter(!kofamAccession%in% day1ColonNegativeModules$kofamAccession)


rbind(day1PositiveModules,
      day1NegativeModules,
      day1DjPositiveModules,
      day1DjNegativeModules,
      day1ColonPositiveModules,
      day1ColonNegativeModules)%>%
  na.omit()%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 1 significant phenotypes')

rbind(day1PositiveModules,
      day1NegativeModules,
      day1DjPositiveModules,
      day1DjNegativeModules,
      day1ColonPositiveModules,
      day1ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness) > 2)%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 1 strong phenotypes')

rbind(day1PositiveModules,
      day1NegativeModules,
      day1DjPositiveModules,
      day1DjNegativeModules,
      day1ColonPositiveModules,
      day1ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 1 Metabolism level 2')

rbind(day1PositiveModules,
      day1NegativeModules,
      day1DjPositiveModules,
      day1DjNegativeModules,
      day1ColonPositiveModules,
      day1ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness) > 2)%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 1 Metabolism level 2 strong phenotypes')

rbind(day1PositiveModules,
      day1NegativeModules,
      day1DjPositiveModules,
      day1DjNegativeModules,
      day1ColonPositiveModules,
      day1ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 1 Brite Hierarchy level 2 ')

rbind(day1PositiveModules,
      day1NegativeModules,
      day1DjPositiveModules,
      day1DjNegativeModules,
      day1ColonPositiveModules,
      day1ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness) > 2)%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 1 Brite Hierarchy Strong Phenotypes')

day3PositiveModules = read_tsv('positiveDay3ModulesRanked.tsv')
day3NegativeModules = read_tsv('negativeDay3ModulesRanked.tsv')
day3DjPositiveModules = read_tsv('day3DjPositiveModulesRanked.tsv')
day3DjNegativeModules = read_tsv('day3DjNegativeModulesRanked.tsv')
day3ColonPositiveModules = read_tsv('day3ColonPositiveModulesRanked.tsv')
day3ColonNegativeModules = read_tsv('day3ColonNegativeModulesRanked.tsv')

day3PositiveModules$tissue = 'both'
day3NegativeModules$tissue = 'both'
day3PositiveModules$selection = 'positive'
day3NegativeModules$selection = 'negative'

day3DjPositiveModules$tissue = 'dj'
day3DjNegativeModules$tissue = 'dj'
day3DjPositiveModules$selection = 'positive'
day3DjNegativeModules$selection = 'negative'

day3ColonPositiveModules$tissue = 'colon'
day3ColonNegativeModules$tissue = 'colon'
day3ColonPositiveModules$selection = 'positive'
day3ColonNegativeModules$selection = 'negative'

day3PositiveModules=day3PositiveModules %>%
  filter(!kofamAccession%in% day3DjPositiveModules$kofamAccession)%>%
  filter(!kofamAccession%in% day3ColonPositiveModules$kofamAccession)

day3NegativeModules=day3NegativeModules %>%
  filter(!kofamAccession%in% day3DjNegativeModules$kofamAccession)%>%
  filter(!kofamAccession%in% day3ColonNegativeModules$kofamAccession)

rbind(day3PositiveModules,
      day3NegativeModules,
      day3DjPositiveModules,
      day3DjNegativeModules,
      day3ColonPositiveModules,
      day3ColonNegativeModules)%>%
  na.omit()%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 3 module level 1 significant')

rbind(day3PositiveModules,
      day3NegativeModules,
      day3DjPositiveModules,
      day3DjNegativeModules,
      day3ColonPositiveModules,
      day3ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness)>2)%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 3 module level 1 strong')

rbind(day3PositiveModules,
      day3NegativeModules,
      day3DjPositiveModules,
      day3DjNegativeModules,
      day3ColonPositiveModules,
      day3ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 3 metabolism level 2 significant')

rbind(day3PositiveModules,
      day3NegativeModules,
      day3DjPositiveModules,
      day3DjNegativeModules,
      day3ColonPositiveModules,
      day3ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness)>2)%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 3 level 2 metabolism strong')

rbind(day3PositiveModules,
      day3NegativeModules,
      day3DjPositiveModules,
      day3DjNegativeModules,
      day3ColonPositiveModules,
      day3ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 3 brite level 2 significant')

rbind(day3PositiveModules,
      day3NegativeModules,
      day3DjPositiveModules,
      day3DjNegativeModules,
      day3ColonPositiveModules,
      day3ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness) > 2)%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 3 brite level 2 strong')

day7PositiveModules = read_tsv('positiveDay7ModulesRanked.tsv')
day7NegativeModules = read_tsv('negativeDay7ModulesRanked.tsv')
day7DjPositiveModules = read_tsv('day7DjPositiveModulesRanked.tsv')
day7DjNegativeModules = read_tsv('day7DjNegativeModulesRanked.tsv')
day7ColonPositiveModules = read_tsv('day7ColonPositiveModulesRanked.tsv')
day7ColonNegativeModules = read_tsv('day7ColonNegativeModulesRanked.tsv')

day7PositiveModules$tissue = 'both'
day7NegativeModules$tissue = 'both'
day7PositiveModules$selection = 'positive'
day7NegativeModules$selection = 'negative'

day7DjPositiveModules$tissue = 'dj'
day7DjNegativeModules$tissue = 'dj'
day7DjPositiveModules$selection = 'positive'
day7DjNegativeModules$selection = 'negative'

day7ColonPositiveModules$tissue = 'colon'
day7ColonNegativeModules$tissue = 'colon'
day7ColonPositiveModules$selection = 'positive'
day7ColonNegativeModules$selection = 'negative'

day7PositiveModules=day7PositiveModules %>%
  filter(!kofamAccession%in% day7DjPositiveModules$kofamAccession)%>%
  filter(!kofamAccession%in% day7ColonPositiveModules$kofamAccession)

day7NegativeModules=day7NegativeModules %>%
  filter(!kofamAccession%in% day7DjNegativeModules$kofamAccession)%>%
  filter(!kofamAccession%in% day7ColonNegativeModules$kofamAccession)

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 significant module level 1')

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness) > 2)%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 strong module level 1')

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 significant metabolism level 2')

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  filter(abs(meanFitness)> 2 )%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 strong metabolism level 2')

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 significant brite level 2')

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness) > 2)%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 significant brite level 2')

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09130 Environmental Information Processing')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 significant environmental information processing level 2')

rbind(day7PositiveModules,
      day7NegativeModules,
      day7DjPositiveModules,
      day7DjNegativeModules,
      day7ColonPositiveModules,
      day7ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness)> 2)%>%
  filter(moduleLevel1 == '09130 Environmental Information Processing')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 7 strong environmental information processing level 2')

day14PositiveModules = read_tsv('positiveDay14ModulesRanked.tsv')
day14NegativeModules = read_tsv('negativeDay14ModulesRanked.tsv')
day14DjPositiveModules = read_tsv('day14DjPositiveModulesRanked.tsv')
day14DjNegativeModules = read_tsv('day14DjNegativeModulesRanked.tsv')
day14ColonPositiveModules = read_tsv('day14ColonPositiveModulesRanked.tsv')
day14ColonNegativeModules = read_tsv('day14ColonNegativeModulesRanked.tsv')

day14PositiveModules$tissue = 'both'
day14NegativeModules$tissue = 'both'
day14PositiveModules$selection = 'positive'
day14NegativeModules$selection = 'negative'

day14DjPositiveModules$tissue = 'dj'
day14DjNegativeModules$tissue = 'dj'
day14DjPositiveModules$selection = 'positive'
day14DjNegativeModules$selection = 'negative'

day14ColonPositiveModules$tissue = 'colon'
day14ColonNegativeModules$tissue = 'colon'
day14ColonPositiveModules$selection = 'positive'
day14ColonNegativeModules$selection = 'negative'

day14PositiveModules=day14PositiveModules %>%
  filter(!kofamAccession%in% day14DjPositiveModules$kofamAccession)%>%
  filter(!kofamAccession%in% day14ColonPositiveModules$kofamAccession)

day14NegativeModules=day14NegativeModules %>%
  filter(!kofamAccession%in% day14DjNegativeModules$kofamAccession)%>%
  filter(!kofamAccession%in% day14ColonNegativeModules$kofamAccession)

rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 significant level 1')


rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  filter(abs(meanFitness)> 2)%>%
  ggplot(aes(x = selection,
             fill = moduleLevel1))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 strong level 1')

rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 significant metabolism level 2')

rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09100 Metabolism')%>%
  filter(abs(meanFitness)>2)%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 strong metabolism level 2')

rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 significant brite level 2')

rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09180 Brite Hierarchies')%>%
  filter(abs(meanFitness) > 2)%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 strong brite level 2')


?rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09130 Environmental Information Processing')%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 significant environmental information processing level 2')

rbind(day14PositiveModules,
      day14NegativeModules,
      day14DjPositiveModules,
      day14DjNegativeModules,
      day14ColonPositiveModules,
      day14ColonNegativeModules)%>%
  na.omit()%>%
  filter(moduleLevel1 == '09130 Environmental Information Processing')%>%
  filter(abs(meanFitness) > 2)%>%
  ggplot(aes(x = selection,
             fill = moduleLevel2))+
  geom_bar()+
  facet_wrap(~tissue)+
  labs(title = 'Day 14 strong environmental information processing level 2')