library (tidyverse)
day1PositiveModules = read_tsv('positiveDay1ModulesRanked.tsv')
day1NegativeModules = read_tsv('positiveDay1ModulesRanked.tsv')
day1DjPositiveModules = read_tsv('day1DjPositiveModulesRanked.tsv')
day1DjNegativeModules = read_tsv('day1DjNegativeModulesRanked.tsv')
day1ColonPositiveModules = read_tsv('dayColonPositiveModulesRanked.tsv')
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

rbind(day1PositiveModules,
      day1NegativeModules,
      day1DjPositiveModules,
      day1DjNegativeModules,
      day1ColonPositiveModules,
      day1ColonNegativeModules)%>%
  na.omit()%>%
  ggplot(aes(x = selection,
             y = n,
             fill = module))+
  geom_col()+
  facet_wrap(~tissue)