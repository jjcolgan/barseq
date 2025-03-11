library(tidyverse)
library(ggpubr)

fitness = read_tsv('barseqAdjustedParams/fit_logratios.tab')
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

metadata$day = factor(metadata$day, levels = c('t0', 'day1','day3', 'day7', 'day14'))
annotations = read_tsv('genesWithAnvioAnnotations.tsv')
fitnessFullLong = fitness%>%
  pivot_longer(cols = c(4:45), names_to = 'sample', values_to = 'fitnessScore')%>%
  left_join(metadata, by = 'sample')

day1MeanTissueFitnessScores=fitnessFullLong %>%
  filter(day == 'day1')%>%
  group_by(locusId, tissue)%>%
  summarise(meanTissueFitness = mean(fitnessScore))%>%
  pivot_wider(id_cols = 'locusId',
              names_from = 'tissue',
              values_from = 'meanTissueFitness')
makeMeanTissueFitnessScatterPlot = function(longFitTab, dayOfInterest){
  meanDayWide=longFitTab %>%
    filter(day == dayOfInterest)%>%
    group_by(locusId, tissue)%>%
    summarise(meanTissueFitness = mean(fitnessScore))%>%
    pivot_wider(id_cols = 'locusId',
                names_from = 'tissue',
                values_from = 'meanTissueFitness')

  meanDayWide$phenotype = 'Deleterious in both populations'
  meanDayWide$phenotype[meanDayWide$dj > 0 & meanDayWide$colon < 0 ] = 'Weak positive Dj weak negative colon'
  meanDayWide$phenotype[meanDayWide$dj >= 2 & meanDayWide$colon < 0 ] = 'Strong positive Dj weak negative colon'
  meanDayWide$phenotype[meanDayWide$dj >=2 & meanDayWide$colon <= -2 ] = 'Strong positive Dj strong negative colon'
  meanDayWide$phenotype[meanDayWide$dj < 0 & meanDayWide$colon > 0 ] = 'Weak negative Dj weak positive colon'
  meanDayWide$phenotype[meanDayWide$dj <= -2 & meanDayWide$colon > 0 ] = 'Strong negative Dj weak positive colon'
  meanDayWide$phenotype[meanDayWide$dj <= -2 & meanDayWide$colon >= 2 ] = 'Strong negative Dj strong positive colon'
  meanDayWide$phenotype[meanDayWide$dj > 0 & meanDayWide$colon >0  ] = 'Weak positive both'
  meanDayWide$phenotype[meanDayWide$dj >= 2 & meanDayWide$colon >=2  ] = 'Strong positive both'
  meanDayWide$phenotype[meanDayWide$dj < 0 & meanDayWide$colon <0  ] = 'Weak negative both'
  meanDayWide$phenotype[meanDayWide$dj <= -2 & meanDayWide$colon <-2  ] = 'Strong negative both'
  p=meanDayWide %>%
    ggplot(aes(x = dj,
               col = phenotype,
               y = colon))+
    geom_point()+
    labs(title = dayOfInterest)+
    xlim(-10,10)+
    ylim(-10,10)+
    theme_bw()
  plot(p)
}

makeMeanTissueFitnessScatterPlot(fitnessFullLong, 'day1')
makeMeanTissueFitnessScatterPlot(fitnessFullLong, 'day3')
makeMeanTissueFitnessScatterPlot(fitnessFullLong, 'day7')
makeMeanTissueFitnessScatterPlot(fitnessFullLong, 'day14')
