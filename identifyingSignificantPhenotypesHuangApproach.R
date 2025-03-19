library(tidyverse)

genSigPhenotypeTab = function(longFitTab, dayOfInterest){
  sigPhenotypeTab = longFitTab %>%
    filter(day == dayOfInterest)%>%
    filter(abs(tscore)> 4)%>%
    filter(abs(fitnessScore) > .5)
  return (sigPhenotypeTab)
}

genSigNegativePhenotab = function(sigPhenoTab){

  negativePhenotypes = sigPhenoTab %>%
    filter(fitnessScore < 0)

  return(negativePhenotypes)
}

genSigPositivePhenotab = function(sigPhenoTab){

  positivePhenotypes = sigPhenoTab %>%
    filter(fitnessScore > 0)

  return(positivePhenotypes)
}

genTissuePheno = function (tissueOfInterest, sigTab){
  sigTab%>%
    filter(tissue == tissueOfInterest)%>%
    return()
}

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

day1SigPhenotypes = genSigPhenotypeTab(mergedTandFitnessScoresLong, dayOfInterest = 'day1')
day1SigNegative= genSigNegativePhenotab(day1SigPhenotypes)
write_tsv(day1SigNegative, file = 'day1SigNegative.tsv')

day1SigPositive = genSigPositivePhenotab(day1SigPhenotypes)
write_tsv(day1SigPositive, file = 'day1SigPositive.tsv')


day1ColonSigNegative = genTissuePheno('colon', day1SigNegative)
day1ColonSigPositive = genTissuePheno('colon', day1SigPositive)
day1DjSigNegative = genTissuePheno('dj', day1SigNegative)
day1DjSigPositive = genTissuePheno('dj', day1SigPositive)

'40 significant positive phenotypes in the colon on the first day'
day1ColonSigPositive%>%
  select(locusId)%>%
  distinct()%>%
  nrow()


'194 signficant negative phenotypes in the colon on the first day'
day1ColonSigNegative%>%
  select(locusId)%>%
  distinct()%>%
  nrow()
'247 significant postive phenotypes in the dj on the first day'
day1DjSigPositive%>%
  select(locusId)%>%
  distinct()%>%
  nrow()
'161 significant negative phenotypes in the dj on the first day'
day1DjSigNegative%>%
  select(locusId)%>%
  distinct()%>%
  nrow()

intersect(day1ColonSigNegative$locusId,
          day1ColonSigPositive$locusId)
'One shared mutants dj with signficant positive and negative score in the SI day1DjSigPositive: BBR_RS20320.
This is actually the PTS identified by the huang lab. Might be a reason to pool the samples from each day'
intersect(day1DjSigNegative$locusId,
          day1DjSigPositive$locusId)
'of 247 significant positive phenotypes in the SI, 219 are unique'
day1DjSigPositive %>%
  filter(!locusId %in% day1ColonSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()%>%
  nrow()

'of 194 significant positive phenotypes in the colon, 12 are unique'
 day1ColonSigPositive%>%
  filter(!locusId %in% day1DjSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()%>%
  nrow()

uniqueColonicPositiveLoci=day1ColonSigPositive%>%
  filter(!locusId %in% day1DjSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

uniqueColonPositiveDay1=day1ColonSigPositive %>%
  filter(locusId %in% uniqueColonicPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(meanFitness))

write_tsv(uniqueColonPositiveDay1, 'uniqueColonPositiveDay1.tsv')

uniqueDjPositiveLoci=day1DjSigPositive%>%
  filter(!locusId %in% day1ColonSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

uniqueDjPositiveDay1=day1DjSigPositive %>%
  filter(locusId %in% uniqueDjPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(meanFitness))
write_tsv(uniqueDjPositiveDay1, 'uniqueDjPositiveDay1.tsv')

sharedPositiveLociDay1=day1DjSigPositive %>%
  filter(!locusId %in% uniqueDjPositiveLoci$locusId & ! locusId %in% uniqueColonicPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(meanFitness))

'84 unique negative loci'
uniqueColonicNegativeLoci=day1ColonSigNegative%>%
  filter(!locusId %in% day1DjSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

uniqueColonNegativeDay1=day1ColonSigNegative %>%
  filter(locusId %in% uniqueColonicNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))

write_tsv(uniqueColonNegativeDay1, 'uniqueColonNegativeDay1.tsv')


'51 unique negative dj loci'
uniqueDjNegativeLoci=day1DjSigNegative%>%
  filter(!locusId %in% day1ColonSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

uniqueDjNegativeDay1=day1DjSigNegative %>%
  filter(locusId %in% uniqueDjNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))

write_tsv(uniqueDjNegativeDay1, 'uniqueDjNegativeDay1.tsv')


'110 shared negative loci on day 1 '
sharedNegativeDay1 = day1SigNegative %>%
  filter(!locusId %in% uniqueColonicNegativeLoci$locusId & !locusId %in% uniqueDjNegativeLoci$locusId)%>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))

day1DjSigNegative %>%
  filter(locusId %in% day1ColonSigPositive$locusId)
'6 loci that are positive in the SI but negative in the colon on day 1'

day1DjSigPositive %>%
  filter(locusId %in% day1ColonSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

day1DjSigPositive %>%
  filter(locusId %in% day1ColonSigNegative$locusId)%>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))


day3SigPhenotypes = genSigPhenotypeTab(mergedTandFitnessScoresLong, dayOfInterest = 'day3')

day3SigNegative= genSigNegativePhenotab(day3SigPhenotypes)
write_tsv(day3SigNegative, 'day3SigNegative.tsv')
day3SigPositive = genSigPositivePhenotab(day3SigPhenotypes)
write_tsv(day3SigPositive, 'day3SigPositive.tsv')
day3ColonSigNegative = genTissuePheno('colon', day3SigNegative)
day3ColonSigPositive = genTissuePheno('colon', day3SigPositive)
day3DjSigNegative = genTissuePheno('dj', day3SigNegative)
day3DjSigPositive = genTissuePheno('dj', day3SigPositive)

uniqueDjNegativeLoci=day3DjSigNegative%>%
  filter(!locusId %in% day3ColonSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

uniqueDjNegativeDay3=day3DjSigNegative %>%
  filter(locusId %in% uniqueDjNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))

write_tsv(uniqueDjNegativeDay3, 'uniqueDjNegativeDay3.tsv')

uniqueDjNegativeLoci=day3DjSigNegative%>%
  filter(!locusId %in% day3ColonSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

'45 unique negative dj'
uniqueDjNegativeDay3=day3DjSigNegative %>%
  filter(locusId %in% uniqueDjNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))

write_tsv(uniqueDjNegativeDay3, 'uniqueDjNegativeDay3.tsv')

uniqueColonNegativeLoci=day3ColonSigNegative%>%
  filter(!locusId %in% day3DjSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

'84 unique negative colonic pheno'
uniqueColonNegativeDay3=day3ColonSigNegative %>%
  filter(locusId %in% uniqueColonNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueColonNegativeDay3, 'uniqueColonNegativeDay3.tsv')

uniqueDjPositiveLoci=day3DjSigPositive%>%
  filter(!locusId %in% day3ColonSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

uniqueDjPositiveDay3=day3DjSigPositive %>%
  filter(locusId %in% uniqueDjPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueDjPositiveDay3, 'uniqueDjPositiveDay3.tsv')

uniqueDjPositiveLoci=day3DjSigPositive%>%
  filter(!locusId %in% day3ColonSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

'19 unique positive dj'
uniqueDjPositiveDay3=day3DjSigPositive %>%
  filter(locusId %in% uniqueDjPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))


uniqueColonPositiveLoci=day3ColonSigPositive%>%
  filter(!locusId %in% day3DjSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

'95 unique positive colon day 3'
uniqueColonPositiveDay3=day3ColonSigPositive %>%
  filter(locusId %in% uniqueColonPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))

write_tsv(uniqueColonPositiveDay3, 'uniqueColonPositiveDay3.tsv')

day7SigPhenotypes = genSigPhenotypeTab(mergedTandFitnessScoresLong, dayOfInterest = 'day7')

day7SigNegative= genSigNegativePhenotab(day7SigPhenotypes)
write_tsv(day7SigNegative, 'day7SigNegative.tsv')

day7SigPositive = genSigPositivePhenotab(day7SigPhenotypes)
write_tsv(day7SigPositive, 'day7SigPositive.tsv')

day7ColonSigNegative = genTissuePheno('colon', day7SigNegative)
day7ColonSigPositive = genTissuePheno('colon', day7SigPositive)
day7DjSigNegative = genTissuePheno('dj', day7SigNegative)
day7DjSigPositive = genTissuePheno('dj', day7SigPositive)

uniqueDjNegativeLoci=day7DjSigNegative%>%
  filter(!locusId %in% day7ColonSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

'18 unique negative phenotypes in the dh'
uniqueDjNegativeDay7=day7DjSigNegative %>%
  filter(locusId %in% uniqueDjNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueDjNegativeDay7, 'uniqueDjNegativeDay7.tsv')

uniqueDjNegativeLoci=day7DjSigNegative%>%
  filter(!locusId %in% day7ColonSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

'45 unique negative dj'
uniqueDjNegativeDay7=day7DjSigNegative %>%
  filter(locusId %in% uniqueDjNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueDjNegativeDay7, 'uniqueDjNegativeDay7.tsv')

uniqueColonNegativeLoci=day7ColonSigNegative%>%
  filter(!locusId %in% day7DjSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()
'48 unique negative colonic pheno'
uniqueColonNegativeDay7=day7ColonSigNegative %>%
  filter(locusId %in% uniqueColonNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueColonNegativeDay7, 'uniqueColonNegativeDay7.tsv')

uniqueDjPositiveLoci=day7DjSigPositive%>%
  filter(!locusId %in% day7ColonSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

uniqueDjPositiveDay7=day7DjSigPositive %>%
  filter(locusId %in% uniqueDjPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueDjPositiveDay7, 'uniqueDjPositiveDay7.tsv')


uniqueDjPositiveLoci=day7DjSigPositive%>%
  filter(!locusId %in% day7ColonSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

'10 unique positive dj'
uniqueDjPositiveDay7=day7DjSigPositive %>%
  filter(locusId %in% uniqueDjPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueColonNegativeDay7, 'uniqueColonNegativeDay7.tsv')


uniqueColonPositiveLoci=day7ColonSigPositive%>%
  filter(!locusId %in% day7DjSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

'10 unique positive colon day 3'
uniqueColonPositiveDay7=day7ColonSigPositive %>%
  filter(locusId %in% uniqueColonPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueColonPositiveDay7, 'uniqueColonPositiveDay7.tsv')
day14SigPhenotypes = genSigPhenotypeTab(mergedTandFitnessScoresLong, dayOfInterest = 'day14')

day14SigNegative= genSigNegativePhenotab(day14SigPhenotypes)
write_tsv(day14SigNegative, 'day14SigNegative.tsv')

day14SigPositive = genSigPositivePhenotab(day14SigPhenotypes)
write_tsv(day14SigPositive, 'day14SigPositive.tsv')

day14ColonSigNegative = genTissuePheno('colon', day14SigNegative)
day14ColonSigPositive = genTissuePheno('colon', day14SigPositive)
day14DjSigNegative = genTissuePheno('dj', day14SigNegative)
day14DjSigPositive = genTissuePheno('dj', day14SigPositive)

uniqueDjNegativeLoci=day14DjSigNegative%>%
  filter(!locusId %in% day14ColonSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

'66 unique negative dj'
uniqueDjNegativeDay14=day14DjSigNegative %>%
  filter(locusId %in% uniqueDjNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueDjNegativeDay14, 'uniqueDjNegativeDay14.tsv')

uniqueColonNegativeLoci=day14ColonSigNegative%>%
  filter(!locusId %in% day14DjSigNegative$locusId)%>%
  select(locusId)%>%
  distinct()

'104 unique negative colonic pheno'
uniqueColonNegativeDay14=day14ColonSigNegative %>%
  filter(locusId %in% uniqueColonNegativeLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueColonNegativeDay14, 'uniqueColonNegativeDay14.tsv')

uniqueDjPositiveLoci=day14DjSigPositive%>%
  filter(!locusId %in% day14ColonSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()
'47 unique dj positive'
uniqueDjPositiveDay14=day14DjSigPositive %>%
  filter(locusId %in% uniqueDjPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueDjPositiveDay14, 'uniqueDjPositiveDay14.tsv')

uniqueColonPositiveLoci=day14ColonSigPositive%>%
  filter(!locusId %in% day14DjSigPositive$locusId)%>%
  select(locusId)%>%
  distinct()

'17 unique positive colon day 7'
uniqueColonPositiveDay14=day14ColonSigPositive %>%
  filter(locusId %in% uniqueColonPositiveLoci$locusId) %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness)))
write_tsv(uniqueColonPositiveDay14, 'uniqueColonPositiveDay14.tsv')

