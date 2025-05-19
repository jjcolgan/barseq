library(tidyverse)
library(KEGGREST)
library(progressr)

# Set up the progress handler
progressr::handlers("txtprogressbar")

# KEGG module fetcher
get_module <- function(accession) {
  keggout <- keggGet(accession)
  if (!is.null(keggout[[1]]$BRITE)) {
    brite <- keggout[[1]]$BRITE
    level1 <- if (length(brite) >= 2) as.character(brite[2]) else NA_character_
    level2 <- if (length(brite) >= 3) as.character(brite[3]) else NA_character_
    return(c(level1, level2))
  } else {
    return(c(NA_character_, NA_character_))
  }
}

annotate_modules <- function(df) {
  df_no_na <- df %>% filter(!is.na(kofamAccession)) %>%
    mutate(moduleLevel1 = NA_character_, moduleLevel2 = NA_character_)

  with_progress({
    p <- progressor(steps = nrow(df_no_na))
    for (i in seq_len(nrow(df_no_na))) {
      p()
      keggOut <- get_module(df_no_na$kofamAccession[i])
      df_no_na$moduleLevel1[i] <- keggOut[1]
      df_no_na$moduleLevel2[i] <- keggOut[2]
    }
  })

  df_annotated <- df %>%
    left_join(df_no_na %>% select(kofamAccession, moduleLevel1, moduleLevel2),
              by = "kofamAccession")

  return(df_annotated)
}
# Load annotations and datasets
keggAnno <- read_tsv('genesWithAnvioAnnotations.tsv') %>%
  select(locusId, kofamAccession) %>%
  distinct()

day3SigNegative <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day3SigNegative.tsv')

negativeDay3Ranked <- day3SigNegative %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

negativeDay3Ranked <- annotate_modules(negativeDay3Ranked)

day3SigPositive <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day3SigPositive.tsv')

positiveDay3Ranked <- day3SigPositive %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

positiveDay3Ranked <- annotate_modules(positiveDay3Ranked)


# Repeat for each dataset
day3SigPositiveDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjPositiveDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigPositiveDj <- annotate_modules(day3SigPositiveDj)


day3SigNegativeDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjNegativeDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigNegativeDj <- annotate_modules(day3SigNegativeDj)


day3SigPositiveColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonPositiveDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigPositiveColon <- annotate_modules(day3SigPositiveColon)


day3SigNegativeColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonNegativeDay3.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day3SigNegativeColon <- annotate_modules(day3SigNegativeColon)


# Example summary
day3SigNegativeColon %>%
  write_tsv('day3ColonNegativeModulesRanked.tsv')

day3SigPositiveColon %>%
  write_tsv('day3ColonPositiveModulesRanked.tsv')

day3SigNegativeDj%>%
  write_tsv('day3DjNegativeModulesRanked.tsv')

day3SigPositiveDj%>%
  write_tsv('day3DjPositiveModulesRanked.tsv')

negativeDay3Ranked%>%
  write_tsv('negativeDay3ModulesRanked.tsv')

positiveDay3Ranked%>%
  write_tsv('positiveDay3ModulesRanked.tsv')