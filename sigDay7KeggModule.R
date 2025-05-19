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

day7SigNegative <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day7SigNegative.tsv')

negativeDay7Ranked <- day7SigNegative %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

negativeDay7Ranked <- annotate_modules(negativeDay7Ranked)

day7SigPositive <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/day7SigPositive.tsv')

positiveDay7Ranked <- day7SigPositive %>%
  group_by(locusId) %>%
  summarise(
    count = n(),
    meanFitness = mean(fitnessScore, na.rm = TRUE)
  ) %>%
  arrange(desc(count), desc(abs(meanFitness))) %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

positiveDay7Ranked <- annotate_modules(positiveDay7Ranked)


# Repeat for each dataset
day7SigPositiveDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjPositiveDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigPositiveDj <- annotate_modules(day7SigPositiveDj)


day7SigNegativeDj <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueDjNegativeDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigNegativeDj <- annotate_modules(day7SigNegativeDj)


day7SigPositiveColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonPositiveDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigPositiveColon <- annotate_modules(day7SigPositiveColon)


day7SigNegativeColon <- read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/B. breve barseq/barseq/uniqueColonNegativeDay7.tsv') %>%
  left_join(keggAnno, by = "locusId") %>%
  rownames_to_column('rank')

day7SigNegativeColon <- annotate_modules(day7SigNegativeColon)


# Example summary
day7SigNegativeColon %>%
  write_tsv('day7ColonNegativeModulesRanked.tsv')

day7SigPositiveColon %>%
  write_tsv('day7ColonPositiveModulesRanked.tsv')

day7SigNegativeDj%>%
  write_tsv('day7DjNegativeModulesRanked.tsv')

day7SigPositiveDj%>%
  write_tsv('day7DjPositiveModulesRanked.tsv')

negativeDay7Ranked%>%
  write_tsv('negativeDay7ModulesRanked.tsv')

positiveDay7Ranked%>%

  write_tsv('positiveDay7ModulesRanked.tsv')