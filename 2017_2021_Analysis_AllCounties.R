# 2017-2021 Analysis: Hennepin, Franklin, and King County:

library(tidyverse)
library(tidycensus)
library(dplyr)
library(cluster)
library(factoextra)
library(stringr)
library(sf)
library(ggplot2)
library(scales)

# Census Key
census_api_key("3c9e13d499f9ed0c1d8f1c848dd6919fc443afad")

# load in crosswalked data
crosswalk <- read_csv("~/Documents/0Class/FINAL PROJECT/crosswalked.csv")

# Separate dataframes for each county:
crosswalk <- crosswalk %>%
  mutate(
    non_som = total - som_total,
    percent_som = som_total / total
  )

hennepin_data <- crosswalk %>%
  # Filter for GEOIDs that start with 27053 (Hennepin County, MN)
  filter(substr(GEOID_2010, 1, 5) == "27053")

frank_data <- crosswalk %>%
  # Filter for GEOIDs that start with 39049 (Franklin County, OH)
  filter(substr(GEOID_2010, 1, 5) == "39049")

king_data <- crosswalk %>%
  # Filter for GEOIDs that start with 53033 (King County, WA)
  filter(substr(GEOID_2010, 1, 5) == "53033")

## Calculating indices for each county + tract
# Hennepin
# dissimilarity index (segregation)
s_2021_hen <- sum(hennepin_data$som_total, na.rm = TRUE)
t_2021_hen <- sum(hennepin_data$non_som, na.rm = TRUE)

hennepin_data$dissimilarity <- abs((hennepin_data$som_total/s_2021_hen) - 
                                (hennepin_data$non_som/t_2021_hen))

hennepin_dissim <- 0.5 * sum(hennepin_data$dissimilarity, na.rm = TRUE)
hennepin_dissim

# isolation index ()
hennepin_data$isolation <- (hennepin_data$som_total/s_2021_hen) * (hennepin_data$som_total/hennepin_data$total)
hennepin_isolation <- sum(hennepin_data$isolation, na.rm = TRUE)
hennepin_isolation

## Proportional weighting of Somali languages for each tract (Proxy estimation)
### Hennepin 2020
hennepin_data <- hennepin_data %>%
  mutate(
    somali_speakersEST = ifelse(som_total>0, other_lang_total*percent_som, 0), # estimate of total somali speakers
    ESTsom_very_well = ifelse(som_total>0, other_very_well*percent_som, 0), # estimate of somali speakers who speak english very well
    ESTsom_less = ifelse(som_total>0, other_less*percent_som, 0), # estimate of somali speakers who speak english less than very well
    pct_somali_EST = ifelse(som_total>0, somali_speakersEST/total, 0), # percent estimate of somali speakers
    EST_pct_very_well = ifelse(somali_speakersEST>0, ESTsom_very_well/somali_speakersEST, 0), # percent estimate of somali speakers who speak english very well
    EST_pct_less = ifelse(somali_speakersEST>0, ESTsom_less/somali_speakersEST, 0) # percent estimate of somali speakers who speak english less than very well
  )

## avg estimate of very well vs less than very well
total_speakers_hen <- sum(hennepin_data$somali_speakersEST, na.rm = TRUE)

total_very_well_hen <- sum(hennepin_data$ESTsom_very_well, na.rm = TRUE)
total_less_hen <- sum(hennepin_data$ESTsom_less, na.rm = TRUE)

# s_2021_hen = total somali population in hennepin county
hennepin_very_well_21 <- (total_very_well_hen/s_2021_hen) * 100
hennepin_less_21 <- (total_less_hen/s_2021_hen) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

## Entropy score formula (multigroup): E = -∑(pi × ln(pi)) 
calc_entropy <- function(props) {
  # Remove NA or 0 values since log(0) is undefined
  props <- props[!is.na(props) & props > 0]
  if (length(props) == 0) return(0)
  -sum(props * log(props))
}

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
hennepin_data <- hennepin_data %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

hennepin_data$local_entropy[is.nan(hennepin_data$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
hen_totals <- hennepin_data %>%
  summarize(
    total_speakers_hen,
    total_very_well_hen,
    total_less_hen,
    # proportions:
    p_very_well = total_very_well_hen / total_speakers_hen,
    p_less = total_less_hen / total_speakers_hen
  )

metro_entropy_hen <- with(hen_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_hen

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
hennepin_data <- hennepin_data %>%
  mutate(metro_entropy_hen = metro_entropy_hen)

# weight by somali pop

hennepin_data <- hennepin_data %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(hennepin_data$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_hen - local_entropy) / metro_entropy_hen,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
hen_entropy_index_org <- sum(hennepin_data$entropy_weight, na.rm = TRUE)
hen_entropy_index_org

## checking negative values in entropy difference (theory contradiction)
# 
negative_tracts <- hennepin_data %>%
  filter(entrop_difference < 0) %>%
  select(GEOID_2010, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy, entrop_difference)

# index -- summation: (adjusted)
hennepin_data <- hennepin_data %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

hen_entropy_index_adj <- sum(hennepin_data$entropy_weight_adj, na.rm = TRUE)
hen_entropy_index_adj

## 1.58% change between original and adjusted

#############################
## Calculating indices for each county + tract
# Franklin
# dissimilarity index (segregation)
s_2021_frank <- sum(frank_data$som_total, na.rm = TRUE)
t_2021_frank <- sum(frank_data$non_som, na.rm = TRUE)

frank_data$dissimilarity <- abs((frank_data$som_total/s_2021_frank) - 
                                     (frank_data$non_som/t_2021_frank))

franklin_dissim <- 0.5 * sum(frank_data$dissimilarity, na.rm = TRUE)
franklin_dissim

# isolation index ()
frank_data$isolation <- (frank_data$som_total/s_2021_frank) * (frank_data$som_total/frank_data$total)
franklin_isolation <- sum(frank_data$isolation, na.rm = TRUE)
franklin_isolation

## Proportional weighting of Somali languages for each tract (Proxy estimation)
### Franklin 2020
frank_data <- frank_data %>%
  mutate(
    somali_speakersEST = ifelse(som_total>0, other_lang_total*percent_som, 0), # estimate of total somali speakers
    ESTsom_very_well = ifelse(som_total>0, other_very_well*percent_som, 0), # estimate of somali speakers who speak english very well
    ESTsom_less = ifelse(som_total>0, other_less*percent_som, 0), # estimate of somali speakers who speak english less than very well
    pct_somali_EST = ifelse(som_total>0, somali_speakersEST/total, 0), # percent estimate of somali speakers
    EST_pct_very_well = ifelse(somali_speakersEST>0, ESTsom_very_well/somali_speakersEST, 0), # percent estimate of somali speakers who speak english very well
    EST_pct_less = ifelse(somali_speakersEST>0, ESTsom_less/somali_speakersEST, 0) # percent estimate of somali speakers who speak english less than very well
  )

## avg estimate of very well vs less than very well
total_speakers_frank <- sum(frank_data$somali_speakersEST, na.rm = TRUE)

total_very_well_frank <- sum(frank_data$ESTsom_very_well, na.rm = TRUE)
total_less_frank <- sum(frank_data$ESTsom_less, na.rm = TRUE)

# s_2021_frank = total somali population in franklin county
franklin_very_well_21 <- (total_speakers_frank/s_2021_frank) * 100
franklin_less_21 <- (total_less_frank/s_2021_frank) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

## Entropy score formula (multigroup): E = -∑(pi × ln(pi)) 
calc_entropy <- function(props) {
  # Remove NA or 0 values since log(0) is undefined
  props <- props[!is.na(props) & props > 0]
  if (length(props) == 0) return(0)
  -sum(props * log(props))
}

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
frank_data <- frank_data %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

frank_data$local_entropy[is.nan(frank_data$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
frank_totals <- frank_data %>%
  summarize(
    total_speakers_frank,
    total_very_well_frank,
    total_less_frank,
    # proportions:
    p_very_well = total_very_well_frank / total_speakers_frank,
    p_less = total_less_frank / total_speakers_frank
  )

metro_entropy_frank <- with(frank_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_frank

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
frank_data <- frank_data %>%
  mutate(metro_entropy_frank = metro_entropy_frank)

# weight by somali pop
frank_data <- frank_data %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(frank_data$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_frank - local_entropy) / metro_entropy_frank,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
frank_entropy_index_org <- sum(frank_data$entropy_weight, na.rm = TRUE)
frank_entropy_index_org

# index -- summation: (adjusted)
frank_data <- frank_data %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

frank_entropy_index_adj <- sum(frank_data$entropy_weight_adj, na.rm = TRUE)
frank_entropy_index_adj

# See how much the adjustment changes the overall index
f.dif <- frank_entropy_index_adj - frank_entropy_index_org
f.p <- (f.dif / frank_entropy_index_org) * 100
f.p

## 4.36% change between original and adjusted





#############################
## Calculating indices for each county + tract
# King
# dissimilarity index (segregation)
s_2021_king <- sum(king_data$som_total, na.rm = TRUE)
t_2021_king <- sum(king_data$non_som, na.rm = TRUE)

king_data$dissimilarity <- abs((king_data$som_total/s_2021_king) - 
                                  (king_data$non_som/t_2021_king))

king_dissim <- 0.5 * sum(king_data$dissimilarity, na.rm = TRUE)
king_dissim

# isolation index ()
king_data$isolation <- (king_data$som_total/s_2021_king) * (king_data$som_total/king_data$total)
king_isolation <- sum(king_data$isolation, na.rm = TRUE)
king_isolation

## Proportional weighting of Somali languages for each tract (Proxy estimation)
### King 2021
king_data <- king_data %>%
  mutate(
    somali_speakersEST = ifelse(som_total>0, other_lang_total*percent_som, 0), # estimate of total somali speakers
    ESTsom_very_well = ifelse(som_total>0, other_very_well*percent_som, 0), # estimate of somali speakers who speak english very well
    ESTsom_less = ifelse(som_total>0, other_less*percent_som, 0), # estimate of somali speakers who speak english less than very well
    pct_somali_EST = ifelse(som_total>0, somali_speakersEST/total, 0), # percent estimate of somali speakers
    EST_pct_very_well = ifelse(somali_speakersEST>0, ESTsom_very_well/somali_speakersEST, 0), # percent estimate of somali speakers who speak english very well
    EST_pct_less = ifelse(somali_speakersEST>0, ESTsom_less/somali_speakersEST, 0) # percent estimate of somali speakers who speak english less than very well
  )

## avg estimate of very well vs less than very well
total_speakers_king <- sum(king_data$somali_speakersEST, na.rm = TRUE)

total_very_well_king <- sum(king_data$ESTsom_very_well, na.rm = TRUE)
total_less_king <- sum(king_data$ESTsom_less, na.rm = TRUE)

# s_2021_king = total somali population in franklin county
king_very_well_21 <- (total_speakers_king/s_2021_king) * 100
king_less_21 <- (total_less_king/s_2021_king) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

## Entropy score formula (multigroup): E = -∑(pi × ln(pi)) 
calc_entropy <- function(props) {
  # Remove NA or 0 values since log(0) is undefined
  props <- props[!is.na(props) & props > 0]
  if (length(props) == 0) return(0)
  -sum(props * log(props))
}

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
king_data <- king_data %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

king_data$local_entropy[is.nan(king_data$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
king_totals <- king_data %>%
  summarize(
    total_speakers_king,
    total_very_well_king,
    total_less_king,
    # proportions:
    p_very_well = total_very_well_king / total_speakers_king,
    p_less = total_less_king / total_speakers_king
  )

metro_entropy_king <- with(king_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_king

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
king_data <- king_data %>%
  mutate(metro_entropy_king = metro_entropy_king)

# weight by somali pop
king_data <- king_data %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(king_data$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_king - local_entropy) / metro_entropy_king,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
king_entropy_index_org <- sum(king_data$entropy_weight, na.rm = TRUE)
king_entropy_index_org

# index -- summation: (adjusted)
king_data <- king_data %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

problem_tracts <- king_data %>%
  filter(entrop_difference < 0) %>%
  select(GEOID_2010, somali_speakersEST, EST_pct_very_well, EST_pct_less, 
         local_entropy, metro_entropy_king, entrop_difference)
problem_tracts

king_entropy_index_adj <- sum(king_data$entropy_weight_adj, na.rm = TRUE)
king_entropy_index_adj

# See how much the adjustment changes the overall index
k.dif <- king_entropy_index_adj - king_entropy_index_org
k.p <- (k.dif / king_entropy_index_org) * 100
k.p

## 0.819% change between original and adjusted

