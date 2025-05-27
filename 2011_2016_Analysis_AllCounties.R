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

# Somali population variables (using ancestry for population)
variables <- c(
  total = "B04006_001E", # Total population
  som_total = "B04006_082E", # Somali alone
  # B-Series (county generalizations)
  afro_lang = "B16001_111E", # Amharic, Somali, or other Afro-Asiatic languages
  afro_very_well  = "B16001_112E", # Speakers who speak English "very well"
  afro_less  = "B16001_113E", # Speakers who speak English less than well
  # C-Series (language proxy)
  other_lang_total = "C16001_036E", # Other & unspecified languages
  other_very_well = "C16001_037E", # Speakers who speak English "very well"
  other_less = "C16001_038E" # Speakers who speak English less than well
)

# Hennepin (2011-2016)
hen <- get_acs(
  geography = "tract",
  variables = variables,
  state = "MN",
  county = "Hennepin",
  year = 2016,
  output = "wide",
  geometry = TRUE
)

# Franklin (2011-2016)
frank <- get_acs(
  geography = "tract",
  variables = variables,
  state = "OH",
  county = "Franklin",
  year = 2016,
  output = "wide",
  geometry = TRUE
)

# King (2011-2016)
king <- get_acs(
  geography = "tract",
  variables = variables,
  state = "WA",
  county = "King",
  year = 2016,
  output = "wide",
  geometry = TRUE
)

# removing moe columns
hen <- hen %>% select(-c(B04006_001M, B04006_082M, B16001_111M, B16001_112M, B16001_113M,
                                   C16001_036M, C16001_037M, C16001_038M))

frank <- frank %>% select(-c(B04006_001M, B04006_082M, B16001_111M, B16001_112M, B16001_113M,
                         C16001_036M, C16001_037M, C16001_038M))

king <- king %>% select(-c(B04006_001M, B04006_082M, B16001_111M, B16001_112M, B16001_113M,
                         C16001_036M, C16001_037M, C16001_038M))

## Calculating indices for each county + tract
# Hennepin
hen <- hen %>%
  mutate(
    non_som = total - som_total,
    percent_som = som_total / total
  )

# dissimilarity index (segregation)
s_2016_hen <- sum(hen$som_total, na.rm = TRUE)
t_2016_hen <- sum(hen$non_som, na.rm = TRUE)

hen$dissimilarity <- abs((hen$som_total/s_2016_hen) - 
                                     (hen$non_som/t_2016_hen))

hennepin_dissim <- 0.5 * sum(hen$dissimilarity, na.rm = TRUE)
hennepin_dissim

# isolation index ()
hen$isolation <- (hen$som_total/s_2016_hen) * (hen$som_total/hen$total)
hennepin_isolation <- sum(hen$isolation, na.rm = TRUE)
hennepin_isolation

## Proportional weighting of Somali languages for each tract (Proxy estimation)
### Hennepin 2020
hen <- hen %>%
  mutate(
    somali_speakersEST = ifelse(som_total>0, other_lang_total*percent_som, 0), # estimate of total somali speakers
    ESTsom_very_well = ifelse(som_total>0, other_very_well*percent_som, 0), # estimate of somali speakers who speak english very well
    ESTsom_less = ifelse(som_total>0, other_less*percent_som, 0), # estimate of somali speakers who speak english less than very well
    pct_somali_EST = ifelse(som_total>0, somali_speakersEST/total, 0), # percent estimate of somali speakers
    EST_pct_very_well = ifelse(somali_speakersEST>0, ESTsom_very_well/somali_speakersEST, 0), # percent estimate of somali speakers who speak english very well
    EST_pct_less = ifelse(somali_speakersEST>0, ESTsom_less/somali_speakersEST, 0) # percent estimate of somali speakers who speak english less than very well
  )

## avg estimate of very well vs less than very well
total_speakers_hen <- sum(hen$somali_speakersEST, na.rm = TRUE)

total_very_well_hen <- sum(hen$ESTsom_very_well, na.rm = TRUE)
total_less_hen <- sum(hen$ESTsom_less, na.rm = TRUE)

# s_2016_hen = total somali population in hennepin county
hennepin_very_well_16 <- (total_very_well_hen/s_2016_hen) * 100
hennepin_less_16 <- (total_less_hen/s_2016_hen) * 100

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
hen <- hen %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

hen$local_entropy[is.nan(hen$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
hen_totals <- hen %>%
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
hen <- hen %>%
  mutate(metro_entropy_hen = metro_entropy_hen)

# weight by somali pop

hen <- hen %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(hen$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_hen - local_entropy) / metro_entropy_hen,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
hen_entropy_index_org <- sum(hen$entropy_weight, na.rm = TRUE)
hen_entropy_index_org

## checking negative values in entropy difference (theory contradiction)
# 
h.negative_tracts <- hen %>%
  filter(entrop_difference < 0) %>%
  select(GEOID, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy_hen, entrop_difference)

# index -- summation: (adjusted)
hen <- hen %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

hen_entropy_index_adj <- sum(hen$entropy_weight_adj, na.rm = TRUE)
hen_entropy_index_adj

## percent difference between original and adjusted:
h.dif <- hen_entropy_index_adj - hen_entropy_index_org
h.p <- (h.dif / hen_entropy_index_org) * 100
h.p

## 0.784% change between original and adjusted



#######################################
# franklin
frank <- frank %>%
  mutate(
    non_som = total - som_total,
    percent_som = som_total / total
  )

# dissimilarity index (segregation)
s_2016_frank <- sum(frank$som_total, na.rm = TRUE)
t_2016_frank <- sum(frank$non_som, na.rm = TRUE)

frank$dissimilarity <- abs((frank$som_total/s_2016_frank) - 
                           (frank$non_som/t_2016_frank))

frank_dissim <- 0.5 * sum(frank$dissimilarity, na.rm = TRUE)
frank_dissim

# isolation index ()
frank$isolation <- (frank$som_total/s_2016_frank) * (frank$som_total/frank$total)
franklin_isolation <- sum(frank$isolation, na.rm = TRUE)
franklin_isolation

## Proportional weighting of Somali languages for each tract (Proxy estimation)
### franklin 2016
frank <- frank %>%
  mutate(
    somali_speakersEST = ifelse(som_total>0, other_lang_total*percent_som, 0), # estimate of total somali speakers
    ESTsom_very_well = ifelse(som_total>0, other_very_well*percent_som, 0), # estimate of somali speakers who speak english very well
    ESTsom_less = ifelse(som_total>0, other_less*percent_som, 0), # estimate of somali speakers who speak english less than very well
    pct_somali_EST = ifelse(som_total>0, somali_speakersEST/total, 0), # percent estimate of somali speakers
    EST_pct_very_well = ifelse(somali_speakersEST>0, ESTsom_very_well/somali_speakersEST, 0), # percent estimate of somali speakers who speak english very well
    EST_pct_less = ifelse(somali_speakersEST>0, ESTsom_less/somali_speakersEST, 0) # percent estimate of somali speakers who speak english less than very well
  )

## avg estimate of very well vs less than very well
total_speakers_frank <- sum(frank$somali_speakersEST, na.rm = TRUE)

total_very_well_frank <- sum(frank$ESTsom_very_well, na.rm = TRUE)
total_less_frank <- sum(frank$ESTsom_less, na.rm = TRUE)

# s_2016_frank = total somali population in franklin county
frank_very_well_16 <- (total_very_well_frank/s_2016_frank) * 100
frank_less_16 <- (total_less_frank/s_2016_frank) * 100

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
frank <- frank %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

frank$local_entropy[is.nan(frank$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
frank_totals <- frank %>%
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
frank <- frank %>%
  mutate(metro_entropy_frank = metro_entropy_frank)

# weight by somali pop

frank <- frank %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(frank$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_frank - local_entropy) / metro_entropy_frank,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
frank_entropy_index_org <- sum(frank$entropy_weight, na.rm = TRUE)
frank_entropy_index_org

## checking negative values in entropy difference (theory contradiction)
# 
f.negative_tracts <- frank %>%
  filter(entrop_difference < 0) %>%
  select(GEOID, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy_frank, entrop_difference)

# index -- summation: (adjusted)
frank <- frank %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

frank_entropy_index_adj <- sum(frank$entropy_weight_adj, na.rm = TRUE)
frank_entropy_index_adj

## percent difference between original and adjusted:
f.dif <- frank_entropy_index_adj - frank_entropy_index_org
f.p <- (f.dif / frank_entropy_index_org) * 100
f.p

## 0.288% change between original and adjusted


#######################################
# king
king <- king %>%
  mutate(
    non_som = total - som_total,
    percent_som = som_total / total
  )

# dissimilarity index (segregation)
s_2016_king <- sum(king$som_total, na.rm = TRUE)
t_2016_king <- sum(king$non_som, na.rm = TRUE)

king$dissimilarity <- abs((king$som_total/s_2016_king) - 
                             (king$non_som/t_2016_king))

king_dissim <- 0.5 * sum(king$dissimilarity, na.rm = TRUE)
king_dissim

# isolation index ()
king$isolation <- (king$som_total/s_2016_king) * (king$som_total/king$total)
king_isolation <- sum(king$isolation, na.rm = TRUE)
king_isolation

## Proportional weighting of Somali languages for each tract (Proxy estimation)
### king 2016
king <- king %>%
  mutate(
    somali_speakersEST = ifelse(som_total>0, other_lang_total*percent_som, 0), # estimate of total somali speakers
    ESTsom_very_well = ifelse(som_total>0, other_very_well*percent_som, 0), # estimate of somali speakers who speak english very well
    ESTsom_less = ifelse(som_total>0, other_less*percent_som, 0), # estimate of somali speakers who speak english less than very well
    pct_somali_EST = ifelse(som_total>0, somali_speakersEST/total, 0), # percent estimate of somali speakers
    EST_pct_very_well = ifelse(somali_speakersEST>0, ESTsom_very_well/somali_speakersEST, 0), # percent estimate of somali speakers who speak english very well
    EST_pct_less = ifelse(somali_speakersEST>0, ESTsom_less/somali_speakersEST, 0) # percent estimate of somali speakers who speak english less than very well
  )

## avg estimate of very well vs less than very well
total_speakers_king <- sum(king$somali_speakersEST, na.rm = TRUE)

total_very_well_king <- sum(king$ESTsom_very_well, na.rm = TRUE)
total_less_king <- sum(king$ESTsom_less, na.rm = TRUE)

# s_2016_king = total somali population in king county
king_very_well_16 <- (total_very_well_king/s_2016_king) * 100
king_less_16 <- (total_less_king/s_2016_king) * 100

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
king <- king %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

king$local_entropy[is.nan(king$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
king_totals <- king %>%
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
king <- king %>%
  mutate(metro_entropy_king = metro_entropy_king)

# weight by somali pop

king <- king %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(king$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_king - local_entropy) / metro_entropy_king,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
king_entropy_index_org <- sum(king$entropy_weight, na.rm = TRUE)
king_entropy_index_org

## checking negative values in entropy difference (theory contradiction)
# 
k.negative_tracts <- king %>%
  filter(entrop_difference < 0) %>%
  select(GEOID, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy_king, entrop_difference)

# index -- summation: (adjusted)
king <- king %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

king_entropy_index_adj <- sum(king$entropy_weight_adj, na.rm = TRUE)
king_entropy_index_adj

## percent difference between original and adjusted:
k.dif <- king_entropy_index_adj - king_entropy_index_org
k.p <- (k.dif / king_entropy_index_org) * 100
k.p

## 0.077% change between original and adjusted





