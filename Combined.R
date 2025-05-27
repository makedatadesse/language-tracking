---
  title: "R Notebook"
output: html_notebook
---
  
  # Libraries and api key 
  ```{r}
library(tidyverse)
library(tidycensus)
library(dplyr)
library(cluster)
library(factoextra)
library(stringr)
library(sf)
library(ggplot2)
library(scales)
library(tibble)

# Census Key
census_api_key("3c9e13d499f9ed0c1d8f1c848dd6919fc443afad")

options(scipen=999)
```

## Somali variables:
```{r}
# Somali population variables (using ancestry for population)
variables <- c(
  poverty = "S1701_C03_046E", # Poverty Rate
  total = "B04006_001E", # Total population
  som_total = "B04006_082E", # Somali alone
  # C-Series (language proxy)
  other_lang_total = "C16001_036E", # Other & unspecified languages
  other_very_well = "C16001_037E", # Speakers who speak English "very well"
  other_less = "C16001_038E" # Speakers who speak English less than well
)
```

## Loading in 2011-2016 county data as separate tables (and cleaning)
```{r}
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
hen <- hen %>% select(-c(B04006_001M, B04006_082M, C16001_036M, C16001_037M, C16001_038M, S1701_C03_046M))

frank <- frank %>% select(-c(B04006_001M, B04006_082M, C16001_036M, C16001_037M, C16001_038M, S1701_C03_046M))

king <- king %>% select(-c(B04006_001M, B04006_082M, C16001_036M, C16001_037M, C16001_038M, S1701_C03_046M))
```

## Loading in 2017-2021 crosswalked data

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
```

## Calculating indices for each county + tract
# (Dissimilarity, Isolation, Entropy, Ratios, etc.)

```{r}
## Entropy score formula (multigroup): E = -∑(pi × ln(pi)) 
calc_entropy <- function(props) {
  # Remove NA or 0 values since log(0) is undefined
  props <- props[!is.na(props) & props > 0]
  if (length(props) == 0) return(0)
  -sum(props * log(props))
}
```

### HENNEPIN 2011-2016
```{r}
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

hennepin_dissim_16 <- 0.5 * sum(hen$dissimilarity, na.rm = TRUE)
hennepin_dissim_16

# isolation index ()
hen$isolation <- (hen$som_total/s_2016_hen) * (hen$som_total/hen$total)
hennepin_isolation_16 <- sum(hen$isolation, na.rm = TRUE)
hennepin_isolation_16

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
total_speakers_hen_16 <- sum(hen$somali_speakersEST, na.rm = TRUE)

total_very_well_hen_16 <- sum(hen$ESTsom_very_well, na.rm = TRUE)
total_less_hen_16 <- sum(hen$ESTsom_less, na.rm = TRUE)

# s_2016_hen = total somali population in hennepin county
hennepin_very_well_16 <- (total_very_well_hen_16/s_2016_hen) * 100
hennepin_less_16 <- (total_less_hen_16/s_2016_hen) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
hen <- hen %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

hen$local_entropy[is.nan(hen$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
hen_totals <- hen %>%
  summarize(
    total_speakers_hen_16,
    total_very_well_hen_16,
    total_less_hen_16,
    # proportions:
    p_very_well = total_very_well_hen_16 / total_speakers_hen_16,
    p_less = total_less_hen_16 / total_speakers_hen_16
  )

metro_entropy_hen_16 <- with(hen_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_hen_16

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
hen <- hen %>%
  mutate(metro_entropy_hen_16 = metro_entropy_hen_16)

# weight by somali pop
hen <- hen %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(hen$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_hen_16 - local_entropy) / metro_entropy_hen_16,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
hen_entropy_index_org_16 <- sum(hen$entropy_weight, na.rm = TRUE)
hen_entropy_index_org_16

## checking negative values in entropy difference (theory contradiction)
h.negative_tracts_16 <- hen %>%
  filter(entrop_difference < 0) %>%
  select(GEOID, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy_hen_16, entrop_difference)

# index -- summation: (adjusted)
hen <- hen %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

hen_entropy_index_adj_16 <- sum(hen$entropy_weight_adj, na.rm = TRUE)
hen_entropy_index_adj_16

## percent difference between original and adjusted:
h.dif.16 <- hen_entropy_index_adj_16 - hen_entropy_index_org_16
h.p.16 <- (h.dif.16 / hen_entropy_index_org_16) * 100
h.p.16

## 0.784% change between original and adjusted

# RELATIVE QUOTIENT Hennepin 2016:
hen_2016 <- hen %>%
  mutate(
    # Relative concentration
    county_pct_somali = s_2016_hen/t_2016_hen,
    relative_quotient = (percent_som / county_pct_somali),
    # Classify tracts based on relative concentration
    concentration_category = case_when(
      relative_quotient >= 2 ~ "High",     # 2x or more county average
      relative_quotient <= 0.5 ~ "Low",    # Half or less county average
      TRUE ~ "Moderate"
    )
  )

proficiency_model <- lm(EST_pct_very_well ~ relative_quotient + poverty + total,
                        data = hen_2016)
summary(proficiency_model)

category_summary <- hen_2016 %>%
  group_by(concentration_category) %>%
  summarize(
    avg_proficiency = mean(EST_pct_very_well, na.rm = TRUE),
    median_proficiency = median(EST_pct_very_well, na.rm = TRUE),
    n_tracts = n()
  )

```

### HENNEPIN 2017-2021
```{r}
# dissimilarity index (segregation)
s_2021_hen <- sum(hennepin_data$som_total, na.rm = TRUE)
t_2021_hen <- sum(hennepin_data$non_som, na.rm = TRUE)

hennepin_data$dissimilarity <- abs((hennepin_data$som_total/s_2021_hen) - 
                                     (hennepin_data$non_som/t_2021_hen))

hennepin_dissim_21 <- 0.5 * sum(hennepin_data$dissimilarity, na.rm = TRUE)
hennepin_dissim_21

# isolation index ()
hennepin_data$isolation <- (hennepin_data$som_total/s_2021_hen) * (hennepin_data$som_total/hennepin_data$total)
hennepin_isolation_21 <- sum(hennepin_data$isolation, na.rm = TRUE)
hennepin_isolation_21

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
total_speakers_hen_21 <- sum(hennepin_data$somali_speakersEST, na.rm = TRUE)

total_very_well_hen_21 <- sum(hennepin_data$ESTsom_very_well, na.rm = TRUE)
total_less_hen_21 <- sum(hennepin_data$ESTsom_less, na.rm = TRUE)

# s_2021_hen = total somali population in hennepin county
hennepin_very_well_21 <- (total_very_well_hen_21/s_2021_hen) * 100
hennepin_less_21 <- (total_less_hen_21/s_2021_hen) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
hennepin_data <- hennepin_data %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

hennepin_data$local_entropy[is.nan(hennepin_data$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
hennepin_totals <- hennepin_data %>%
  summarize(
    total_speakers_hen_21,
    total_very_well_hen_21,
    total_less_hen_21,
    # proportions:
    p_very_well = total_very_well_hen_21 / total_speakers_hen_21,
    p_less = total_less_hen_21 / total_speakers_hen_21
  )

metro_entropy_hen_21 <- with(hennepin_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_hen_21

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
hennepin_data <- hennepin_data %>%
  mutate(metro_entropy_hen_21 = metro_entropy_hen_21)

# weight by somali pop
hennepin_data <- hennepin_data %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(hennepin_data$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_hen_21 - local_entropy) / metro_entropy_hen_21,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
hen_entropy_index_org_21 <- sum(hennepin_data$entropy_weight, na.rm = TRUE)
hen_entropy_index_org_21

## checking negative values in entropy difference (theory contradiction)
h.negative_tracts.21 <- hennepin_data %>%
  filter(entrop_difference < 0) %>%
  select(GEOID_2010, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy_hen_21, entrop_difference)

# index -- summation: (adjusted)
hennepin_data <- hennepin_data %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

hen_entropy_index_adj_21 <- sum(hennepin_data$entropy_weight_adj, na.rm = TRUE)
hen_entropy_index_adj_21

## 1.58% change between original and adjusted

# RELATIVE QUOTIENT AND UPDATED STRING: Hennepin 2021:
hen_2021 <- hennepin_data %>%
  mutate(
    # Relative concentration
    county_pct_somali = s_2021_hen/t_2021_hen,
    relative_quotient = (percent_som / county_pct_somali),
    # Classify tracts based on relative concentration
    concentration_category = case_when(
      relative_quotient >= 2 ~ "High",     # 2x or more county average
      relative_quotient <= 0.5 ~ "Low",    # Half or less county average
      TRUE ~ "Moderate"
    )
  )

proficiency_model_21 <- lm(EST_pct_very_well ~ relative_quotient + poverty + total,
                           data = hen_2021)
summary(proficiency_model_21)

hennepin_data$proficiency_change <- hennepin_data$EST_pct_very_well - hen$EST_pct_very_well

model.test <- lm(proficiency_change ~ percent_som + isolation + local_entropy + poverty, data = hennepin_data)
summary(model.test)

cor(hennepin_data$isolation, hennepin_data$pct_somali_EST, use = "complete.obs", method = "pearson")

```

### FRANKLIN 2011-2016
```{r}
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

frank_dissim_16 <- 0.5 * sum(frank$dissimilarity, na.rm = TRUE)
frank_dissim_16

# isolation index ()
frank$isolation <- (frank$som_total/s_2016_frank) * (frank$som_total/frank$total)
franklin_isolation_16 <- sum(frank$isolation, na.rm = TRUE)
franklin_isolation_16

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
total_speakers_frank_16 <- sum(frank$somali_speakersEST, na.rm = TRUE)

total_very_well_frank_16 <- sum(frank$ESTsom_very_well, na.rm = TRUE)
total_less_frank_16 <- sum(frank$ESTsom_less, na.rm = TRUE)

# s_2016_frank = total somali population in franklin county
frank_very_well_16 <- (total_very_well_frank_16/s_2016_frank) * 100
frank_less_16 <- (total_less_frank_16/s_2016_frank) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
frank <- frank %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

frank$local_entropy[is.nan(frank$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
frank_totals <- frank %>%
  summarize(
    total_speakers_frank_16,
    total_very_well_frank_16,
    total_less_frank_16,
    # proportions:
    p_very_well = total_very_well_frank_16 / total_speakers_frank_16,
    p_less = total_less_frank_16 / total_speakers_frank_16
  )

metro_entropy_frank_16 <- with(frank_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_frank_16

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
frank <- frank %>%
  mutate(metro_entropy_frank_16 = metro_entropy_frank_16)

# weight by somali pop
frank <- frank %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(frank$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_frank_16 - local_entropy) / metro_entropy_frank_16,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
frank_entropy_index_org_16 <- sum(frank$entropy_weight, na.rm = TRUE)
frank_entropy_index_org_16

## checking negative values in entropy difference (theory contradiction)
f.negative_tracts.16 <- frank %>%
  filter(entrop_difference < 0) %>%
  select(GEOID, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy_frank_16, entrop_difference)

# index -- summation: (adjusted)
frank <- frank %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

frank_entropy_index_adj_16 <- sum(frank$entropy_weight_adj, na.rm = TRUE)
frank_entropy_index_adj_16

## percent difference between original and adjusted:
f.dif.16 <- frank_entropy_index_adj_16 - frank_entropy_index_org_16
f.p.16 <- (f.dif.16 / frank_entropy_index_org_16) * 100
f.p.16

## 0.288% change between original and adjusted

# RELATIVE QUOTIENT AND UPDATED STRING: Frank 2016:
frank_2016 <- frank %>%
  mutate(
    # Relative concentration
    county_pct_somali = s_2016_frank/t_2016_frank,
    relative_quotient = (percent_som / county_pct_somali),
    # Classify tracts based on relative concentration
    concentration_category = case_when(
      relative_quotient >= 2 ~ "High",     # 2x or more county average
      relative_quotient <= 0.5 ~ "Low",    # Half or less county average
      TRUE ~ "Moderate"
    )
  )

f.prof_model_16 <- lm(EST_pct_very_well ~ relative_quotient + poverty + total,
                      data = frank_2016)
summary(f.prof_model_16)


```

### FRANKLIN 2017-2021
```{r}
# dissimilarity index (segregation)
s_2021_frank <- sum(frank_data$som_total, na.rm = TRUE)
t_2021_frank <- sum(frank_data$non_som, na.rm = TRUE)

frank_data$dissimilarity <- abs((frank_data$som_total/s_2021_frank) - 
                                  (frank_data$non_som/t_2021_frank))

franklin_dissim_21 <- 0.5 * sum(frank_data$dissimilarity, na.rm = TRUE)
franklin_dissim_21

# isolation index ()
frank_data$isolation <- (frank_data$som_total/s_2021_frank) * (frank_data$som_total/frank_data$total)
franklin_isolation_21 <- sum(frank_data$isolation, na.rm = TRUE)
franklin_isolation_21

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
total_speakers_frank_21 <- sum(frank_data$somali_speakersEST, na.rm = TRUE)

total_very_well_frank_21 <- sum(frank_data$ESTsom_very_well, na.rm = TRUE)
total_less_frank_21 <- sum(frank_data$ESTsom_less, na.rm = TRUE)

# s_2021_frank = total somali population in franklin county
franklin_very_well_21 <- (total_speakers_frank_21/s_2021_frank) * 100
franklin_less_21 <- (total_less_frank_21/s_2021_frank) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
frank_data <- frank_data %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

frank_data$local_entropy[is.nan(frank_data$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
franklin_totals <- frank_data %>%
  summarize(
    total_speakers_frank_21,
    total_very_well_frank_21,
    total_less_frank_21,
    # proportions:
    p_very_well = total_very_well_frank_21 / total_speakers_frank_21,
    p_less = total_less_frank_21 / total_speakers_frank_21
  )

metro_entropy_frank_21 <- with(franklin_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_frank_21

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
frank_data <- frank_data %>%
  mutate(metro_entropy_frank_21 = metro_entropy_frank_21)

# weight by somali pop
frank_data <- frank_data %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(frank_data$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_frank_21 - local_entropy) / metro_entropy_frank_21,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
frank_entropy_index_org_21 <- sum(frank_data$entropy_weight, na.rm = TRUE)
frank_entropy_index_org_21

# index -- summation: (adjusted)
frank_data <- frank_data %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

frank_entropy_index_adj_21 <- sum(frank_data$entropy_weight_adj, na.rm = TRUE)
frank_entropy_index_adj_21

# See how much the adjustment changes the overall index
f.dif.21 <- frank_entropy_index_adj_21 - frank_entropy_index_org_21
f.p.21 <- (f.dif.21 / frank_entropy_index_org_21) * 100
f.p.21

## 4.36% change between original and adjusted

# RELATIVE QUOTIENT AND UPDATED STRING: Frank 2021:
frank_2021 <- frank_data %>%
  mutate(
    # Relative concentration
    county_pct_somali = s_2021_frank/t_2021_frank,
    relative_quotient = (percent_som / county_pct_somali),
    # Classify tracts based on relative concentration
    concentration_category = case_when(
      relative_quotient >= 2 ~ "High",     # 2x or more county average
      relative_quotient <= 0.5 ~ "Low",    # Half or less county average
      TRUE ~ "Moderate"
    )
  )

f.prof_model_21 <- lm(EST_pct_very_well ~ relative_quotient + poverty + total,
                      data = frank_2021)
summary(f.prof_model_21)


```

### KING 2011-2016
```{r}
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

king_dissim_16 <- 0.5 * sum(king$dissimilarity, na.rm = TRUE)
king_dissim_16

# isolation index ()
king$isolation <- (king$som_total/s_2016_king) * (king$som_total/king$total)
king_isolation_16 <- sum(king$isolation, na.rm = TRUE)
king_isolation_16

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
total_speakers_king_16 <- sum(king$somali_speakersEST, na.rm = TRUE)

total_very_well_king_16 <- sum(king$ESTsom_very_well, na.rm = TRUE)
total_less_king_16 <- sum(king$ESTsom_less, na.rm = TRUE)

# s_2016_king = total somali population in king county
king_very_well_16 <- (total_very_well_king_16/s_2016_king) * 100
king_less_16 <- (total_less_king_16/s_2016_king) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
king <- king %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

king$local_entropy[is.nan(king$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
k_totals <- king %>%
  summarize(
    total_speakers_king_16,
    total_very_well_king_16,
    total_less_king_16,
    # proportions:
    p_very_well = total_very_well_king_16 / total_speakers_king_16,
    p_less = total_less_king_16 / total_speakers_king_16
  )

metro_entropy_king_16 <- with(k_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_king_16

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
king <- king %>%
  mutate(metro_entropy_king_16 = metro_entropy_king_16)

# weight by somali pop
king <- king %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(king$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_king_16 - local_entropy) / metro_entropy_king_16,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
king_entropy_index_org_16 <- sum(king$entropy_weight, na.rm = TRUE)
king_entropy_index_org_16

## checking negative values in entropy difference (theory contradiction)
k.negative_tracts.16 <- king %>%
  filter(entrop_difference < 0) %>%
  select(GEOID, somali_speakersEST, EST_pct_very_well, EST_pct_less, local_entropy, metro_entropy_king_16, entrop_difference)

# index -- summation: (adjusted)
king <- king %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

king_entropy_index_adj_16 <- sum(king$entropy_weight_adj, na.rm = TRUE)
king_entropy_index_adj_16

## percent difference between original and adjusted:
k.dif.16 <- king_entropy_index_adj_16 - king_entropy_index_org_16
k.p.16 <- (k.dif.16 / king_entropy_index_org_16) * 100
k.p.16

## 0.077% change between original and adjusted

# RELATIVE QUOTIENT AND UPDATED STRING: King 2016:
king_2016 <- king %>%
  mutate(
    # Relative concentration
    county_pct_somali = s_2016_king/t_2016_king,
    relative_quotient = (percent_som / county_pct_somali),
    # Classify tracts based on relative concentration
    concentration_category = case_when(
      relative_quotient >= 2 ~ "High",     # 2x or more county average
      relative_quotient <= 0.5 ~ "Low",    # Half or less county average
      TRUE ~ "Moderate"
    )
  )

k.prof_model_16 <- lm(EST_pct_very_well ~ relative_quotient + poverty + total,
                      data = king_2016)
summary(k.prof_model_16)

```

### KING 2017-2021
```{r}
# dissimilarity index (segregation)
s_2021_king <- sum(king_data$som_total, na.rm = TRUE)
t_2021_king <- sum(king_data$non_som, na.rm = TRUE)

king_data$dissimilarity <- abs((king_data$som_total/s_2021_king) - 
                                 (king_data$non_som/t_2021_king))

king_dissim_21 <- 0.5 * sum(king_data$dissimilarity, na.rm = TRUE)
king_dissim_21

# isolation index ()
king_data$isolation <- (king_data$som_total/s_2021_king) * (king_data$som_total/king_data$total)
king_isolation_21 <- sum(king_data$isolation, na.rm = TRUE)
king_isolation_21

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
total_speakers_king_21 <- sum(king_data$somali_speakersEST, na.rm = TRUE)

total_very_well_king_21 <- sum(king_data$ESTsom_very_well, na.rm = TRUE)
total_less_king_21 <- sum(king_data$ESTsom_less, na.rm = TRUE)

# s_2021_king = total somali population in franklin county
king_very_well_21 <- (total_speakers_king_21/s_2021_king) * 100
king_less_21 <- (total_less_king_21/s_2021_king) * 100

## Check 2: entropy index along linguistic variables (english proficency)
# "The Dimensions of Residential Segregation" Denton + Massey
# Census: "Housing Patterns: Appendix B: Measures of Residential Segregation"
# "The Multigroup Entropy Index" Iceland

### Local [tract] level entropy (binary formula): Ei = (pi)log[ 1/pi] + (1 — pi)1og[1/(1— pi)]
king_data <- king_data %>%
  mutate(
    local_entropy = EST_pct_very_well * log(1/EST_pct_very_well) + EST_pct_less * log(1/EST_pct_less)
  )

king_data$local_entropy[is.nan(king_data$local_entropy)] <- 0.00000000 # changing NaN to 0s

### Metro-wide entropy: E = (P)log[1/P] + (1— P)log[1/(1 — P)]
king_totals <- king_data %>%
  summarize(
    total_speakers_king_21,
    total_very_well_king_21,
    total_less_king_21,
    # proportions:
    p_very_well = total_very_well_king_21 / total_speakers_king_21,
    p_less = total_less_king_21 / total_speakers_king_21
  )

metro_entropy_king_21 <- with(king_totals, calc_entropy(c(p_very_well, p_less)))
metro_entropy_king_21

### Entropy Index (H = ∑[(ti(E - Ei))/ET] or H = ∑[(ti/T) × (E - Ei)/E])
king_data <- king_data %>%
  mutate(metro_entropy_king_21 = metro_entropy_king_21)

# weight by somali pop
king_data <- king_data %>%
  mutate(
    # ti/T: "ti = the total population of area i, T= the sum of all ti"
    weight_tract = somali_speakersEST / sum(king_data$somali_speakersEST),
    # (E - Ei)/E: E = metro-wide entropy score, Ei=tract level entropy score
    entrop_difference = (metro_entropy_king_21 - local_entropy) / metro_entropy_king_21,
    # form result (weighted)
    entropy_weight = weight_tract * entrop_difference
  )

# index -- summation: (original)
king_entropy_index_org_21 <- sum(king_data$entropy_weight, na.rm = TRUE)
king_entropy_index_org_21

# index -- summation: (adjusted)
king_data <- king_data %>%
  mutate(
    entrop_difference_adj = pmax(0, entrop_difference),
    entropy_weight_adj = weight_tract * entrop_difference_adj
  )

k.negative_tracts.21 <- king_data %>%
  filter(entrop_difference < 0) %>%
  select(GEOID_2010, somali_speakersEST, EST_pct_very_well, EST_pct_less, 
         local_entropy, metro_entropy_king_21, entrop_difference)
k.negative_tracts.21

king_entropy_index_adj_21 <- sum(king_data$entropy_weight_adj, na.rm = TRUE)
king_entropy_index_adj_21

# See how much the adjustment changes the overall index
k.dif.21 <- king_entropy_index_adj_21 - king_entropy_index_org_21
k.p.21 <- (k.dif.21 / king_entropy_index_org_21) * 100
k.p.21

## 0.819% change between original and adjusted

# RELATIVE QUOTIENT AND UPDATED STRING: King 2021:
king_2021 <- king_data %>%
  mutate(
    # Relative concentration
    county_pct_somali = s_2021_king/t_2021_king,
    relative_quotient = (percent_som / county_pct_somali),
    # Classify tracts based on relative concentration
    concentration_category = case_when(
      relative_quotient >= 2 ~ "High",     # 2x or more county average
      relative_quotient <= 0.5 ~ "Low",    # Half or less county average
      TRUE ~ "Moderate"
    )
  )

k.prof_model_21 <- lm(EST_pct_very_well ~ relative_quotient + poverty + total,
                      data = king_2021)
summary(k.prof_model_21)

```

### Summary Results Table + ggplot 1 ("Average Isolation Index vs Change in English Proficiency")
```{r}
df <- tibble(
  county = c("Hennepin", "Franklin", "King"),
  average_isolation = c((hennepin_isolation_16 + hennepin_isolation_21) / 2, (franklin_isolation_16 + franklin_isolation_21) / 2, (king_isolation_16 + king_isolation_21) / 2),
  language_shift = c((hennepin_totals$p_very_well - hen_totals$p_very_well) * 100, (franklin_totals$p_very_well - frank_totals$p_very_well) * 100, (king_totals$p_very_well - k_totals$p_very_well) * 100)
)

df.detailed <- data_frame(
  Period = c("2011-2016", "2017-2021", "2011-2016", "2017-2021", "2011-2016", "2017-2021"),
  County = c("Hennepin", "Hennepin", "Franklin", "Franklin", "King", "King"),
  Somali_Ratio = c((s_2016_hen/t_2016_hen)*100, (s_2021_hen/t_2021_hen)*100, 
                   (s_2016_frank/t_2016_frank)*100, (s_2021_frank/t_2021_frank)*100,
                   (s_2016_king/t_2016_king)*100, (s_2021_king/t_2021_king)*100),
  Dissimilarity_Index = c(hennepin_dissim_16, hennepin_dissim_21, frank_dissim_16, 
                          franklin_dissim_21, king_dissim_16, king_dissim_21),
  Isolation_index = c(hennepin_isolation_16, hennepin_isolation_21, franklin_isolation_16, 
                      franklin_isolation_21, king_isolation_16, king_isolation_21),
  Proficiency_change = c(" ", (hennepin_totals$p_very_well - hen_totals$p_very_well)*100,
                         " ", (franklin_totals$p_very_well - frank_totals$p_very_well)*100,
                         " ", (king_totals$p_very_well - k_totals$p_very_well)*100)
)

df.detailed <- df.detailed %>%
  mutate(
    concentration_score = (scale(Dissimilarity_index) + scale(Isolation_index) +
                             scale(Somali_Ratio))/3,
    concentration_category = case_when(
      concentration_score > 0.5 ~ "High",
      concentration_score < -0.2 ~ "Low",
      TRUE ~ "Moderate"
    )
  )

ggplot(df, aes(average_isolation, language_shift, label=county)) +
  geom_point() + geom_smooth(method="lm", se=FALSE) +
  geom_text(vjust=-1) +
  labs(x="Average Isolation Index (2011-2021)", y="English “Very Well” Δ (pp)", 
       title="County-Level Concentration vs. Language Shift")

```

### Joined Tract Table
```{r}
## Rename GEOID_2010 to GEOID:
hennepin_data <- hennepin_data %>% 
  rename(
    GEOID = GEOID_2010
  )
hennepin_data$GEOID <- as.character(hennepin_data$GEOID)

frank_data <- frank_data %>% 
  rename(
    GEOID = GEOID_2010
  )                                               
frank_data$GEOID <- as.character(frank_data$GEOID)


king_data <- king_data %>% 
  rename(
    GEOID = GEOID_2010
  )
king_data$GEOID <- as.character(king_data$GEOID)

prep_county <- function(df16, df21, county_name){
  df16 <- df16 %>%
    transmute(
      GEOID,
      county = county_name,
      pct_som_16 = percent_som,
      isol_16 = isolation,
      diss_16 = dissimilarity,
      pov_16 = poverty,
      ent_16 = local_entropy,
      prof16 = EST_pct_very_well,
      quot16 = relative_quotient
    )
  
  df21 <- df21 %>%
    transmute(
      GEOID,
      pct_som_21 = percent_som,
      isol_21 = isolation,
      diss_21 = dissimilarity,
      pov_21 = poverty,
      ent_21 = local_entropy,
      prof21 = EST_pct_very_well,
      quot21 = relative_quotient
    )
  
  df16 %>%
    inner_join(df21, by = "GEOID") %>%
    mutate(
      proficiency_change = prof21 - prof16
    ) %>%
    select(
      GEOID, county,
      pct_som_16, pct_som_21,
      isol_16, isol_21,
      diss_16, diss_21,
      pov_16, pov_21,
      ent_16, ent_21,
      prof16, prof21,
      quot16, quot21,
      proficiency_change
    )
}

hennepin_ct <- prep_county(hen_2016, hen_2021, "Hennepin")
franklin_ct <- prep_county(frank_2016, frank_2021, "Franklin")
king_ct <- prep_county(king_2016, king_2021, "King")

# combining tracts
all_tracts <- bind_rows(hennepin_ct, franklin_ct, king_ct)

model_pooled <- lm(proficiency_change ~ pct_som_15 + isol_15 + pov_15 + ent_15 + county,
                   data = all_tracts)
summary(model_pooled)



```


