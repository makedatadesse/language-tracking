## Crosswalking 2020s tracts to 2010 tracts (for 2017-2021 data):

library(tidycensus)
library(tidyverse)
library(dplyr)
library(sf)

# Census Key
census_api_key("YOUR-API-KEY-HERE")

# Somali population variables (using ancestry for population)
variables <- c(
  Povertyrate = "S1701_C03_046E", # Poverty Rate
  total = "B04006_001E", # Total population
  som_total = "B04006_082E", # Somali alone
# C-Series (language proxy)
  other_lang_total = "C16001_036E", # Other & unspecified languages
  other_very_well = "C16001_037E", # Speakers who speak English "very well"
  other_less = "C16001_038E" # Speakers who speak English less than well
)

# Hennepin (2017-2021)
hen <- get_acs(
  geography = "tract",
  variables = variables,
  state = "MN",
  county = "Hennepin",
  year = 2021,
  output = "wide",
  geometry = TRUE
)

# Franklin (2017-2021)
frank <- get_acs(
  geography = "tract",
  variables = variables,
  state = "OH",
  county = "Franklin",
  year = 2021,
  output = "wide",
  geometry = TRUE
)

# King (2017-2021)
kin <- get_acs(
  geography = "tract",
  variables = variables,
  state = "WA",
  county = "King",
  year = 2021,
  output = "wide",
  geometry = TRUE
)

# combining counties for crosswalking:

# 2017-2021 Tracts:
tract_data_2021 <- bind_rows(
  hen,
  frank,
  kin
)

# removing moe columns
tract_data_2021 <- tract_data_2021 %>% select(-c(B04006_001M, B04006_082M, C16001_036M, C16001_037M, C16001_038M, S1701_C03_046M))


## Loading nhgis crosswalk data: 2020 Tracts to 2010 Tracts (ALL US)
crosswalk <- read.csv("nhgis_tr2020_tr2010.csv")



crosswalk <- crosswalk %>%
  # matching codes
  mutate(
    # Extract state, county, tract portions and recombine to match Census GEOID format
    GEOID_2020 = paste0(
      substr(tr2020gj, 2, 3),  # State code
      substr(tr2020gj, 5, 7),  # County code
      substr(tr2020gj, 9, 14)  # Tract code
    ),
    GEOID_2010 = paste0(
      substr(tr2010gj, 2, 3),  # State code
      substr(tr2010gj, 5, 7),  # County code
      substr(tr2010gj, 9, 14)  # Tract code
    )
  )

joined_data <- tract_data_2021 %>%
  left_join(crosswalk, by = c("GEOID" = "GEOID_2020"))

crosswalked_data <- joined_data %>%
  group_by(GEOID_2010) 


%>%
  summarize(
    poverty = sum(Povertyrate * wt_pop, na.rm = TRUE),
    total = sum(total * wt_pop, na.rm = TRUE),
    som_total = sum(som_total * wt_pop, na.rm = TRUE),
    other_lang_total = sum(other_lang_total * wt_pop, na.rm = TRUE),
    other_very_well = sum(other_very_well * wt_pop, na.rm = TRUE),
    other_less = sum(other_less * wt_pop, na.rm = TRUE)
  )

## formatting for readability:
options(scipen = 999)

crosswalked_data$GEOID_2010 <- as.character(crosswalked_data$GEOID_2010)
crosswalked_data <- crosswalked_data %>%
  mutate(across(where(is.numeric), ~ ifelse(. == 0, 0, round(., 4))))

## save + export
write_csv(crosswalked_data, "crosswalked.csv")



