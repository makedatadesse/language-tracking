# language-tracking
## Tracking Immigrant Acculturation Patterns Through Ethnic Concentrations and English Language Adoption

## Description:
The goal of this project is to visualize and quantify language adoption as an aspect of acculturation by spatially and statistically exploring the relationship between ethnic enclaves and English language proficiency trends by focusing on Somali immigrants in the United States (using ACS time series data).

## Rationale:
Classic spatial segregation theories suggest that ethnic enclaves prevent or slow language proficiency; there is limited research that focuses on African immigrant populations at the tract level, and studies their communities through socio-cultural trends. This study, though hyper-focused to one specific ethnic and language group, challenges these theories by directly testing whether tightly clustered Somali populations demonstrate slower or faster growth in English proficiency.

## Workflow/Required Libraries:
Almost all the preliminary workflow was done in RStudio (API authenticating, crosswalking census tracts, creating custom data frames and cleaning ACS variables, and exporting csv files). Request an API key from the [Census](https://api.census.gov/data/key_signup.html)  
***Many libraries were used in this project only as supplements, the most notable and important to run this code:***  
```
  tidyverse (to run tidycensus)  
  tidycensus (to access census api)  
  dplyr  
  lme4 (for linear mixed modeling)  
```

ArcGIS Pro was utilized to run geoprocessing tools (hotspots) to identify enclaves and then create maps based on cleaned data (as well as bridge GEOIDs and FIPS codes)

_Thiels H entropy, dissimilarity, and isolation indices were all used at the tract and county level in addition to regression for analysis_

## Further uses:
The statistical analysis of this project, though illuminating, does face some limitations. To bypass issues with data granularity, ACS language variables were proxied. This project is scalable and a possible base model for investigating other immigrant language proficiency and acculturation patterns at larger scales and across longer time periods.
This project may be beneficial for non-profit resource allocation, community outreach program direction, cultural heritage institutions, and socio-cultural anthropology researchers and artists.

## Files list:
2011to2016.R – analysis of tract and county level data for 2011-2016  
2017to2021.R – analysis of tract and county level data for 2017-2021  
Combined.R – combining both of the above into codebook to run trend analysis  
Crosswalking.R – crosswalking of needed variables and tracts using NHGIS file  
FullCodeBook.Rmd – full cleaned codebook for analysis and visualization in its entirety. Pulled from for ArcGIS  
crosswalked.csv – updated crosswalked tracts with needed variables  
nhgis_tr2020_tr2010.csv – csv file from [NHGIS](https://www.nhgis.org/geographic-crosswalks)  

## ArcGIS Hotspot Map 
_county specific maps can be found in published paper_

![National Tract Level Hotspot Map where Somali Concentration exceeds the county average. National map is predominantly gray with tracts that have a greater somali concentration that their countys average are highlighted in yellow to red, red being the highest. To the left are three county's that are zoomed in and are most if not fully red. In order these county's are Hennepin County, Franklin County, and King County.](https://github.com/user-attachments/assets/9a77a877-16be-4757-8e2e-e3f3bfe99c1f)
