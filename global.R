# Brandon Ng
# 
# Global fields for registering files into database
# 

library(plyr)
library(tidyverse)
library(plotly)
library(readxl)
library(gridExtra)
library(ggseqlogo)
library(DT)
library(runjags)
library(pracma)
library(shiny)
library(data.table)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)

# Global fields  ----
proteins <- c("Protein A", "Protein B", "Protein D", "Protein C") %>% sort()
cells <- c("THP-1", "EL4", "RA1", "CT26", "A375") %>% sort()
scientists <- c("Brandon", "Gramoz", "Tom", "Vanessa") %>% sort()
drugs <- c("123","456","12424") %>% sort()
fields <- c("dataFile","date", "scientist", "protein", "cell", "description","drug", "eln")
empty_fields <- c("dataFile","scientist", "protein", "cell", "description","drug","eln")