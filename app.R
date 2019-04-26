# Brandon Ng
# 
# Shiny app for registering files into a database
# 


# Loads Packages  ----
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
source("server.R")
source("ui.R")
source("global.R")


# Runs the app   ----
shinyApp(ui=ui,server=server)