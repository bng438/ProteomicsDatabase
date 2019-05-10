# Brandon Ng
# 
# Shiny app for registering files into a database
# 


# Loads Packages  ----
library(plyr)
library(broom)
library(colorspace)
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
library(xfun)
library(shinycssloaders)
library(rsconnect)
source("server.R")
source("ui.R")
source("global.R")


# Runs the app   ----
shinyApp(ui=ui,server=server)