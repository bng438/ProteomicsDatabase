# Brandon Ng
# 
# Shiny app for registering files into a database
# 


# Installs Packages  ----
packages <- c("plyr","tidyverse","plotly", 
              "gridExtra","ggseqlogo","DT",
              "RDocumentation","runjags","pracma",
              "shiny","data.table","shinydashboard",
              "shinyjs", "shinyWidgets")
for (i in seq_along(packages))
{
  if(!requireNamespace(packages[i]))
  {
    install.packages(packages[i])
  }
}


# Clears Global Workspace  ----
rm(list = ls())


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