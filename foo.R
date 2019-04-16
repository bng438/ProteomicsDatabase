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

db <- data.table(c("a","b","c","d"))
db <- data.table(c("form_data","date", "scientist", "protein", "cell", "description","drug"))

write.csv(db,"./Proteomics Database/db.csv")

db_csv <- read.csv("./Proteomics Database/db.csv")
db_csv <- db_csv[-1] %>% setDT()
foo <- c("melon","melon","melon","melon","melon","jp")

db_csv[, f := foo]

t_db <- t(db) %>% as.data.frame()

for (i in 1:6)
{
  colnames(t_db)[i] <- t_db[1,i]
}



foo <- t_db[,"date"]
