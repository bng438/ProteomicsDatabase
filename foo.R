# Installs Packages  ----
packages <- c("plyr","broom","colorspace","tidyverse",
              "plotly", "gridExtra",
              "ggseqlogo", "DT", "RDocumentation",
              "runjags","pracma", "shiny",
              "data.table","shinydashboard",
              "shinyjs", "shinyWidgets", "rsconnect",
              "BiocManager", "installr", "shinycssloaders",
              "reactlog")
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
library(rsconnect)
library(shinycssloaders)
library(reactlog)

options(shiny.reactlog = TRUE)
library(installr)
updateR()
BiocManager::install()

db <- data.table("sub_date","dataFile","date", "scientist", "protein", "cell", "description","drug","eln")
foo <- data.table("melon", "melon","melon","melon","melon","melon","jp","foo", "melon")
f <- bind_rows(db,foo)

write.csv(f,"./Proteomics Database/db.csv")

db_csv <- read.table("./Proteomics Database/db.csv",
                     header=TRUE,
                     sep=",",
                     colClasses="character"
                     )

db_csv <- db_csv[-1] %>% setDT()





db_csv[, f := foo]

t_db <- t(db) %>% as.data.frame()

for (i in 1:6)
{
  colnames(t_db)[i] <- t_db[1,i]
}



foo <- t_db[,"date"]
