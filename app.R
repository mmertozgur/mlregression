# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# library(BiocManager)
# options(repos = BiocManager::repositories())
# install.packages("dplyr")
# BiocManager::install("edgeR")
# install.packages("tidyverse")
# install.packages("glmnet")
# BiocManager::install("msigdb")
# install.packages("shiny")
# install.packages("shinythemes")
# install.packages("DT")
# install.packages("broom")
# install.packages("ggplot2")
# install.packages("plotly")
library(dplyr)
library(edgeR)
library(tidyverse)
library(glmnet)
library(msigdbr)
library(shiny)
library(shinythemes)
library(DT)
library(broom)
library(ggplot2)
library(plotly)

source("TE_ml.regression_modules_final.R")
ui = fluidPage(
  theme= shinytheme("flatly"),
  tcgaexplorer_ml_UI("ml")
)
server= function(input,output,session) {
  dataprepml_Serverz("ml")
  
  df = dataprepml_Serverz("ml")
  regression_ml_Server("ml",regress_data  = df)
  
}

shinyApp(ui, server)