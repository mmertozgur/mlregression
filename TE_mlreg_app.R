source("TE_mlreg_modules.R")

ui = fluidPage(
  theme= shinytheme("flatly"),
  tcgaexplorer_ml_UI("ml")
)
server= function(input,output,session) {
  dataprepml_Server("ml") 
  df =dataprepml_Server("ml")
  regularization_ml_Server("ml",regress_data  = df)
  
}

shinyApp(ui, server)