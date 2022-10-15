library(dplyr)
library(edgeR)
library(tidyverse)
library(glmnet)
library(msigdbr)
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(cluster)
library(ClusterR)
library(factoextra)
library(broom)


BRCA_metadata = readRDS("~/Desktop/miRegularization/miregularization/BRCA_metadata.RDS")
hallmark_sets = as.data.frame(msigdbr(species = "Homo sapiens", category = "H" ))
hallmark_sets = hallmark_sets %>% select(c(gs_name,gene_symbol))
hsets_names4freq = as.data.frame(table(hallmark_sets %>% select(gs_name)))


zero_adjuster = function(mdata = BRCA_metadata,index= 3:59000, max_zero_percent) {
  mdata_genes = mdata[index]
  zero_freqmaatrix = t(as.matrix(lapply(mdata_genes, function(x){ length(which(x==0))/length(x) *100 })))
  available_genes = zero_freqmaatrix[,zero_freqmaatrix[]  < max_zero_percent]
  available_genelist = names(available_genes)
  
  return(as.data.frame(available_genelist))
}


tcgaexplorer_ml_UI = function(id) {
  ns = NS(id)
  navbarPage(
    "TCGA Explorer ML",
    tabPanel(
      "Data&Prep",
      sidebarPanel(
        selectizeInput(NS(id,"tcga_subtype"),"Select a TCGA project",choices = c("BRCA_metadata.RDS")),
        selectizeInput(NS(id,"predictor_set"), "Select a predictor set", choices = c("miRNA_BRCA")),
        selectizeInput(ns("response_type"), "Select a data preparation pipeline for the response set",
                       choices = c( "GSEA Hallmark Gene Sets" ="hallmark", "Gene List" = "selgene")),
        conditionalPanel(
          condition = "input.response_type == 'hallmark'",ns = ns,
          selectizeInput(NS(id,"hm_name"), "Select a gene set from GSEA hallmark sets as a response set",
                         choices = hsets_names4freq$gs_name)
        ),
        conditionalPanel(
          condition = "input.response_type == 'selgene' ", ns=ns,
          sliderInput(NS(id,"max_zero_percent"),"Set the maximum percentage of zero count",
                      min = 0, max = 100, value = 100, step = 1),
          selectizeInput(NS(id,"gene_list"),"Select gene(s)", choices = c(""), multiple = TRUE)
        )
      ),
      mainPanel(
        textOutput(NS(id,"selected_genes")),
        dataTableOutput(NS(id,"response_genes")),
        dataTableOutput(NS(id,"reg_df")),
        conditionalPanel(
          condition = "input.response_type == 'selgene' ", ns=ns,
          dataTableOutput(NS(id,"ava_genes_table"))
        )
      )
    ),
    tabPanel(
      "Regularized Regression",
      sidebarPanel(
        sliderInput(NS(id,"user_alpha"),"Enter an alpha parameter to set the regression technique", min = 0, max = 1, value = 1, step = 0.2),
        selectizeInput(NS(id,"lambda_for_coef"), "Select lambda value at which model coefficients are being displayed: ", 
                       choices = c("lambda.min","lambda.1se")),
        actionButton("analyze_button", "Analyze!")
      ),
      mainPanel(
        fluidRow(
          column(6,plotOutput(NS(id,"lambda_error_plot"))),
          column(6,plotOutput(NS(id,"coef_lambda_plot")))
        ),
        fluidRow(
          column(4,
                 dataTableOutput(NS(id,"coef_data"))),
          column(2,
                 textOutput(NS(id,"lambda_value_min")),
                 textOutput(NS(id,"lambda_value_1se"))
                 #textInput(NS(id,"cluster_number"), label = "Enter the number of centroids for the KM Clustering:")
          ),
          column(6,
                 textOutput(NS(id,"cf_summary"))),
          plotOutput(NS(id,"try_plot"))
        )
      )
    ) #insert comma here when adding the tab!
    # tabPanel(
    #   "work",
    #   sidebarPanel(),
    #   mainPanel(
    #     dataTableOutput(NS(id,"work_table"))
    #   )
    # )
  )
}

dataprepml_Server = function(id) {
  moduleServer(id,function(input,output,session){
    
    
    gene_list_hallmark = reactive({filter(hallmark_sets, gs_name == input$hm_name) %>% select(gene_symbol)})
    gene_list_selected_df = reactive({as.data.frame(input$gene_list)})
    
    available_genelist = reactive({zero_adjuster(index = 3:5000, max_zero_percent =  input$max_zero_percent)})
    
    observeEvent(input$max_zero_percent, {
      percent_info = reactive({as.numeric(input$max_zero_percent)})
      updateSelectizeInput(inputId = "gene_list", choices = available_genelist())
    })
    
    #Returns a (user-def.) gene-list of the response variable determinants.  
    glist_r = reactive({
      if(input$response_type == 'hallmark') {
        glist = gene_list_hallmark()
        
      } else if(input$response_type == 'selgene') {
        glist  = gene_list_selected_df()
        output$ava_genes_table = renderDataTable({available_genelist()})  
      }
      glist
      
    })
    
    #Data for regularized regression
    reg_data = reactive({select(BRCA_metadata, glist_r()[,1], starts_with("hsa")) %>%
        mutate(response = select(., glist_r()[,1]) %>% rowMeans()) %>%
        select(response,glist_r()[,1],starts_with("hsa")) %>%
        na.omit()})
    
    
    #Outputs
    
    output$response_genes <- renderDataTable({
      DT::datatable(glist_r(), options = list(scrollX = T), fillContainer = TRUE) 
    })
    
    output$reg_df = renderDataTable({
      DT::datatable(select(reg_data(),response,glist_r()[,1]), options = list(scrollY = T), fillContainer = TRUE) 
    })
    
    output$selected_genes = renderText({ 
      paste("Selected Data&Prep technique is", input$response_type)})
    
    #Return Data
    
    return(reg_data)
    
  } )
}

regularization_ml_Server = function(id,regress_data) {
  moduleServer(id,function(input,output,session) {
    
    #Regularization
    response_var = reactive({as.matrix(select(regress_data(), response)) })
    predictor_var = reactive({ as.matrix(select(regress_data()[,colMeans(regress_data()) > 0 ], starts_with("hsa")))})
    
    foldid = reactive({sample(1:10,size = length(response_var()), replace = TRUE)})
    cvfit = reactive({cv.glmnet(predictor_var(), response_var(), foldid = foldid(), alpha = input$user_alpha)})
    
    output$lambda_error_plot = renderPlot({plot(cvfit(),xvar = "lambda", label = TRUE)})
    output$coef_lambda_plot = renderPlot({plot(cvfit()$glmnet.fit ,xvar = "lambda", label = TRUE)})
    
    #Coefficients
    coef_data = reactive({
      c = as.matrix(coef(cvfit(), s = input$lambda_for_coef))
      as.data.frame(c)})
    
    output$coef_data = renderDataTable({
      coef_data()
    })
    
    output$cf_summary = renderText({
      paste("", summary(coef_data()) )
    })
    
    #Lambda Values
    output$lambda_value_min = renderText({
      paste("lambda.min: ", cvfit()$lambda.min) })
    output$lambda_value_1se = renderText({
      paste("lambda.1se: ", cvfit()$lambda.1se)})
    
    
    
    # output$try_plot = renderPlot({
    #   km = kmeans(coef_data(), centers = 5, nstart = 100) #kmeans clustering.
    # 
    #   fviz_cluster(km, data = coef_data())
    #    
    # })
    
  })
}


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



# workmodule_Server = function(id,df) {
#   
#   moduleServer(id,function(input,output,session) {
# 
#     output$work_table = renderDataTable({df()})
#     
#   })
# }
