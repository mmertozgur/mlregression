library(dplyr)
library(edgeR)
library(tidyverse)
library(glmnet)
library(msigdbr)
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(plotly)
#Functions
msigdb_handler <-  function(category,subcategory) {
  hallmark_sets <- as.data.frame(msigdbr(species = "Homo sapiens",category = category, subcategory = subcategory))
  hallmark_sets <-hallmark_sets %>% select(c(gs_name,gene_symbol))
  hsets_names4freq <- as.data.frame(table(hallmark_sets %>% select(gs_name)))
  alist <- list("data" = hallmark_sets, "names" = hsets_names4freq)
  return(alist)
}
zero_adjuster = function(mdata, max_zero_percent, subgroup = "") {
  
  
  
  if (subgroup == "gene") {
    #mdata = mdata[-1]
    mdata <- select(mdata,-c("patient"))
    mdata <- mdata[,1:60000]
    zerofreq_list <- sapply(mdata, function(x){ length(which(x==0))/length(x) *100 })
    available <- names(zerofreq_list[zerofreq_list[] < max_zero_percent])
  } else if (subgroup == "miRNA") {
    mdata <- select(mdata,starts_with("hsa"))
    zerofreq_list <- sapply(mdata, function(x){ length(which(x==0))/length(x) *100 })
    available <- names(zerofreq_list[zerofreq_list[] < max_zero_percent])
    
  } else if (subgroup == "all") {
    
    mdata <- select(mdata, -c(starts_with("META"),starts_with("patient")))
    zerofreq_list <- sapply(mdata, function(x){ length(which(x==0))/length(x) *100 })
    available <- names(zerofreq_list[zerofreq_list[] < max_zero_percent])
  }
  
  
  return(as.data.frame(available))
}

#UI
dataprepInputControl_UI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Metadata"),
    selectizeInput(NS(id,"tcga_metadata"),"Select a TCGA metadata to analyze",choices = c("SKCM")),
    sliderInput(NS(id,"max_zero_percent"),"Set the maximum percentage of acceptable zero count (for manual selections)",
                min = 0, max = 100, value = 0, step = 1),
    h3("Response Variable"),
    radioButtons(NS(id,"response_prep_method"), "Data preparation pipeline: ",
                 c("Create response set from gene sets available on MSigDB" = "msigdb",
                   "Create response set manually" = "gene_list")
    ),
    conditionalPanel(condition = "input.response_prep_method == 'msigdb' ", ns = ns,
                     selectizeInput(NS(id,"msigdb_setnames_response"), "MSigDB Human Collections", choices = c("hallmark gene sets" = "H","ontology gene sets" = "C5" ,"oncogenic gene sets" = "C6",
                                                                                                               "immunologic gene sets" = "C7")),
                     conditionalPanel(condition = "input.msigdb_setnames_response == 'C5'", ns = ns, 
                                      radioButtons(NS(id,"go_term_response"),"Select go term" ,choices = c("GO:BP","GO:CC","GO:MF"))
                     ),
                     conditionalPanel(condition = "input.msigdb_setnames_response == 'C7'", ns = ns,
                                      radioButtons(NS(id,"immuno_response"), "Select subcategory", choices = c("IMMUNESIGDB","VAX"))
                     ),
                     selectizeInput(NS(id,"msigdb_gene_set_response"), "Select a gene set as a response set",
                                    choices = c(""))
    ),
    conditionalPanel(condition = "input.response_prep_method == 'gene_list' ", ns = ns,
                     radioButtons(NS(id,"obtain_response"), "Create response set by",c("Upload TXT file" = "txt_upload", "Enter gene" = "gene_enter" )),
                     conditionalPanel(condition = "input.obtain_response == 'gene_enter' ", ns = ns,
                                      selectizeInput(NS(id,"gene_list_response"),"Select gene(s)", choices = NULL, multiple = TRUE)
                     ),
                     conditionalPanel(condition = "input.obtain_response == 'txt_upload' ", ns = ns,
                                      fileInput(NS(id,"response_set_file"), "Upload File",
                                                accept =  c(".txt",".csv"))
                     )
    ),
    h3("Predictor Variables"),
    radioButtons(NS(id,"predictor_prep_method"), "Data preparation pipeline: ",
                 c("Create predictor set from gene sets available on MSigDB" = "msigdb",
                   "Create predictor set manually" = "gene_list")
    ),
    conditionalPanel(condition = "input.predictor_prep_method == 'msigdb' ", ns = ns,
                     selectizeInput(NS(id,"msigdb_setnames_predictor"), "MSigDB Human Collections", choices = c("hallmark gene sets" = "H","ontology gene sets" = "C5" ,
                                                                                                                "oncogenic gene sets" = "C6","immunologic gene sets" = "C7")),
                     conditionalPanel(condition = "input.msigdb_setnames_predictor == 'C5'", ns = ns,
                                      radioButtons(NS(id,"go_term_predictor"),"Select go term" ,choices = c("GO:BP","GO:CC","GO:MF"))
                     ),
                     conditionalPanel(condition = "input.msigdb_setnames_predictor == 'C7'", ns = ns,
                                      radioButtons(NS(id,"immuno_predictor"), "Select subcategory", choices = c("IMMUNESIGDB","VAX"))
                     ),
                     selectizeInput(NS(id,"msigdb_gene_set_predictor"), "Select a gene set as a response set",
                                    choices = c(""))
    ),
    conditionalPanel(condition = "input.predictor_prep_method == 'gene_list' ", ns = ns,
                     radioButtons(NS(id,"obtain_predictor"), "Create predictor set by",c("Upload TXT file" = "txt_upload", "Enter gene" = "gene_enter",
                                                                                         "Use all mRNAs available on the selected data(!)" = "allmRNA_aspredictor",
                                                                                         "Use all miRNAs available on the selected data" = "allmiRNA_aspredictor")),
                     conditionalPanel(condition = "input.obtain_predictor == 'gene_enter' ", ns = ns,
                                      selectizeInput(NS(id,"gene_list_predictor"),"Select gene(s)", choices = c(""), multiple = TRUE)
                     ),
                     conditionalPanel(condition = "input.obtain_predictor == 'txt_upload' ", ns = ns,
                                      fileInput(NS(id,"predictor_set_file"), "Upload File",
                                                accept =  c(".txt","csv"))
                     )
    )
    
    
  )
}
regression_sidecontrols <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Workflow"),
    radioButtons(NS(id,"regression_workflow"), "Select a workflow type to train the model",
                 c("Use 100% of the data to train a model" = "create_model",
                   "Split data into train and test sets and evaluate model accuracy" = "ctest_model")
    ),
    conditionalPanel(
      condition = "input.regression_workflow == 'ctest_model' ", ns = ns, 
      sliderInput(NS(id,"train_percentage"),
                  "Set the training set percentage (If 30% selected, 30% of the data will be used as training set)",
                  min = 10, max = 90, value = 10, step = 10),
      selectizeInput(NS(id,"lambda_for_prediction"), "Select lambda value for prediction: ", 
                     choices = c("lambda.min","lambda.1se"))
      
    ),
    h3("Regression: Ridge/ElasticNet/Lasso"),
    fluidRow(
      column(12,
             sliderInput(NS(id,"user_alpha"),"Enter an alpha parameter (shrinkage penalty term) to set the regression technique", min = 0, max = 1, value = 1, step = 0.2)
      )
      # column(4,
      #        actionButton(NS(id,"analyze_button"), "Analyze"))
    ),
    selectizeInput(NS(id,"lambda_for_coef"), "Select lambda value at which model coefficients are being displayed: ", 
                   choices = c("lambda.min","lambda.1se")),
  )
}
tcgaexplorer_ml_UI <- function(id) {
  ns <- NS(id)
  navbarPage(
    "TCGA Explorer ML",
    tabPanel(
      "Data&Prep",
      sidebarPanel(
        dataprepInputControl_UI("ml")
      ),
      mainPanel(
        fluidRow(
          column(6,h3("Response")),
          column(6,h3("Predictor"))
        ),
        fluidRow(
          
          column(6,dataTableOutput(NS(id,"response_set"))),
          
          column(6,dataTableOutput(NS(id,"predictor_set")))
        ),
        
        verbatimTextOutput(NS(id,"validation_message")),
        verbatimTextOutput(NS(id,"test"))
        
      )
    ),
    tabPanel(
      "Regression R/EN/L",
      sidebarPanel(
        regression_sidecontrols("ml")
      ),
      mainPanel(
        fluidRow(
          column(6,plotOutput(NS(id,"lambda_error_plot"))),
          column(6,plotlyOutput(NS(id,"coef_lambda_plot")))
        ),
        fluidRow(
          column(6,
                 dataTableOutput(NS(id,"coef_data"))),
          column(6,
                 textOutput(NS(id,"lambda_value_min")),
                 textOutput(NS(id,"lambda_value_1se")),
                 verbatimTextOutput(NS(id,"model_error")),
                 downloadButton(NS(id,"download_coef"),"Download Coefficients")
          )
        )
      )
    )
  )
}


#Server Modules
dataprepml_Serverz <- function(id) {
  moduleServer(id,function(input,output,session){
    
    metadata <- reactive({readRDS("~/Desktop/MLRegression/mlregression/skcm_merged_l_clean.rds")})
    
    
    #Response, obtain from msigdb
    response_msigdb_data <- reactive({
      if (input$msigdb_setnames_response == "C5") {
        msigdb_response <- msigdb_handler("C5",input$go_term_response)
        updateSelectizeInput(session, 'msigdb_gene_set_response', choices = msigdb_response$names$gs_name,server = TRUE)
        
        response_m <- as.data.frame(msigdb_response$data)
      } else if (input$msigdb_setnames_response == "C7") {
        msigdb_response <- msigdb_handler("C7",input$immuno_response)
        updateSelectizeInput(session, 'msigdb_gene_set_response', choices = msigdb_response$names$gs_name,server = TRUE)
        
        response_m <- as.data.frame(msigdb_response$data)
      } else {
        msigdb_response <- msigdb_handler(input$msigdb_setnames_response,"")
        updateSelectizeInput(session, 'msigdb_gene_set_response', choices = msigdb_response$names$gs_name,server = TRUE)
        
        response_m <- as.data.frame(msigdb_response$data)
      }
      response_m
    })
    
    response_msigdb_genes <- reactive({
      filter(response_msigdb_data(), gs_name == input$msigdb_gene_set_response) %>% select(gene_symbol)
    })
    
    
    
    #Predictor, obtain from msigdb.
    predictor_msigdb_data <- reactive({
      if (input$msigdb_setnames_predictor == "C5") {
        msigdb_predictor <- msigdb_handler("C5",input$go_term_predictor)
        updateSelectizeInput(session, 'msigdb_gene_set_predictor', choices = msigdb_predictor$names$gs_name,server = TRUE)
        
        predictor_m <- as.data.frame(msigdb_predictor$data)
      } else if (input$msigdb_setnames_predictor == "C7") {
        msigdb_predictor <- msigdb_handler("C7",input$immuno_predictor)
        updateSelectizeInput(session,'msigdb_gene_set_predictor', choices = msigdb_predictor$names$gs_name,server = TRUE)
        
        predictor_m <- as.data.frame(msigdb_predictor$data)
      } else {
        msigdb_predictor <- msigdb_handler(input$msigdb_setnames_predictor,"")
        updateSelectizeInput(session,'msigdb_gene_set_predictor', choices = msigdb_predictor$names$gs_name,server = TRUE)
        
        predictor_m <- as.data.frame(msigdb_predictor$data)
      }
      predictor_m
    })
    
    predictor_msigdb_genes <- reactive({
      filter(predictor_msigdb_data(), gs_name == input$msigdb_gene_set_predictor) %>% select(gene_symbol)
    })
    
    
    # create a list.
    
    available_genelist <- reactive({zero_adjuster(mdata = metadata(),max_zero_percent = input$max_zero_percent, subgroup = "all")})
    
    
    observeEvent(input$max_zero_percent, {
      
      genelist <- available_genelist()$available
      updateSelectizeInput(session,'gene_list_response', choices = genelist, server = TRUE)
      updateSelectizeInput(session, 'gene_list_predictor', choices = genelist, server = TRUE)
      
      output$test <- renderText({
        print(input$max_zero_percent)
      })
      
    })
    
    
    
    gene_list_selected_df_response <- reactive({as.data.frame(input$gene_list_response)})
    
    gene_list_selected_df_predictor <- reactive({as.data.frame(input$gene_list_predictor)})
    
    #
    response_det <- reactive({
      if (input$response_prep_method == "msigdb") {
        list_r = response_msigdb_genes()
      } else if (input$response_prep_method == "gene_list") {
        if (input$obtain_response == "gene_enter") {
          list_r = gene_list_selected_df_response()
        } else {
          file <- input$response_set_file
          ext <- tools::file_ext(file$datapath)
          req(file)
          validate(need(ext == "txt", "Please upload a txt file"))
          list_r = read.table(file$datapath, header = TRUE, sep = "", dec = ".")
          colnames(list_p) = "gene_symbol"
        }
        
      }
      list_r
    })
    
    predictor_det <- reactive({
      if (input$predictor_prep_method == "msigdb") {
        list_p <- predictor_msigdb_genes()
      } else if (input$predictor_prep_method == "gene_list") {
        if (input$obtain_predictor == "gene_enter") {
          list_p <- gene_list_selected_df_predictor()
        } else if (input$obtain_predictor == "allmiRNA_aspredictor") {
          mir_list <- colnames(select(metadata(),starts_with("hsa")))
          list_p <- as.data.frame(mir_list)
          
        } else {
          file <- input$predictor_set_file
          ext <- tools::file_ext(file$datapath)
          req(file)
          validate(need(ext == "txt", "Please upload a txt file"))
          list_p <- read.table(file$datapath, header = FALSE, sep = "", dec = ".")
          colnames(list_p) <- "gene_symbol"
        }
      }
      list_p
    })
    
    #
    output$response_set <- renderDataTable({response_det() }) #response determinants 
    output$predictor_set <- renderDataTable({predictor_det() }) #predictor determinants 
    
    output$validation_message <- renderText({
      genelist <- predictor_det()$gene_symbol
      is.exist <- genelist %in% colnames(metadata())
      x <- predictor_det()[which(is.exist == "FALSE"),]
      x
    })
    
    
    
    reg_data <- reactive({
      if (input$obtain_predictor == "allmiRNA_aspredictor") {
        select(metadata(),response_det()[,1], starts_with("hsa")) %>% 
          mutate(response = select(., response_det()[,1]) %>% rowMeans()) %>% 
          select(response,starts_with("hsa")) %>% 
          na.omit()
        
      } else {
        select(metadata(), response_det()[,1], predictor_det()[,1]) %>%
          mutate(response = select(., response_det()[,1]) %>% rowMeans()) %>%
          select(response, predictor_det()[,1]) %>% # ,response_det()[,1] --> predictor data
          na.omit()
      }
      
    })
    
    # output$test = renderText({
    #   print(colMeans(reg_data()[1:3]))
    # })
    
    
    return(reg_data)
  } )
}
regression_ml_Server <- function(id,regress_data) {
  moduleServer(id,function(input,output,session) {
    
    trows <- reactive({
      nrowd <- nrow(regress_data())
      set.seed(30)
      trowsx <- sample(1:nrowd, (input$train_percentage*nrowd)/100 ,replace = FALSE)
      trowsx
    })
    
    
    cvfit <- reactive({
      if (input$regression_workflow == "create_model") {
        response_var <- as.matrix(select(regress_data(), response)) 
        predictor_var <-  as.matrix(select(regress_data(), -c("response"))) 
        foldid <- sample(1:10,size = length(response_var), replace = TRUE)
        fitty <- cv.glmnet(predictor_var, response_var, foldid = foldid, alpha = input$user_alpha)
      } else if (input$regression_workflow == "ctest_model") {
        train_regress_data <- regress_data()[trows(),]
        response_var_train <- as.matrix(select(train_regress_data, response)) 
        predictor_var_train <-  as.matrix(select(train_regress_data, -c("response"))) 
        foldid <- sample(1:10,size = length(response_var_train), replace = TRUE)
        fitty <- cv.glmnet(predictor_var_train, response_var_train, foldid = foldid, alpha = input$user_alpha)
      }
      fitty
    })
    
    
    prediction_error <- reactive({
      if(input$regression_workflow == "ctest_model"){
        
        test_regress_data <- regress_data()[-trows(),]
        response_var_test <- as.matrix(select(test_regress_data, response)) 
        predictor_var_test <-  as.matrix(select(test_regress_data, -c("response"))) #[,colMeans(test_regress_data) > 0 ]
        prediction <- predict(cvfit(), s = input$lambda_for_prediction, newx = predictor_var_test)
        
        error <- mean((response_var_test - prediction)^2)
        
      }
      error
    })
    
    coef_data <- reactive({
      c <- as.matrix(coef(cvfit(), s = input$lambda_for_coef))
      c <-  c[2:length(rownames(c)), ]
      as.data.frame(c)})
    
    
    # observeEvent(input$analyze_button,
    #              renderPlot({plot(cvfit(),xvar = "lambda", label = TRUE)}),
    #              renderPlot({plot(cvfit()$glmnet.fit ,xvar = "lambda", label = TRUE)}),
    #              renderDataTable({
    #                coef_data()
    #              })
    #              )
    
    
    output$lambda_error_plot <- renderPlot({plot(cvfit(),xvar = "lambda", label = TRUE)})
    
    
    
    output$coef_lambda_plot <- renderPlotly({
      
      fit <- cvfit()$glmnet.fit
      l1se <- cvfit()$lambda.1se
      lmin <- cvfit()$lambda.min
      
      lam <- fit$lambda %>%
        as.data.frame() %>%
        mutate(penalty = fit$a0 %>% names()) %>%
        rename(lambda = ".")
      results <- fit$beta %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        gather(penalty, coefficients, -rowname) %>%
        left_join(lam)
      
      
      
      result_labels <- results %>%
        filter(lambda == min(lambda))
      g=ggplot() +
        geom_line(data = results, aes(lambda, coefficients, group = rowname, color = rowname), show.legend = FALSE) +
        scale_x_log10() +
        geom_text(data = result_labels, aes(lambda, coefficients, label = rowname, color = rowname), nudge_x = .15,nudge_y = 0.0035,
                  show.legend = FALSE)+
        geom_vline(xintercept = lmin) +
        #annotate("text", x=lmin, y=-2.5, label="lambda.min", angle=0, color="blue") +
        geom_vline(xintercept = l1se) 
      #annotate("text", x=l1se, y=2.5, label="lambda.1se", angle=0, color="red")
      
      g <- g + theme(legend.title = element_blank())
      g <- g + theme(legend.position='none')
      
      
      ggplotly(g)
      
      
      
      
      
    })
    
    output$coef_data <- renderDataTable({
      coef_data()
    })
    
    output$lambda_value_min <- renderText({
      paste("lambda.min: ", cvfit()$lambda.min) })
    
    output$lambda_value_1se <- renderText({
      paste("lambda.1se: ", cvfit()$lambda.1se)})
    
    output$model_error <- renderPrint({
      paste(" MSE Error: ", prediction_error() ) 
    })
    
    output$download_coef <- downloadHandler(
      filename <-  function() {
        paste(input$lambda_for_coef, ".csv", sep = "")
      },
      content <-  function(file) {
        write.csv(coef_data(), file, row.names = FALSE)
      }
    )
    
  })
}


#App-run
ui <- fluidPage(
  theme = shinytheme("flatly"),
  tcgaexplorer_ml_UI("ml")
)
server<- function(input,output,session) {
  dataprepml_Serverz("ml")
  
  df <- dataprepml_Serverz("ml")
  regression_ml_Server("ml",regress_data  = df)
  
}

shinyApp(ui, server)

