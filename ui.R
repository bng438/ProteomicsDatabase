# Brandon Ng
# 
# UI for registering files into a database
# 


# Defines UI  -----
ui <- dashboardPage(
  
  # Header  ----
  dashboardHeader(title=HTML("<i class='fab fa-fort-awesome'> Proteomics</i>")),
  
  
  # Sidebar  ----
  dashboardSidebar(
    sidebarMenu(
      # Home tab in sidebar
      menuItem(tags$b("Home"), tabName="home",icon=icon("home")),
      
      # Volcano plotter tab in sidebar
      menuItem(tags$b("Volcano"), tabName="volcano", icon=icon("cookie-bite"))
    )
  ),
  
  
  # Body  ----
  dashboardBody(
    useShinyjs(),
    
    # Contains the bodies for all tabs  ----
    tabItems(
      
      # Body for home tab  ----
      tabItem(
        tabName="home",
        
        # Row containing boxes of values  ----
        fluidRow(
          valueBox(0,"Searches",icon=icon("coffee")),
          valueBox(1,"Widgets", icon=icon("laptop-code"), color="purple"),
          valueBox(0,"Files archived",icon=icon("archive"), color="yellow")
        ),
        
        # Row containing form to upload file  ----
        fluidRow(
          box(title=HTML("<i class='fas fa-clipboard-list'> Upload</i>"),
              background="olive",
              solidHeader=TRUE,
              collapsible=TRUE,
              width="12",
              
              # Contains all ui inputs of form ----
              div(
                id="form",
                # Input for uploading data file  ----
                fileInput("dataFile", label="Data File"),
                
                fluidRow(
                  # Input to specify date experiment was done  ----
                  column(4,dateInput("date",
                                     label="Date",
                                     format="yyyy-mm-dd")),
                  
                  # Input to specify scientist who conducted experiment  ----
                  column(4,selectizeInput("scientist",
                                          label="Scientist",
                                          choices=scientists,
                                          options=list(create=TRUE,
                                                       placeholder="Select a scientist",
                                                       onInitialize = I('function() { this.setValue(""); }')))),
                  
                  # Input to specify protein of interest in experiment  ----
                  column(4,selectizeInput("protein",
                                          label="Protein",
                                          choices=proteins,
                                          options=list(create=TRUE,
                                                       placeholder="Select a protein",
                                                       onInitialize = I('function() { this.setValue(""); }'))))
                ),
                
                fluidRow(
                  # Input to specify cell line used in experiment  ----
                  column(4,selectizeInput("cell",
                                          label="Cell Line",
                                          choices=cells,
                                          options=list(create=TRUE,
                                                       placeholder="Select a cell line",
                                                       onInitialize = I('function() { this.setValue(""); }')))),
                  
                  # Input to specify drug used in experiment  ----
                  column(4,selectizeInput("drug",
                                          label="Drug",
                                          choices=drugs,
                                          options=list(create=TRUE,
                                                       placeholder="Select a drug",
                                                       onInitialize = I('function() { this.setValue(""); }')))),
                  
                  # Input to specify which ELN used for experiment  ----
                  column(4, textInput("eln", label="ELN"))
                ),
                
                # Input for brief description of experiment & any comments  ----
                textInput("description", label="Description"),
                
                # Submition button  ----
                actionButton("submit", label=tags$b("Submit"), icon=icon("file-upload")),
                
                # Form reset button  ----
                actionButton("reset", label=tags$b("Reset"), icon=icon("redo-alt"))
              )
          )
        ),
        
        # Row containing search  ----
        fluidRow(
          box(title=HTML("<i class='fas fa-search'> Search</i>"),
              status="info",
              solidHeader=TRUE,
              collapsible=TRUE,
              width="12",
              
              # Contains all ui inputs of search ----
              div(
                id="search",
                DTOutput("db")
              )
          )
        )
      ),
      
      
      # Body for volcano tab  ----
      tabItem(
        tabName="volcano",
        
        # Creates sidebar layout with input and output definitions
        sidebarLayout(
          
          # Sidebar panel for intputs ----
          sidebarPanel(
            
            # Stores excel data file
            fileInput("volc_data", h4("Data File")),
            
            # Stores tmt label map
            fileInput("tmt", h4("TMT Map")),
            
            # Stores number of groups in experiment
            numericInput("num_groups", h5("Number of Groups"),value=3,min=1),
            
            # Stores number of replicates per group
            numericInput("num_replicates", h5("Replicates per Group"),value=3,min=1,max=9),
            
            # Stores p-val threshold
            numericInput("pval_threshold", h5("Pval Threshold"), value=.05, min=0),
            
            # Stores fold-change threshold
            numericInput("fc_threshold", h5("Fold-change Threshold"), value=1.414, min=0),
            
            # Text box detailing limitations to script
            helpText("Data file should not contain abundance ratios")
          ),
          
          # Main panel for displaying outputs  ----
          mainPanel(
            
            # Sets up different tabs
            tabsetPanel(type="tabs",
                        
                        # Displays raw protein abundance values
                        tabPanel("Raw Data",
                                 DTOutput("data_tidy")
                                 ),
                        
                        # Displays normalized protein abundance values
                        tabPanel("Normalized",
                                 dataTableOutput("data_normalized")
                                 ),
                        
                        # Displays fold-change and pvals for selected comparison group
                        tabPanel("Fold-change & Pval",
                                 uiOutput("comparison_group"),
                                 dataTableOutput("comparison_group_data"),
                                 uiOutput("fc_pval_button")
                                 ),
                        
                        # Displays volcano plot for all combination of groups
                        tabPanel("Volcano",
                                 plotlyOutput("volcano")
                                 ),
                        
                        # Displays significant proteins
                        tabPanel("Significant Proteins",
                                 uiOutput("sig_comparison_group"),
                                 dataTableOutput("sig_protein_data"),
                                 uiOutput("sig_prot_button")
                                 ),
                        
                        # Displays miscellaneous statistical values
                        tabPanel("Misc",
                                 dataTableOutput("percent_median")
                                 )
            )
          )
        )
      )
    )
  )
)