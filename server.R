# Brandon Ng
# 
# Server for registering files into a database
# 

source("volcano_helper.R")


# Defines server logic  ----
server <- function(input, output, session) 
{
  # Reads database file for use  ----
  renderDB <- function()
  {
    db <- read.table("db.csv",
                     header=TRUE,
                     sep=",",
                     colClasses="character")
    db <- db[-1] %>% setDT()
  }
  
  # Aggregates all form inputs and adds them to the database  ----
  submitForm <- function()
  {
    db <- renderDB()
    
    # Aggregates all inputs  ----
    # First column records submission time and date, thus i+1
    entry <- as.data.table(matrix("0",nrow=1,ncol=(length(fields)+1)))
    
    entry[,1] <- as.character(Sys.time())
    for (i in 1:length(fields))
    {
      if (fields[i] == "dataFile")
        entry[,i+1] <- input$dataFile$name
      
      else if (fields[i] == "date")
        entry[,i+1] <- as.character(input$date)
      
      else
        entry[,i+1] <- input[[fields[i]]]
    }
    
    # Creates new row, containing all new entries, with current time as row name  ----
    db <- bind_rows(db,entry)
    
    # Saves updated database  ----
    write.csv(db,"db.csv")
    
    reset("form")
    sendSweetAlert(session=session,
                   title="Success!",
                   text="Successful Submition!",
                   type="success")
  }
  
  # Submits inputs to database when submit button pressed  ----
  observeEvent(input$submit, submitForm())
  
  
  # Resets form's inputs when reset button pressed  ----
  observeEvent(input$reset, reset("form"))
  
  
  # Allows submit button to be pressed when all empty fields have a value  ---- 
  observe({
    # Determines if all empty fields have a value  ----
    fieldsFilled <- vapply(empty_fields, 
                           function(x) 
                           {!is.null(input[[x]]) && input[[x]] != ""}, 
                           FALSE)
    fieldsFilled <- all(fieldsFilled)
    
    # Toggles submit button allowing it to be pressed based on condition  ----
    toggleState("submit", condition=fieldsFilled)
  })
  
  
  # Renders database  ----
  output$database <- renderDT({
    datatable(
      renderDB()
    )
  })
  
  # Searches database based on user search inputs  ----
  searchDB <- reactive({
    
    db <- renderDB()
    params <- input$search_input
    
    index <- list()
    for (i in 1:length(params))
    {
      for (j in 1:ncol(db))
      {
        if (str_detect(db[[j]],params[[i]]))
        {
          db <- filter(db,str_detect(db[[j]],params[[i]]))
        }
      }
    }
    db
  })
  
  
  output$search_output <- renderDT({
    
    if (is.null(input$search_input))
      return(NULL)
    else
      return(datatable(
        searchDB()
      ))
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # # VOLCANO SERVER  ----
  # #
  # #
  # #
  # #
  # #
  # #
  # Determines getPlex size of experiment ----
  getPlex <- reactive({
    as.numeric(input$num_groups) * as.numeric(input$num_replicates)
  })
  
  
  # Gets sample names and labels user provided  ----
  getSamples <- reactive({
    as.data.frame(read_excel(input$tmt$datapath))
  })
  
  
  # General housekeeping tasks to prepare dataset for analysis  ----
  #   - Replaces missing values with 0
  #   - Removes proteins with less PSMs than specified by user
  #   - Removes all kerotin proteins
  #   - Labels abundance column with corresponding sample name
  #   - Displays only PSMs, # Unique peptides, Gene name, & Protein abundances
  tidyData <- reactive({
    # Replaces missing values in dataset with 0  ----
    dat <- input$volc_data$datapath %>% read_excel() %>% 
      mutate_all(funs(replace(.,is.na(.),0)))
    samples <- getSamples()
    
    # Removes all proteins with less or equal to PSMs than specified by user   ----
    names(dat)[colnames(dat) == "# Unique Peptides"] <- "UniquePeptides"
    names(dat)[colnames(dat) == "# PSMs"] <- "PSMs"
    dat <- subset(dat, UniquePeptides > input$uniq_peptide_cutoff)
    
    # Removes all kerotin proteins found in dataset  ----
    description <- dat["Description"]
    keratin <- 0
    for (i in 1:nrow(dat))
    {
      if(grepl("keratin", description[i,1],ignore.case=TRUE))
      {
        keratin <- c(i,keratin)
      }
    }
    dat <- dat[-keratin,]
    
    # Renames column names in dataset with corresponding sample name  ----
    for (i in 1:ncol(dat))
    {
      for (j in 1:nrow(samples))
      {
        if (grepl(samples[j,2], names(dat)[i], ignore.case=TRUE))
        {
          names(dat)[i] <- samples[j,1]
        }
      }
    }
    
    # Selects protein abundances from dataset   ----
    selected <- as.data.frame(matrix(0,nrow(dat),getPlex()))
    for (i in 1:nrow(samples))
    {
      names(selected)[i] <- samples[i,1]
      selected[i] <- dat[samples[i,1]]
    }
    
    # Extracts gene name from protein description  ----
    description <- dat["Description"]
    geneName <- as.data.frame(matrix(0,nrow(description),1))
    names(geneName)[1] <- "Gene Name"
    
    for (i in 1:nrow(geneName))
    {
      geneName[i,1] <- getGeneName(description[i,1])
    }
    
    # Creates datatable with gene name, PSMs, & # unique peptides alongside abundances ----
    cbind(geneName, dat["PSMs"], dat["UniquePeptides"], selected)
  })
  
  
  # Normalizes data about median  ----
  data_normalized <- reactive({
    # Gene name located in column 1, PSMs located in column 2,
    # # unique peptides located in column 3, thus -1,-2,-3
    dat <- tidyData()
    abundances <- dat[-c(1,2,3)]
    
    sums <- colSums(abundances)
    median <- median(sums)
    percent_median <- median / sums
    data_normalized <- as.data.frame(matrix(0,nrow(abundances),getPlex()))
    
    for(i in 1:ncol(abundances))
    {
      data_normalized[i] <- abundances[,i] * percent_median[i]
      names(data_normalized)[i] <- names(abundances)[i]
    }
    data_normalized <- roundValues(data_normalized,3)
    cbind(dat["Gene Name"], dat["PSMs"], dat["UniquePeptides"], data_normalized)
  })
  
  
  # Performs t-test among all groups  ----
  data_pval <- reactive({
    toAllGroups(data_normalized()[-c(1,2,3)], input$num_replicates, "pVal") %>%
      removeNan() %>% roundValues(.,5)
  })
  
  
  # Calculates negative log10 of p-values  ----
  data_log_pval <- reactive({
    -log10(data_pval())
  })
  
  
  # Calculates fold-change among all groups  ----
  data_fc <- reactive({
    toAllGroups(data_normalized()[-c(1,2,3)], input$num_replicates, "foldChange") %>%
      removeNan() %>% roundValues(.,5)
  })
  
  
  # Calculates log2 of fold-change  ----
  data_log_fc <- reactive({
    log2(data_fc())
  })
  
  
  # Creates volcano plot  ----
  volcano <- reactive({
    dat <- cbind(data_log_fc(),data_log_pval())
    og_data <- tidyData()
    pval <- -log10(input$pval_threshold)
    fc <- log2(input$fc_threshold)
    
    # Determines number of group comparisons
    num_comparisons <- ncol(dat) / 2
    
    # Attaches protein descriptions to the merged dataset that
    # contains log transformed fold-change and pval
    dat <- cbind(og_data["Gene Name"], dat)
    
    
    # Creates individual volcano plot  ----
    createPlot <- function(index)
    {
      # plot_ly: Plots fold-change on x-axis and pval on y-axis
      # text: Specifies that protein description appear when a
      #       datapoint is hovered over
      # name: Specifies the name that appears on the legend for
      #       that specific dataset
      plot_ly(x=dat[[index+1]],
              y=dat[[index+num_comparisons+1]],
              text=dat[[1]],
              name=substr(names(dat)[index+1],12,nchar(names(dat)[index+1]))) %>%
        
        # Specifies that each datapoint is a hollow circle
        add_markers(symbol=I(1)) %>%
        
        # Creates a horiztonal line layered over the origianl
        # plot, representing the p-val threshold
        add_lines(y=pval,
                  line=list(color="pink"),
                  showlegend=FALSE,
                  hoverinfo= "text",
                  text="P-val Threshold") %>%
        
        # Creates a vertical line layered over the original
        # plot, representing positive fold-change threshold
        add_lines(x=fc,
                  line=list(color="pink"),
                  showlegend=FALSE,
                  hoverinfo="text",
                  text="Fold-change Threshold") %>%
        
        # Creates a vertical line layered over the original
        # plot, representing negative fold-change threshold
        add_lines(x=-fc,
                  line=list(color="pink"),
                  showlegend=FALSE,
                  hoverinfo="text",
                  text="Fold-change Threshold") %>%
        
        # Specifies x and y axis titles
        layout(xaxis=list(title="Log2 fold-change"),
               yaxis=list(title="-Log10 p-val"))
    }
    
    # Creates a list of plots for group comparisons
    plotList <- function(nplot)
    {
      lapply(seq_len(nplot), createPlot)
    }
    
    # Creates single plot containing all group comparison plots
    subplot(plotList(num_comparisons),
            shareX=TRUE,
            shareY=TRUE,
            titleX=TRUE,
            titleY=TRUE)
  })
  
  
  # Determines names of all comparison groups  ----
  getComparisonNames <- reactive({
    names <- colnames(data_fc())
    # Removes "fold-change" originally found in column names
    names <- lapply(names, function(x) {substr(x,12,nchar(x))})
    names
  })
  
  
  # Produces prompt to select pval & fold-change comparison group to view  ----
  #   -- Includes "all" option --
  comparison_group_prompt <- reactive({
    names <- cbind(getComparisonNames(),"all")
    selectizeInput("compare_group","Comparison Group:",
                   choices=names,selected=names[1],multiple=FALSE)
  })
  
  
  # Produces prompt to select significant protein comparison group to view  ----
  #   -- Doesn't include "all" option --
  sig_comparison_prompt <- reactive({
    names <- getComparisonNames()
    selectizeInput("sig_compare_group","Comparison Group:",
                   choices=names,selected=names[1],multiple=FALSE)
  })
  
  
  # Produces data table of pval & fc of selected comparison group  ----
  comparison_group_data <- reactive({
    
    # Selected comparison group
    comp_grp <- input$compare_group
    
    dat <- tidyData()
    fc_data <- data_fc()
    pval_data <- data_pval()
    
    if (input$compare_group == "all")
    {
      # Merge function puts corresponding pvals next to corresponding fold-changes
      merged <- merge(fc_data,pval_data) 
      comp_group_data <- cbind(dat["Gene Name"], dat["PSMs"], dat["UniquePeptides"], merged)
    }
    else
    {
      # Determines which comparison group user selected
      for (i in 1:ncol(fc_data))
      {
        if (grepl(comp_grp, names(fc_data)[i], ignore.case=TRUE))
        {
          merged <- cbind(fc_data[i],pval_data[i]) 
          comp_group_data <- cbind(dat["Gene Name"], dat["PSMs"], dat["UniquePeptides"], merged)
        }
      }
    }
    comp_group_data
  })
  
  
  # Produces data table of significant proteins  ----
  sig_protein_data <- reactive({
    
    # Selected comparison group
    comp_grp <- input$sig_compare_group
    
    dat <- tidyData()
    fc_data <- data_fc()
    pval_data <- data_pval()
    
    # Determines which comparison group user selected
    for (i in 1:ncol(fc_data))
    {
      if (grepl(comp_grp, names(fc_data)[i], ignore.case=TRUE))
      {
        sig_protein <- getSigProt(cbind(dat["Gene Name"], 
                                        dat["PSMs"],
                                        dat["UniquePeptides"],
                                        fc_data[i],
                                        pval_data[i]))
      }
    }
    sig_protein
  })
  
  
  # Determines proteins beyond pval and fold-change thresholds  ----
  getSigProt <- function(dat)
  {
    fc <- input$fc_threshold
    pval <- input$pval_threshold
    
    # Vector keeping track of rows to include
    keep <- 0
    
    # Gene name is at position 1, PSMs is at position 2, Unique peptides is at position 3,
    # Fold-change is at position 4 & Pval is at position 5
    for (i in 1:nrow(dat))
    {
      if (dat[i,5] < pval & (dat[i,4] > fc | dat[i,4] < 1/fc) & (dat[i,4] != 0))
      {
        keep <- c(keep,i)
      }
    }
    return(dat[keep,])
  }
  
  
  # Misc.: Determines percent median of all samples  ----
  percent_median <- reactive({
    dat <- tidyData()[-c(1,2,3)]
    names <- colnames(dat)
    sums <- colSums(dat)
    median <- median(sums)
    percent_median <- median / sums
    as.data.frame(cbind(names,percent_median))
  })
  
  
  # Renders all outputs  ----
  
  # Delays rendering until datatset and TMT map has be uploaded  ----
  coverUp <- function(x)
  {
    if (is.null(input$volc_data) || is.null(input$tmt))
      return(NULL)
    else
      return(x)
    
  }
  
  
  # Renders dataset with only protein abundances
  output$tidy_data <- renderDT({
    coverUp(
      datatable(
        tidyData(),
        rownames=FALSE
      )
    )
  })
  
  # Renders dataset of protein abundances after being normalized about median
  output$data_normalized <- renderDT({
    coverUp(
      datatable(
        data_normalized(),
        rownames=FALSE
      )
    )
  })
  
  # Renders prompt to select comparison group
  output$comparison_group <- renderUI({
    coverUp(
      comparison_group_prompt()
    )
  })
  
  # Renders table of fold-changes & pvals for selected comparison group
  output$comparison_group_data <- renderDT({
    coverUp(
      datatable(
        comparison_group_data(),
        rownames=FALSE
      )
    )
  })
  
  # Renders button to download fc & pval data  ----
  output$fc_pval_button <- renderUI({
    coverUp(
      downloadButton("fc_pval_dwnld", label=tags$b("Download"))
    )
  })
  
  # Downloads fc & pval data once download button has been pressed  ----
  output$fc_pval_dwnld <- downloadHandler(
    filename = function() {
      paste("data-",Sys.Date(),".csv",sep="")
    },
    contentType="csv",
    content = function(file)
    {
      write.csv(comparison_group_data(),file)
    }
  )
  
  # Renders volcano plots for all combinations of groups
  output$volcano <- renderPlotly({
    coverUp(
      volcano()
    )
  })
  
  # Renders prompt to select comparison group for significant protein
  output$sig_comparison_group <- renderUI({
    coverUp(
      sig_comparison_prompt()
    )
  })
  
  # Renders table of significant proteins
  output$sig_protein_data <- renderDT({
    coverUp(
      datatable(
        sig_protein_data(),
        rownames=FALSE
      )
    )
  })
  
  # Renders button to download significant protein data  ----
  output$sig_prot_button <- renderUI({
    coverUp(
      downloadButton("sig_prot_dwnld", label=tags$b("Download"))
    )
  })
  
  # Downloads significant protein data once download button has been pressed  ----
  output$sig_prot_dwnld <- downloadHandler(
    filename = function() {
      paste("data-",Sys.Date(),".csv",sep="")
    },
    contentType="csv",
    content = function(file)
    {
      write.csv(sig_protein_data(),file)
    }
  )
  
  # Renders table of pecent medians of all samples
  output$percent_median <- renderDT({
    coverUp(
      datatable(
        percent_median(),
        rownames=FALSE
      )
    )
  })
  
  # Renders prompt to select comparison group for p-value histogram
  output$sig_comparison_group2 <- renderUI({
    coverUp(
      comparison_prompt2()
    )
  })
  
  
  # Renders histogram of p-values
  output$pvalHisto <- renderPlot({
    coverUp(
      createPvalHisto()
    )
  })
  
  
  # Renders button to download p-value histogram  ----
  output$pvalHisto_button <- renderUI({
    coverUp(
      downloadButton("pvalHisto_dwnld", label=tags$b("Download"))
    )
  })
  
  # Produces prompt to select which comparison group to view w/o "all" option ----
  comparison_prompt2 <- reactive({
    names <- getComparisonNames()
    selectizeInput("compare_group2","Comparison Group:",
                   choices=names,selected=names[1],multiple=FALSE)
  })
  
  
  # Produces histogram of p-value of selected comparison group  ----
  createPvalHisto <- reactive({
    
    # Selected comparison group
    comp_grp <- input$compare_group2
    
    # Data set containing pvals for all comparison gorups
    pval_data <- data_pval()
    
    # Determines which comparison group user selected and selects corresponding pval
    for (i in 1:ncol(pval_data))
    {
      if (grepl(comp_grp, names(pval_data)[i], ignore.case=TRUE))
      {
        index <- i
      }
    }
    pval_data <- pval_data[index] %>% unlist() %>% as.numeric()
    hist(pval_data,
         xlab="p-value",
         main=paste(input$compare_group2," p-value Histogram"))
  })
  
}