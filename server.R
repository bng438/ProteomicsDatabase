# Brandon Ng
# 
# Server for registering files into a database
# 

source("volcano_helper.R")


# Defines server logic  ----
server <- function(input, output, session) 
{
  # # Determines time in readable format
  # humanTime <- function() 
  # {
  #   return(format(Sys.time(), "%Y%m%d-%H%M%OS"))
  # }
  
  
  # Aggregates all form inputs and adds them to the database  ----
  submitForm <- function()
  {
    db <- read.table("db.csv",
                     header=TRUE,
                     sep=",",
                     colClasses="character"
    )
    db <- db[-1] %>% setDT()
    
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
    
    # Creates new column, containing all new entries, with current time as name  ----
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
  
  
  dbRender <- reactive({
    db <- read.csv("db.csv")
    db <- db[-1] %>% setDT()
    db
  })
  
  # Renders database
  output$db <- renderDT({
    datatable(
      dbRender(),
      extensions="Responsive"
    )
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # VOLCANO SERVER  ----
  # 
  # 
  # 
  # 
  # 
  # 
  # Determines plex size of experiment ---- 
  plex <- reactive({
    as.numeric(input$num_groups) * as.numeric(input$num_replicates)
  })
  
  
  # Gets sample names and labels user provided  ----
  getSamples <- reactive({
    samples <- input$tmt$datapath %>% read_excel()
    as.data.frame(samples)
  })
  
  
  # Replaces missing values in dataset with 0  ----
  data_rm0 <- reactive({
    input$volc_data$datapath %>%
      read_excel() %>%
      mutate_all(funs(replace(.,is.na(.),0)))
  })
  
  
  # Renames column names with corresponding sample name  ----
  data_rename <- reactive({
    samples <- getSamples()
    dat <- data_rm0()
    
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
    dat
  })
  
  
  # Removes all columns except protein abundances  ----
  data_selected <- reactive({
    samples <- getSamples()
    dat <- data_rename()
    selected <- as.data.frame(matrix(0,nrow(dat),plex()))
    
    for (i in 1:nrow(samples))
    {
      names(selected)[i] <- samples[i,1]
      selected[i] <- dat[samples[i,1]]
    }
    selected
  })
  
  
  # Normalizes data about median  ----
  data_normalized <- reactive({
    dat <- data_selected()
    sums <- colSums(dat)
    median <- median(sums)
    percent_median <- median / sums
    data_normalized <- as.data.frame(matrix(0,nrow(dat),plex()))
    
    for(i in 1:ncol(dat))
    {
      data_normalized[i] <- dat[,i] * percent_median[i]
      names(data_normalized)[i] <- names(dat)[i]
    }
    data_normalized
  })
  
  
  # Performs t-test among all groups  ----
  data_pval <- reactive({
    toAllGroups(data_normalized(), input$num_replicates, "pVal")
  })
  
  
  # Calculates negative log10 of p-values  ----
  data_log_pval <- reactive({ 
    -log10(data_pval())
  })
  
  
  # Calculates fold-change among all groups  ----
  data_fc <- reactive({
    toAllGroups(data_normalized(), input$num_replicates, "foldChange")
  })
  
  
  # Calculates log2 of fold-change  ----
  data_log_fc <- reactive({
    log2(data_fc())
  })
  
  
  # Creates volcano plot  ----
  volcano <- reactive({
    dat <- cbind(data_log_fc(),data_log_pval()) %>% removeNan()
    og_data <- data_rm0()
    pval <- -log10(input$pval_threshold)
    fc <- log2(input$fc_threshold)
    
    # Determines number of group comparisons
    num_comparisons <- ncol(dat) / 2
    
    # Attaches protein descriptions to the merged dataset
    # which contains fold-change and pval
    dat <- cbind(og_data["Description"], dat)
    
    
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
        layout(xaxis=list(title="Log2 Fold-change"),
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
  
  
  # Produces prompt to select which comparison group to view includes "all" option  ----
  comparison_group_prompt <- reactive({
    names <- cbind(getComparisonNames(),"all")
    selectInput("compare_group","Comparison Group:",
                choices=names)
  })
  
  
  # Produces prompt to select which comparison group to view w/o "all" option ----
  sig_comparison_prompt <- reactive({
    names <- getComparisonNames()
    selectInput("sig_compare_group","Comparison Group:",
                choices=names)
  })
  
  
  # Produces data table of pval & fc of selected comparison group  ----
  comparison_group_data <- reactive({
    
    # Selected comparison group
    comp_grp <- input$compare_group
    
    # Data set containing protein description
    raw_data <- data_rm0()
    
    fc_data <- data_fc()
    pval_data <- data_pval()
    
    if (input$compare_group == "all")
    {
      # Merge function puts corresponding pvals next to corresponding fold-changes
      merged <- merge(fc_data,pval_data) %>% removeNan() %>% roundValues(.,4)
      comp_group_data <- cbind(raw_data["Description"],merged)
    }
    else
    {
      # Determines which comparison group user selected
      for (i in 1:ncol(fc_data))
      {
        if (grepl(comp_grp, names(fc_data)[i], ignore.case=TRUE))
        {
          merged <- cbind(fc_data[i],pval_data[i]) %>% removeNan() %>% roundValues(.,4)
          comp_group_data <- cbind(raw_data["Description"],merged)
        }
      }
    }
    comp_group_data
  })
  
  
  # Produces data table of significant proteins  ----
  sig_protein_data <- reactive({
    
    # Selected comparison group
    comp_grp <- input$sig_compare_group
    
    # Data set containing protein description
    raw_data <- data_rm0()
    
    fc_data <- data_fc() %>% removeNan() %>% roundValues(.,4)
    pval_data <- data_pval() %>% removeNan() %>% roundValues(.,4)
    
    # Determines which comparison group user selected
    for (i in 1:ncol(fc_data))
    {
      if (grepl(comp_grp, names(fc_data)[i], ignore.case=TRUE))
      {
        sig_protein <- getSigProt(cbind(raw_data["Description"],fc_data[i],pval_data[i]))
      }
    }
    # write.csv(sig_protein,file=paste(comp_grp,".csv"))
    sig_protein
  })
  
  
  # Determines proteins beyond pval and fold-change thresholds  ----
  getSigProt <- function(dat)
  {
    fc <- input$fc_threshold
    pval <- input$pval_threshold
    keep <- 0
    
    for (i in 1:nrow(dat))
    {
      if (dat[i,3] < pval & (dat[i,2] > fc | dat[i,2] < 1/fc) & (dat[i,2] != 0))
      {
        keep <- c(keep,i)
      }
    }
    return(dat[keep,])
  }
  
  
  # Misc.: Determines percent median of all samples  ----
  percent_median <- reactive({
    dat <- data_selected()
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
  output$data_tidy <- renderDT({
    coverUp(
      datatable(
        data_selected(),
        rownames=FALSE,
        options=list(columnDefs = list(list(className="dt-ceter")))
      )
    )
  })
  
  # Renders dataset of protein abundances after being normalized about median
  output$data_normalized <- renderDataTable({
    coverUp(
      roundValues(data_normalized(),3)
    )
  })
  
  # Renders prompt to select comparison group
  output$comparison_group <- renderUI({
    coverUp(
      comparison_group_prompt()
    )
  })
  
  # Renders table of fold-changes & pvals for selected comparison group
  output$comparison_group_data <- renderDataTable({
    coverUp(
      comparison_group_data()
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
  output$sig_protein_data <- renderDataTable({
    coverUp(
      sig_protein_data()
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
  output$percent_median <- renderDataTable({
    coverUp(
      percent_median()
    )
  })
  
}