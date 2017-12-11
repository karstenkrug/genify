#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(pacman)
p_load(shiny)
p_load(shinydashboard)
p_load(magrittr)

source('src/gene_collapse.r')
source('src/gct-io.r')


# Define UI for application
ui <- fluidPage(
   
    tags$head(tags$link(rel="stylesheet", href="css/styles.css", type="text/css"),
                tags$script(src="getdata.js")),
    
    fluidRow(h3('Genify - generate gene-centric expression tables', align='left')),
    hr(),
    br(),
    uiOutput('dataDropper'),
    uiOutput('selectGeneCol'), 
   ## uiOutput('selectAnno'),
    uiOutput('selectExprsCol'),
    uiOutput('selectCombMeth'),
    fluidRow(uiOutput('download')),
    tableOutput('tables')
    
)

# Define server logic 
server <- function(input, output) {
  
  global.data <- reactiveValues()
  
  # ################################################
  # UI: data droper
  output$dataDropper <- renderUI({
    if(!is.null(global.data$data)) return()
  
    list(
    fluidRow(column(12, h4('Drop your dataset in the box below (txt, csv, tsv)'))),
    
    fluidRow(
        column(3),
        column(6, title = 'Drag and drop here', 
                  div(class="col-xl-12", id="drop-area", ondragover="dragOver(event)", 
                  ondrop="dropData(event)")
               ),
        column(3)
      
    )
    )
  })
  # ########################################
  # UI: select column containing gene names
  output$selectGeneCol <- renderUI({
    
      if(is.null(global.data$data)) return()
      if(!is.null(global.data$gene_col)) return()
      list(
        h4('Select column that contains gene symbols.'),
        selectizeInput('gene_select', label=NULL, choices = global.data$cn, selected=grep('geneSymbol', global.data$cn, value=T)  ),
        h4('Select second annotation column (GCT v1.2).'),
        selectizeInput('anno_select', label=NULL, choices = global.data$cn, selected=grep('entry_name', global.data$cn, value=T)  ),
        actionButton('gene_select_go', 'Next')
      )
    
  })
  
  # ########################################
  # UI: select expression data
  output$selectExprsCol <- renderUI({
    if(is.null(global.data$gene_col)) return()
    if(!is.null(global.data$exprs_col)) return()
      list(
        h4('Select columns containing expression data.'),
        actionButton('exprs_select_go', 'Next'),
        checkboxGroupInput('exprs_select', label=NULL, choices=global.data$cn, selected=grep(paste(126:131, collapse = '|'), global.data$cn, value=T))
        
      )
    
  })
  # ##########################################
  # UI: select method to combine expression data
  output$selectCombMeth <- renderUI({
    if(is.null(global.data$exprs_col)) return()
    if(!is.null(global.data$comb_meth)) return()
    list(
      h4('Select how to combine expression data.'),
      radioButtons('select_comb_meth', label=NULL, choices=c('Median', 'Mean'), selected='Median'),
      actionButton('GO', 'GO')
    )
  })
  

  # ##########################################
  # OBS: select expression columns
  observeEvent(input$exprs_select_go,{
    global.data$exprs_col <- input$exprs_select
  })
  
  # #################################################
  # OBS: select gene name column
  observeEvent(input$gene_select_go, {
    global.data$gene_col <- input$gene_select
    global.data$anno_col <- input$anno_select
  })
  # #################################################
  # OBS: upload table
  observeEvent(input$mydata, {
    
    SEPARATOR=c('\t', ',', ';')
    
    ## fn is a text stream, NOT a file
    fn=input$mydata[[1]]
    withProgress(message = 'Importing dataset:',
      tab.tmp <- read.csv(text = fn, stringsAsFactors = F, header = F)
    )
    # #############################
    # determine the separator
    if(ncol(tab.tmp) == 1){
         tab.sep=NULL
        # try to figure out the separator, DON'T USE THE HEADER FOR THAT
        # use the fourth row instead (should be data)
        for(s in SEPARATOR){
          tab <- strsplit(as.character(tab.tmp[4,1]), s) %>% unlist()
          #View(tab)
          if(length(tab) > 1){
            global.data$tabsep <- s
            break;
          }
        }
   
         withProgress(message = 'Importing dataset:',
                    tab <- read.csv(text = fn, stringsAsFactors = F, header = F, sep=global.data$tabsep)
         )    
      # turn list into a matrix
      #tab.list <- lapply(tab.tmp, function(x) unlist(strsplit(x, global.data$tabsep)) )
      #tab <- unlist(tab.list) %>% matrix(nrow=nrow(tab.tmp), byrow=T)
    
      # generate header
      cn <- make.names(tab[1, ])
      colnames(tab) <- cn 
      tab <- data.frame(tab[-c(1),], stringsAsFactors = F)
      global.data$cn <- cn
    
    } else {
      tab <- tab.tmp
      cn <- make.names(tab[1, ])
      colnames(tab) <- cn 
      tab <- data.frame(tab[-c(1),], stringsAsFactors = F)
      global.data$cn <- cn
    
    }

    global.data$data <- tab
  
    })
  
  # ######################################################
  # display table
  output$tables <- renderTable({
    
    if(is.null(global.data$gc)) return()
    
      head(global.data$gc)

    
  })
  # ##################################
  # UI: download
  output$download <- renderUI({
    if(is.null(global.data$gc)) return()
    list(
      h4(paste('Collapsed ', nrow(global.data$data), ' entries into ', nrow(global.data$gc), ' genes-centric rows using ', global.data$comb_meth, ' expression.', sep='')),
      br(),
      downloadButton('download_gct', label='Download'),
      br(),
      h4('Preview:')
      )
  })

  # ########################################
  # UI: select column containing gene names
  output$selectAnno <- renderUI({
    if(is.null(global.data$data)) return()
    if(!is.null(global.data$gene_col)) return()
    
    list(
   )
  
})

  
  
  # ##########################################
  # download table
  output$download_gct <- downloadHandler(
      filename =  paste( 'gene-centric_', nrow(global.data$gc), 'x', length(global.data$exprs_col),'.gct', sep=''),
      content = function(file){
          withProgress({
            setProgress(message = 'Export')
            
            write.gct(global.data$gct, file)
            
            #write.table(global.data$gc, file, sep='\t', quote=F, na='', row.names = F)
            })
      }
  )
  
  
  # ##########################################
  # OBS: GO
  observeEvent(input$GO,{
    
    global.data$comb_meth <- input$select_comb_meth
    
    tab <- global.data$data
    exprs <- tab[, global.data$exprs_col]
    exprs[exprs==''] <- NA
    anno <- tab[, setdiff(global.data$cn, global.data$exprs_col) ]
    
    withProgress(message = 'Rolling up to gene symbols',
      gc <- gene.collapse(exprs = exprs, anno = anno, id.col=global.data$gene_col, meth=tolower(global.data$comb_meth))
    )
    # gct
    mat <- data.matrix( gc[, global.data$exprs_col] )
    rdesc <- data.frame( Description=gc[, global.data$anno_col], stringsAsFactors = F)
    rownames(rdesc) <- gc[, global.data$gene_col]
    #View(rdesc)
    gct <- new('GCT')
    gct@mat <- mat
    gct@rid <- rownames(rdesc)
    gct@cid <- colnames(mat)
    gct@rdesc <- rdesc
    gct@version <- '#1.2'
    
    
    ##
    gc <- data.frame(id=gct@rid, Description=rdesc, mat, stringsAsFactors = F)
    global.data$gc <- gc
    global.data$gct <- gct
    
    
  })
  
  
  }

# Run the application 
shinyApp(ui = ui, server = server)

