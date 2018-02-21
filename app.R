### mass spec app ## Wojciech Rymaszewski 2017_03_21
library(shiny)
library(XLConnect)
library(gplots)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


ui = fluidPage (
  includeCSS("styles.css"),
  headerPanel("Peptide lists comparison"),
  sidebarLayout(
    sidebarPanel(
      fileInput('uploadFile', "Upload file"),
      textOutput('variants'),
      tags$hr(),
      textInput('chosen_variants', 'Choose variants (separated by ";")'),
      textInput('zero', 'Choose mock'),
      textInput('excluded', 'Excluded expressions (separated by ";")'),
      textInput('included', 'Included expressions (separated by ";")'),
      actionButton('submit', "Submit")
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Venn diagram", plotOutput('venn', height='500px', width = '500px')),
        tabPanel("All peptides",
                 downloadButton('download_all', 'Download CSV'),
                 tableOutput('tab_all')),
        tabPanel("Venn intersections",
                 downloadButton('download_intr', 'Download CSV'),
                 tableOutput('tab_intersections')),
        tabPanel("Common peptides",
                 downloadButton('download_dups', 'Download CSV'),
                 tableOutput('dups'))
      )
    )
  )
)

server = function(input,output) {
 wb = reactive({
     if (is.null(input$uploadFile)){
       return (NULL)
     }
     loadWorkbook(input$uploadFile$datapath)
     })
 sheet.list = reactive({
     if (is.null(wb())) {
       return (NULL)
     }
     getSheets(wb())
     })
 output$variants = renderText({
     if(is.null(sheet.list())){
       return (NULL)
     }
     paste('Available variants: ', paste(sheet.list(), collapse=', '), sep='')
   })
 table = eventReactive(input$submit, {
       table<-data.frame(gi = character(),
                         match = character(),
                         variant = character(),
                         stringsAsFactors = F)
       sheet.list = sheet.list()
       
       for(i in 1:length(sheet.list()))
       {
         sheet = readWorksheet(wb(), sheet.list[i], header=FALSE)
         sheet = cbind(sheet, rep(sheet.list[i],times=nrow(sheet)))
         colnames(sheet) = c('gi', 'match','variant')
         table<-rbind(table,sheet)
       }
       table$variant = as.character(table$variant)
       table = table[!table$gi %in% table[table$variant==as.character(input$zero), 'gi'],]
       variants = strsplit(as.character(input$chosen_variants), split=';')[[1]]
       table = table[table$variant %in% variants,]
       table$gi = trim(table$gi)
       table = table[!duplicated(table[,c('gi', 'variant')]),]
       excluded = strsplit(input$excluded, split=';')[[1]]
       included = strsplit(input$included, split=';')[[1]]
       for (i in excluded){
         table = table[!grepl(i, table$match),]
       }
       for (j in included){
         table = table[grepl(j, table$match),]
       }
       table
 })
 
 dups = reactive({
     if (is.null(table())) {
       return(NULL)
     }
     table = table()
     dups = table
     for (i in unique(dups$gi)) {
       vars = dups[dups$gi==i,'variant']
       dups[dups$gi==i, 'variants'] = paste(vars, collapse ='|')
     }
     dups$variant = NULL
     dups = dups[order(dups$variants),]
     dups = dups[!duplicated(dups),]
     dups = dups[!dups$variants %in% unique(table$variant),]
     dups
   })
 output$dups = renderTable({
     dups()
   })
 
 diagram = reactive({
     venn_data = list()
     if (is.null(table())) {
       return(NULL)
     }
     table = table()
     for (i in unique(table$variant)) {
       venn_data[[i]] = table[table$variant==i,'gi']
     }
  
     diag = venn(venn_data, show.plot=F)
     diag
   })
 
 output$tab_all = renderTable({var
   if (is.null(table())) {
     return(NULL)
   }
   table()
 })
 
 output$venn = renderPlot({
     if (is.null(diagram())) {
       return(NULL)
     }
     plot(diagram(), res=100)
   })
 
   intr = reactive({
       intr = attr(diagram(), 'intersections')
       table = table()
       intersections = data.frame()
       it = 0
       for(i in names(intr)) {
         subset = intr[[i]]
         for (j in subset) {
           it=it+1
           intersections[it, 'gi'] = j
           intersections[it, 'match'] = table[table$gi==j, 'match'][1]
           intersections[it, 'variants'] = i
         }
       }
       intersections[order(intersections$variants),]
   })
   
   
  output$tab_intersections = renderTable({
      intr()
  })
  
  output$download_all = downloadHandler(
    filename = function() {'all.csv'},
    content = function(file){
      write.csv(table(), file)
        }
  )
  output$download_dups = downloadHandler(
    filename = function() {'common.csv'},
    content = function(file){
      write.csv(dups(), file)
    }
  )
  output$download_intr = downloadHandler(
    filename = function() {'intersections.csv'},
    content = function(file){
      write.csv(intr(), file)
    }
  )
}
shinyApp(ui=ui, server=server)