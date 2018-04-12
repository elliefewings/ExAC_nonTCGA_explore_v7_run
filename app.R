#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plyr)
library(dplyr)
library(DT)
library(reshape)
library(ggplot2)

#Load data
load("ExAC_nonTCGA_csqfiltered_shiny_PAV_v5_carrierfreq_20180323.Rdata")



# Define UI for application that draws a histogram
ui <- fluidPage( 
  titlePanel("Explore ExAC non-TCGA release 1.0 data"),
  h4("Protein Affecting Variants: Loss of function, inframe indels, and predicted deleterious and damaging missense variants (by SIFT and PolyPhen respectively)"),
  h4("Loss of function Variants: stop gained, stop lost, start lost, frameshift variant, splice donor and splice acceptor variants"),
  fluidRow(h4("")),
  # Create a new Row in the UI for selectInputs
  fluidRow(column(3, 
                  wellPanel(
                    textInput("gene", "Gene:", ""),
                    actionButton("search", "Search"),
                    checkboxInput("match", "Match gene name exactly", value=FALSE)),
                  wellPanel(
                    selectInput("lof","Variant impact:",c("Protein Affecting", "Loss of Function")),
                    selectInput("af", "Variant Allele Frequency:", c("All", "<0.05", "<0.01")),
                    downloadButton("downloadData", "Download"))),
        column(5, align="center", textOutput("text"), textOutput("text2"), tableOutput("tablesum"),
               tags$style(type="text/css", "#tablesum tr:last-child {font-weight:bold;}")),
        column(4, align="center", plotOutput("summary", width="100%", height="500px"))),
  
  fluidRow(
    DT::dataTableOutput("table"))
  )


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #All data for start
  fullInput <- reactive({
    data <- vv
    if (input$lof != "Protein Affecting") {
      data <- data %>% filter(IMPACT == "HIGH")
    }
    if (input$af == "<0.05") {
      data <- data %>% filter(AF < 0.05)
    }
    if (input$af == "<0.01") {
      data <- data %>% filter(AF < 0.01)
    }
    data
  })
  
  #Filter on gene
  geneInput <- eventReactive(input$search, {
    data <- vv
    if (input$gene != "" & input$match == FALSE) {
      data <- data %>% filter(grepl(toupper(input$gene), SYMBOL))
    }
    if (input$gene != "" & input$match == TRUE) {
      data <- data %>% filter(SYMBOL == toupper(input$gene))
    }
  data
  })
  
  #Set title
  geneout <- eventReactive(input$search, {
    if (input$gene != "") {
      genen <- toupper(input$gene)
      genen}
  })
  

  # Filter data based on selections
  datasetInput <- reactive({
    data <- geneInput()
    if (input$lof != "Protein Affecting") {
      data <- data %>% filter(IMPACT == "HIGH")
    }
    if (input$af == "<0.05") {
      data <- data %>% filter(AF < 0.05)
    }
    if (input$af == "<0.01") {
      data <- data %>% filter(AF < 0.01)
    }
    data})
  
  output$summary <- renderPlot({
    data <- geneInput()
    out <- data.frame(impact=c("LoF", "Missense"), 
                        common=c(sum(data[data$IMPACT == "HIGH" & data$AF >= 0.01,]$AC, na.rm=TRUE), sum(data[grepl("missense", data$Consequence)& data$AF >= 0.01,]$AC, na.rm=TRUE)),
                        rare5=c(sum(data[data$IMPACT == "HIGH" & data$AF < 0.01,]$AC, na.rm=TRUE), sum(data[grepl("missense", data$Consequence)& data$AF < 0.01,]$AC, na.rm=TRUE)))
      out.melt <- melt(out)
      title <- geneout()
      ggplot(out.melt, aes(impact, value, fill=variable, label=value)) + geom_bar(stat="identity") + 
        labs(title=title, x="Variant consequence", y="Allele Count") +
        scale_fill_manual(name="Variant Allele Frequency", labels=c("AF >= 0.01", "AF < 0.01"), values = c("orange2", "deepskyblue3")) +
        geom_text(data=subset(out.melt, value != 0), size = 5, position = position_stack(vjust = 0.5)) +
        theme(axis.text = element_text(size=12), axis.title = element_text(size=14), legend.text = element_text(size = 14), legend.title = element_text(size=14), 
              title=element_text(size=16, face="bold"))
  })
  
  output$tablesum <- renderTable({
    data <- geneInput()
      if (input$af == "<0.05") {
        data <- data %>% filter(AF < 0.05)
      }
      if (input$af == "<0.01") {
        data <- data %>% filter(AF < 0.01)
      }
      out <- data.frame(LoF=c(sum(data[data$IMPACT == "HIGH",]$AC_AFR, na.rm=TRUE),
                              sum(data[data$IMPACT == "HIGH",]$AC_EAS, na.rm=TRUE),
                              sum(data[data$IMPACT == "HIGH",]$AC_FIN, na.rm=TRUE),
                              sum(data[data$IMPACT == "HIGH",]$AC_NFE, na.rm=TRUE),
                              sum(data[data$IMPACT == "HIGH",]$AC_AMR, na.rm=TRUE),
                              sum(data[data$IMPACT == "HIGH",]$AC_SAS, na.rm=TRUE),
                              sum(data[data$IMPACT == "HIGH",]$AC_OTH, na.rm=TRUE),
                              sum(data[data$IMPACT == "HIGH",]$AC, na.rm=TRUE)),
                        Miss=c(sum(data[grepl("missense", data$Consequence),]$AC_AFR, na.rm=TRUE),
                               sum(data[grepl("missense", data$Consequence),]$AC_EAS, na.rm=TRUE),
                               sum(data[grepl("missense", data$Consequence),]$AC_FIN, na.rm=TRUE),
                               sum(data[grepl("missense", data$Consequence),]$AC_NFE, na.rm=TRUE),
                               sum(data[grepl("missense", data$Consequence),]$AC_AMR, na.rm=TRUE),
                               sum(data[grepl("missense", data$Consequence),]$AC_SAS, na.rm=TRUE),
                               sum(data[grepl("missense", data$Consequence),]$AC_OTH, na.rm=TRUE),
                               sum(data[grepl("missense", data$Consequence),]$AC, na.rm=TRUE)),
                        an=c(mean(data$AN_AFR, na.rm=TRUE),
                             mean(data$AN_EAS, na.rm=TRUE),
                             mean(data$AN_FIN, na.rm=TRUE),
                             mean(data$AN_NFE, na.rm=TRUE),
                             mean(data$AN_AMR, na.rm=TRUE),
                             mean(data$AN_SAS, na.rm=TRUE),
                             mean(data$AN_OTH, na.rm=TRUE),
                             mean(data$AN, na.rm=TRUE)))      
      colnames(out) <- c("LoF AC", "Missense AC", "Average AN")
      row.names(out) <- c("African", "East Asian", "European (Finnish)", "European (Non-Finnish)", "Latino", "South Asian", "Other", "Total")
      out
  }, spacing = "l", rownames=TRUE, digits=0, hover=TRUE)
  
  #Set datatable column name explanations 
  sketch = htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th('', title = 'Row Names'),
        th('CHROM', title = 'Chromosome'),
        th('POS', title = 'Genomic Position'),
        th('REF', title = 'Reference allele'),
        th('ALT', title = 'Alternative allele'),
        th('Carrier_count', title = 'Count of individuals carrying alternative allele'),
        th('Carrier_frequency', title = 'Frequency of individuals carrying alternative allele'),
        th('AC', title = 'Allele Count (Het=1, Hom=2)'),
        th('AF', title = 'Allele Frequency'),
        th('AN', title = 'Allele Number'),
        th('Consequence', title = 'Variant consequence as predicted by VEP'),
        th('IMPACT', title = 'Variant impact as predicted by VEP'),
        th('SYMBOL', title = 'Gene symbol'),
        th('cDNA_position', title = 'Variant position in cDNA'),
        th('CDS_position', title = 'Variant position in coding sequence'),
        th('Protein_position', title = 'Position of affected amino acid'),
        th('Amino_acids', title = 'Amino acid produced by reference and alternative'),
        th('Existing_variation', title = 'Variant rsID'),
        th('HGVSc', title = 'Variant nomenclature according to HGVS'),
        th('SIFT', title = 'Prediction of variant effect by SIFT'),
        th('PolyPhen', title = 'Prediction of variant effect by PolyPhen'),
        th('CLIN_SIG', title = 'Prediction of variant clinical significance'),
        th('AC_FEMALE', title = 'Allele count in females'),
        th('AC_MALE', title = 'Allele count in males'),
        th('AC_AFR', title = 'Allele count in African/African American population'),
        th('AC_AMR', title = 'Allele count in Latino population'),
        th('AC_EAS', title = 'Allele count in East Asian population'),
        th('AC_FIN', title = 'Allele count in Finnish population'),
        th('AC_NFE', title = 'Allele count in Non-Finnish European population'),
        th('AC_SAS', title = 'Allele count in South Asian population'),
        th('AC_OTH', title = 'Allele count in other populations'),
        th('AN_FEMALE', title = 'Allele number in females'),
        th('AN_MALE', title = 'Allele number in males'),
        th('AN_AFR', title = 'Allele number in African/African American population'),
        th('AN_AMR', title = 'Allele number in Latino population'),
        th('AN_EAS', title = 'Allele number in East Asian population'),
        th('AN_FIN', title = 'Allele number in Finnish population'),
        th('AN_NFE', title = 'Allele number in Non-Finnish European population'),
        th('AN_SAS', title = 'Allele number in South Asian population'),
        th('AN_OTH', title = 'Allele number in other populations'),
        th('AF_FEMALE', title = 'Allele frequency in females'),
        th('AF_MALE', title = 'Allele frequency in males'),
        th('AF_AFR', title = 'Allele frequency in African/African American population'),
        th('AF_AMR', title = 'Allele frequency in Latino population'),
        th('AF_EAS', title = 'Allele frequency in East Asian population'),
        th('AF_FIN', title = 'Allele frequency in Finnish population'),
        th('AF_NFE', title = 'Allele frequency in Non-Finnish European population'),
        th('AF_SAS', title = 'Allele frequency in South Asian population'),
        th('AF_OTH', title = 'Allele frequency in other populations')
      )
    )
  ))
  ###Set main table based on whether genes are searched
  
  output$table <- DT::renderDataTable({
    DT::datatable(fullInput(), container=sketch, options = list(lengthMenu = c(10, 100, 200), pageLength = 10, columnDefs = list(list(
      targets=c(3, 4),
      render = JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data.length > 5 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, 5) + '...</span>' : data;",
        "}")
    ),
    list(
      targets=10,
      render = JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data.indexOf('&') > 0 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, data.indexOf('&')) + '...</span>' : data;",
        "}")))))
  })
  
  observeEvent(input$search, {
    output$table <- DT::renderDataTable({
      DT::datatable(datasetInput(), container=sketch, options = list(lengthMenu = c(10, 100, 200), pageLength = 10, columnDefs = list(list(
        targets=c(3, 4),
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 5 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 5) + '...</span>' : data;",
          "}")
      ),
      list(
        targets=10,
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.indexOf('&') > 0 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, data.indexOf('&')) + '...</span>' : data;",
          "}")))))
    })
  })
  
  
  ###Set text for headers
  
  observeEvent(input$search, output$text <-renderText({
    title <- paste("Occurences of variants in gene: ", geneout(), sep="")
    title
  }))
  
  observeEvent(input$search, {
    output$text2 <- renderText({
      taf <- "All variants selected"
      if (input$af != "All") {
        taf <- paste("Variants with AF ", input$af, " selected", sep="")}
      else {
        taf <- "All variants selected" 
      }
      taf
    })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ExAC_nonTCGA_", input$gene, "_", input$lof, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

