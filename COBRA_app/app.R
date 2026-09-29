library(shiny)
library(bslib)
library(ggplot2)
library(gggenes)
library(stringr)
library(shinyFeedback)
library(data.table)
library(tidyr)
library(dplyr)
library(patchwork)
library(svglite)


options(shiny.maxRequestSize = 5000 * 1024^2)

downloadButton <- function(...) {
    tag <- shiny::downloadButton(...)
    tag$attribs$download <- NULL
    tag
}

ui <- fluidPage(
    shinyFeedback::useShinyFeedback(),
    div(
        style = "margin-top: 15px; margin-bottom: 15px; margin-left: 5px; margin-right: 5px;display: inline-block;",
        tags$img(src = "COBRALogo.png", height = "125px")
    ),
    sidebarLayout(sidebarPanel(
        fileInput(
            "pileupfile",
            "Pileup file with 100 bp windows/binsize",
            accept = NULL,
            buttonLabel = "Browse...",
            placeholder = "No file selected"),
        fileInput(
            "gfffile",
            "Gff file associated with the contigs in the uploaded pileup file",
            accept = NULL,
            buttonLabel = "Browse...",
            placeholder = "No file selected"),
        div(
            helpText(strong(
                "Pileup and gff files can be large and may
                take a minute or two to load in. Please be patient.")
            ), 
            style = "margin-bottom: 20px;" 
        ),
        div(
            helpText(
                "The plot on the bottom is a subset of the larger plot. Annotations 
                are only visible on the contig subset when viewing a <30 kbp contig range."
            ), 
            style = "margin-bottom: 20px;" 
        ),
        selectizeInput(
            "contig",
            "Select contig to display",
            choices = NULL,
            options = list(maxOptions = 5000)
        ),
        numericInput("start", "Start range (bp)", value = 1),
        numericInput("stop", "End range (bp)", value = 25000),
        selectInput(
            "searchCol",
            "Select gff column to display",
            choices = NULL
        ),
        textInput("keywords", label = "Optional: Highlight a gene annotation(s) in the subset plot", value = NULL),
        radioButtons(
            "radio",
            "Select option",
            choices = list("Raw read coverage" = 1, "Log10 read coverage" = 2),
            selected = 1
        ),
        tags$a(href = "https://en.wikipedia.org/wiki/General_feature_format", "What is a gff file?", target = "_blank"),
        tags$br(), 
        tags$a(href = "https://bbmap.org/tools/pileup", "What is a pileup file?", target = "_blank")
    ),mainPanel(
        plotOutput("plots", width = "100%", height = "800px"),
        downloadButton("download", "Download plot.svg")
    ),
    ),
    theme = bslib::bs_theme(preset = "cerulean")
    
)


server <- function(input, output, session) {
    parsed_data <- reactive({
        req(input$pileupfile) 
        pileup <- read.delim(input$pileupfile$datapath, header=FALSE, comment.char="#")
        colnames(pileup) <- c("contigName", "coverage", "position", "tmp")
        contigs <- unique(pileup[,1])
        updateSelectizeInput(session, "contig", choices = contigs, options = list(maxOptions = 5000), server=TRUE)
        return(pileup)
    })
    
    
    parsed_data2 <- reactive({
        req(input$gfffile) 
        id <- showNotification("Processing gff file, please be patient...", duration = NULL, closeButton = FALSE)
        on.exit(removeNotification(id), add = TRUE)
        gff <- read_gff_full(input$gfffile$datapath)
        gff <- as.data.frame(gff)
        return(gff)
    })
    
    observeEvent(parsed_data2(), {
        gff <- parsed_data2()
        columns <- colnames(gff[,c(9:ncol(gff))])
        updateSelectInput(session, "searchCol", choices = columns)
    })
    
    observeEvent(input$contig, {
        updateNumericInput(session, "start", value = 1)
        updateNumericInput(session, "stop", value = 25000)
    })
    
    parsed_data3 <- reactive({
        req(input$contig)
        req(input$searchCol)
        req(parsed_data2())
        gff <- parsed_data2()
        pileup <- parsed_data()
        shinyFeedback::feedbackWarning("gfffile", length(unique(gff[,1]) %in% unique(pileup[,1])==TRUE) != length(unique(gff[,1])),"There are contigs present in the gff file that are not present in the pileup. 
                                           Are you sure the provided gff and pileup are associated with each other?")
        gff <- gff[which(gff$seqid == input$contig), ]
        shinyFeedback::feedbackWarning("gfffile", nrow(gff)==0, "Contig does not exist in gff file.")
        colIdx <- which(colnames(gff) == input$searchCol)
        shinyFeedback::feedbackWarning("stop", input$stop - input$start > 30000, "Choose a range that is <30 kbp to see gene annotations")
        if (input$stop - input$start < 30000) {
            geneAnnotMatches <- gff[which(gff$start %in% seq(input$start, input$stop) &
                                              gff$end %in% seq(input$start, input$stop)),]
            geneAnnotMatches$strand <- ifelse(geneAnnotMatches$strand=="+", TRUE, FALSE)
            if (input$keywords != ""){
                matchIdxs <- str_which(gff[,input$searchCol], regex(paste(input$keywords), ignore_case = TRUE))
                geneAnnotMatchesSubset <- gff[matchIdxs,]
            } else {
                geneAnnotMatchesSubset <- gff[-c(1:nrow(gff)),]
            }
        } else {
            geneAnnotMatches <- gff[-c(1:nrow(gff)),]
            geneAnnotMatchesSubset <- gff[-c(1:nrow(gff)),]
        }
        return(list(geneAnnotMatches, geneAnnotMatchesSubset))
    })
    
    plot2 <- reactive({
        req(parsed_data3())
        my_list <- parsed_data3() 
        geneAnnotMatches <- my_list[[1]]
        geneAnnotMatchesSubset <- my_list[[2]]
        pileup <- parsed_data()
        subset <- pileup[which(pileup[,1]==input$contig),]
        shinyFeedback::feedbackWarning("pileupfile", nrow(subset)==0, "Contig does not exist in pileup file.")
        if(input$radio == 2){
            subset$coverageTypeViral <- abs(log10(subset[, 2]))
            subset[subset == Inf] <- 0
            yaxis <- "Log10 read coverage"
        } else {
            colnames(subset)[2] <- "coverageTypeViral"
            yaxis <- "Read coverage"
        }
        familyplus_pallete <- c( "#00bfff","#e9967a", "#dc143c","#dda0dd","#148890","#ff69b4", "#2f4f4f", "#adff2f","#ff6347","#ff1493","#008080","#562",
                                 "#DCBB26", "#8fbc8f","#2A1C5A","#F9BECA","#591D17","#ffd700","#ffa500","#48d1cc","#9932cc","#824","#6495ed",
                                 "#da70d6", "#008000","#b03060", "#808","#111190","#aa7", "black","#dcdcdc")
        
        geneAnnotMatcheshyp <- geneAnnotMatches[grep("hypothetical", geneAnnotMatches[,input$searchCol], ignore.case = TRUE), ]
        geneAnnotMatchesnohyp <- geneAnnotMatches[grep("hypothetical", geneAnnotMatches[,input$searchCol], ignore.case = TRUE, invert=TRUE), ]
        
        if(nrow(geneAnnotMatcheshyp)==0){
            legendtitle <- input$searchCol
        } else {
            legendtitle <- NULL
        }
        start <- input$start
        stop <- input$stop
        ggplot() +
            geom_area(data = subset, aes(x = position, y = coverageTypeViral), fill = "#009E73") +
            geom_rect(data = geneAnnotMatchesSubset,aes(xmin = start, xmax = end, ymin = 0, ymax = max(subset$coverageTypeViral)),fill = "khaki", alpha = 0.4)+
            geom_gene_arrow(data=geneAnnotMatcheshyp, aes(xmin=start, xmax=end, y=max(subset$coverageTypeViral), color=get(input$searchCol), forward=strand))+
            geom_gene_arrow(data=geneAnnotMatchesnohyp, aes(xmin=start, xmax=end, y=max(subset$coverageTypeViral), fill=get(input$searchCol), forward=strand))+
            labs(
                x = "Basepair position",
                y = yaxis,
                fill="annotation")+
            scale_x_continuous(expand = c(0, 0), limits=c(start,stop)) +
            scale_fill_manual(name=legendtitle, values=familyplus_pallete)+
            scale_color_manual(name=input$searchCol,values="black")+
            guides(
                color = guide_legend(
                    title = input$searchCol, 
                    order = 1), 
                fill = guide_legend(
                    title = legendtitle, 
                    order = 2))+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  text = element_text(size = 15),
                  plot.margin = margin(
                      t = 0,
                      r = 10,
                      b = 0,
                      l = 2),
                  legend.margin = margin(0, 0, 0, 0),
                  legend.spacing.y = unit(0.03, "cm"),
                  legend.title = element_text(size = 12), legend.text= element_text(size = 11))
        
    })
    
    plot1 <- reactive({
        req(parsed_data())
        pileup <- parsed_data()
        if(is.null(input$gfffile)){
            start <- 0
            stop <- 0
        } else {
            start <- input$start
            stop <- input$stop
        }
        subset <- pileup[which(pileup[,1]==input$contig),]
        if(input$radio == 2){
            subset$logcoverage <- abs(log10(subset[, 2]))
            subset[subset == Inf] <- 0
            coverageTypeViral <- subset$logcoverage
            yaxis <- "Log10 read coverage"
        } else {
            colnames(subset)[2] <- "coverageTypeViral"
            yaxis <- "Read coverage"
            
        }
        ggplot() +
            geom_area(data = subset, aes(x = position, y = coverageTypeViral), fill = "#009E73") +
            annotate(
                "rect", 
                xmin = start, xmax = stop, 
                ymin = 0, ymax =max(subset$coverageTypeViral), 
                fill = "lightblue", 
                alpha = 0.2
            )+
            labs(title = input$contig,
                 x = "Basepair position",
                 y = yaxis)+
            scale_x_continuous(expand = c(0, 0)) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  text = element_text(size = 15),
                  plot.margin = margin(
                      t = 0,
                      r = 10,
                      b = 0,
                      l = 2
                  ))
    })
    outputplot <- reactive({
        req(plot1())
        if(!is.null(input$gfffile)){
            plot1 <- plot1()
            plot2 <- plot2()
            return(plot1/plot_spacer()/free(plot2) + plot_layout(heights = c(1, 0.2, 1)))
        } else {
            plot1 <- plot1()
            return(plot1)
        }
    })
    output$plots <- renderPlot({
        outputplot()
    })
    
    output$download <- downloadHandler(
        filename = function() {
            paste0(input$contig, "_COBRAplot.svg")
        },
        content = function(file) {
            ggsave(file, outputplot(), device=svglite::svglite, width=15, height=9, units="in")
        }
    )
}



read_gff_full <- function(file_path) {
    read.delim(
        file_path,
        comment.char = "#",
        header = FALSE,
        sep = "\t",
        col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"),
        stringsAsFactors = FALSE,
        quote = ""
    ) %>%
        mutate(
            start = as.numeric(start),
            end = as.numeric(end),
            score = suppressWarnings(as.numeric(score))
        ) %>%
        separate_longer_delim(attributes, delim = ";") %>%
        filter(attributes != "") %>%
        separate(attributes, into = c("key", "value"), sep = "=", extra = "merge", fill = "right") %>%
        mutate(key = trimws(key), value = trimws(value)) %>%
        pivot_wider(names_from = key, values_from = value, values_fn = ~ paste(.x, collapse = ","))
}


shinyApp(ui = ui, server = server)


