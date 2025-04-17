#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(readr)
library(ggplot2)
library(plotly)
library(ggbeeswarm)
library(ggdist)
test_small <- read_csv("sampled_results.csv")
test_small <- test_small |> mutate(log10p = -log10(pvalue))
test_small <- test_small |> mutate(FC = 2^abs(log2FoldChange) * ifelse(log2FoldChange > 0, 1, -1))
test_small$key <- test_small$gene
expression_small <- read_csv("small_expression.csv")
candidates <- test_small |>
  filter(padj < 0.05) |>
  filter(log2FoldChange >= 1 | log2FoldChange <=-1)
covariates <- read_csv("covariates.csv")
genes <- test_small$gene

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Explore Expression Data"),
    
    p("Click a point to visualize the normalized data"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(position= "right",
        sidebarPanel(
          
          plotOutput("expression_plot")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("plotly2")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  
  output$plotly2 <- renderPlotly({
    plot_small <- test_small |>
      ggplot(aes(y=log10p, x=log2FoldChange, pval=pvalue, key=key)) +
      geom_point(alpha=0.2) +
      geom_point(aes(y=log10p, x=log2FoldChange, pval=pvalue, key=key), data=candidates, color="blue") +
      geom_hline(yintercept = -log10(0.05)) +
      geom_vline(xintercept = -1) +
      geom_vline(xintercept = 1)
    
    out <- ggplotly(plot_small) 
      event_register(out,'plotly_click')
      out
    }
  )
  
  filtered_data <- reactive({
    row <- event_data("plotly_click")
    row2 <- test_small[test_small$key == row$key,]
    gene <- row2$gene

    out <- expression_small |>
      filter(.feature == gene) |>
      left_join(y=covariates, by=c(sample="geo_accession")) 
    
    out
  })
  
  test_row <- reactive({
    row <- event_data("plotly_click")
    row2 <- test_small[test_small$key == row$key,]

    row2
  })
  
  make_expression_plot <- function(dat, testdat){
    if(is.null(dat)){
      gene <- "Frmd8"
    }
    
    print(testdat)
    
    gene <- dat$.feature[1]

    pval <- testdat$padj
    
    FC <- testdat$log2FoldChange
      
    FC <- 2^abs(FC)
    
    out_plot <-  dat |>
      #filter(.feature == gene) |>
      #filter(time %in% c("Day0", "Day8")) |>
      ggplot() +
      aes(y=expression, x=time) +
      geom_boxplot(aes(fill=time)) +
      geom_beeswarm(cex=5) +
      # stat_halfeye(
      #   aes(fill=time),
      #   # adjust bandwidth
      #   adjust = 0.5,
      #   # move to the right
      #   justification = -0.5,
      #   # remove the slub interval
      #   .width = 0,
      #   point_colour = "black"
      # ) +
      #   facet_wrap(~.feature, scales="free")+
      labs(title=gene, 
           subtitle = paste0( "p=", signif(pval,digits = 2),  
                              ", FC=", signif(FC, digits=2))) + 
      theme_minimal() +
      theme(legend.position = "none") 
    
    out_plot
  }
  
  output$expression_plot <- renderPlot({
    testdat <- test_row()
    make_expression_plot(filtered_data(), testdat)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
