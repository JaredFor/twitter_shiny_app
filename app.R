library(shiny)
library(shinythemes)
library(memoise)
library(wordcloud)
library(tm)

shinyApp(
  ui = tagList(
    shinythemes::themeSelector(),
    navbarPage(
      theme = shinytheme("sandstone"),  # <--- To use a theme, uncomment this sandstone, simplex, slate, spacelab, superhero +, united, yeti
      
      # Title
      "Genetics Company Public Persona",

##--------------------------------------------------------------------------------------------------------------------
      # Tab 4 UI
      
      tabPanel("Navbar 1 - Information/ReadMe",

                 verbatimTextOutput("txtout_verbatim")#,

      ),      
##--------------------------------------------------------------------------------------------------------------------
      # Tab 1 UI
      tabPanel("Navbar 2 - Wordcloud",
               sidebarPanel(
                 radioButtons(inputId = "wordcloud_radio", label = "Company Tweets",
                             choices = list('23andMe'='tw23andMe', 'Ancestry' = 'Ancestry', 'illumina' = 'illumina', 
                                         'Oxford Nanopore' = 'nanopore', 'SOPHiA GENETICS'='SOPHiAGENETICS',
                                         'Veritas Genetics' = 'VeritasGenetics'), selected = 'Ancestry'),
                 numericInput(inputId = "wordcloud_ngram", label = "N-gram Size",
                             value = 1, min = 1, max = 4, step = 1),
                 selectInput(inputId = "wordcloud_weight", label = "N-gram Weight",
                            choices = c('weightTf', 'weightTfIdf', 'weightBin', 'weightSMART'), multiple = FALSE) ,
                 
                 sliderInput(inputId = "wordcloud_slider", label = "Maximum Number of Words:", min = 100, max = 2000,
                             value = 500, step = 50),
                 
                 submitButton("Update")#,

               ),
               
               mainPanel(
                 plotOutput(outputId = "wordcloud_plot")
                 
               )
      ),
      
##--------------------------------------------------------------------------------------------------------------------
      # Tab 2 UI
      tabPanel("Navbar 3 - Commonality and Comparison Wordclouds",
               sidebarPanel(
                 checkboxGroupInput(inputId = "gene_tweets_checkbox", label = "Company Tweets",
                              choices = c('23andMe'='tw23andMe', 'Ancestry' = 'Ancestry', 'illumina', 
                                          'Oxford Nanopore' = 'nanopore', 'SOPHiA GENETICS'='SOPHiAGENETICS',
                                          'Veritas Genetics' = 'VeritasGenetics'), 
                              selected = c('23andMe'='tw23andMe', 'Ancestry' = 'Ancestry')),
                 numericInput(inputId = "commonality_comparison_ngram", label = "N-gram Size",
                              value = 1, min = 1, max = 4, step = 1),
                 selectInput(inputId = "commonality_comparison_weight", label = "N-gram Weight",
                             choices = c('weightTf', 'weightTfIdf', 'weightBin', 'weightSMART'), multiple = FALSE),
                 sliderInput(inputId = "commonality_comparison_slider", label = "Maximum Number of Words:", min = 100, max = 2000,
                             value = 500, step = 50),
                 submitButton("Update")#,
                 
               ),
              
              mainPanel(
                tabsetPanel(
                  tabPanel("Tab 1 - Commonality",
                           plotOutput("Commonality_Wordcloud_plot")),
                  tabPanel("Tab 2 - Comparison",
                           plotOutput("Comparison_Wordcloud_plot")) #,
                )
              )
              ),
      
##--------------------------------------------------------------------------------------------------------------------
      # Tab 3 UI
            
      tabPanel("Navbar 4 - Text Sentiment and Emotion",
               sidebarPanel(
                 radioButtons(inputId = "sentiment_emotion_radio", label = "Company Tweets",
                              choices = c('23andMe'='tw23andMe', 'Ancestry', 'illumina', 
                                          'Oxford Nanopore' = 'nanopore', 'SOPHiA GENETICS'='SOPHiAGENETICS',
                                          'Veritas Genetics' = 'VeritasGenetics')),

                 selectInput(inputId = "sentiment_emotion_lexicon", label = "Lexicon",
                             choices = c("afinn", "bing", "nrc", "loughran"), selected = 'bing', multiple = FALSE),
                 
                 submitButton("Update")#, 

               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Tab 1 - Emotion Plot",
                            plotOutput("sentiment_emotion_radar")#,
                            
                   ),
                   tabPanel("Tab 2 - Sentiment Plot",
                            plotOutput("sentiment_emotion_plot")#,
                 
                   ),
                   tabPanel("Tab 3 - Sentiment Highlighted - HTML",
                            uiOutput("sentiment_emotion_ui"),
                            h5("The HTML will open shortly"))#,
                 
               )
              )
            )
      

    )
  ),
##--------------------------------------------------------------------------------------------------------------------
      

  server = function(input, output, session) {
##--------------------------------------------------------------------------------------------------------------------
    # Tab 1 Server#
    #########
    
    # Tab 1 plot wordcloud
    output$wordcloud_plot <- renderPlot({make_wordcloud(gene_tweets = get(input$wordcloud_radio),
                                                       ngram = (input$wordcloud_ngram),
                                                       weight_type = as.list(input$wordcloud_weight)
   )})

   
##--------------------------------------------------------------------------------------------------------------------
    # Tab 2 Server#
    #########
    
    # Tab 2 plot Commonality wordcloud
    output$Commonality_Wordcloud_plot <- renderPlot({make_commonality_wordcloud(
                    gene_tweets_list = lapply(input$gene_tweets_checkbox, get),
                    list_names = as.vector(input$gene_tweets_checkbox),
                    ngram = (input$commonality_comparison_ngram),
                    weight_type = as.list(input$commonality_comparison_weight)
    )})
    
    # Tab 2 plot comparison wordcloud
    output$Comparison_Wordcloud_plot <- renderPlot({make_comparison_wordcloud(
                    gene_tweets_list = lapply(input$gene_tweets_checkbox, get),
                    list_names = as.vector(input$gene_tweets_checkbox),
                    ngram = (input$commonality_comparison_ngram),
                    weight_type = as.list(input$commonality_comparison_weight)
    )})
 

##--------------------------------------------------------------------------------------------------------------------
    # Tab 3 Server # 
    #########
    output$sentiment_emotion_radar <- renderPlot({make_emo_radar(
                    gene_tweets = get(input$sentiment_emotion_radio)
    )})
    
    # Tab 3 plot Sentiment histogram
    output$sentiment_emotion_plot <- renderPlot({make_sentiment_hist(
                    gene_tweets = get(input$sentiment_emotion_radio),
                    lexicon = (input$sentiment_emotion_lexicon)
    )})
    
    # Tab 3 display html with tweet analysis
    output$sentiment_emotion_ui <- renderUI({make_sentiment(
                    gene_tweets = get(input$sentiment_emotion_radio)
    )})
##--------------------------------------------------------------------------------------------------------------------
    # Tab 4 Server#
    #########

    
    output$txtout_verbatim <- renderPrint({
      cat(paste(" This work pools 6 genetics companies. Five were feature in the article 'The 5 Smartest Companies Analyzing Your DNA' 
        https://www.technologyreview.com/s/608569/the-5-smartest-companies-analyzing-your-dna/,along with Ancestry, as they are
                in a similar market to 23andMe.","
              
              
                Warning!!! Some of these processes take a while.The 'Select Theme' widget was left in so you
                  can see if something is processing or not. 
                (If you change the theme setting and nothing happens, an output is being generated)","
                
                
                Areas for near-term continued improvement:
                -more comments in the code!
                -fix how the emotion radar outputs to RStudio and not the Shiny app
                -give greater control of the hyperparameters within the app, such as the max number of 
                      words and minimum word count needed to be included.
                -the aesthetic values need to either adapt automatically to the n-gram size, etc. (preferred),
                      or need to be controlled from within the app 
                -add a company logo splash page","
                
                Long-term ambition: connect a text input to the app to search the tweet history of any twitter
                      userID and use that within the app.",
                sep = "\n"))
    })
   
  }
)