setwd("C:\\Users\\Dolomite\\Desktop\\Analytics\\902\\project")

library(readtext)
library(qdap)
library(readr)
library(tm)
library(wordcloud)
library(wordcloud2)
library(RWeka)
library(stringr)
library(dplyr)
library(shiny)
library(stringr)
library(sentimentr)
library(lexicon)
library(tidyr)
library(tidytext)
library(radarchart)
library(memoise)
library(ggplot2)

### Read the tweets csv. To make a word cloud, put this in the make_wordcloud function for 'gene_tweets'
gene_tweets <- read_csv("genetics_companies_2000each.csv", col_names = TRUE, trim_ws = TRUE)
tw23andMe <- read_csv("genetics_companies_2000tw23andMe.csv", col_names = TRUE, trim_ws = TRUE)
illumina <- read_csv("genetics_companies_2000illumina.csv", col_names = TRUE, trim_ws = TRUE)
nanopore <- read_csv("genetics_companies_2000nanopore.csv", col_names = TRUE, trim_ws = TRUE)
SOPHiAGENETICS <- read_csv("genetics_companies_2000SOPHiAGENETICS.csv", col_names = TRUE, trim_ws = TRUE)
VeritasGenetics <- read_csv("genetics_companies_2000VeritasGenetics.csv", col_names = TRUE, trim_ws = TRUE)
Ancestry <- read_csv("genetics_companies_2000Ancestry.csv", col_names = TRUE, trim_ws = TRUE)


# Convert tweet text format
# gene_tweets$text <- iconv(gene_tweets$text, from = "UTF-8", to = "ASCII", sub = "")
# tw23andMe$text <- iconv(tw23andMe$text, from = "UTF-8", to = "ASCII", sub = "")
# illumina$text <- iconv(illumina$text, from = "UTF-8", to = "ASCII", sub = "")
# nanopore$text <- iconv(nanopore$text, from = "UTF-8", to = "ASCII", sub = "")
# SOPHiAGENETICS$text <- iconv(SOPHiAGENETICS$text, from = "UTF-8", to = "ASCII", sub = "")
# VeritasGenetics$text <- iconv(VeritasGenetics$text, from = "UTF-8", to = "ASCII", sub = "")
# Ancestry$text <- iconv(Ancestry$text, from = "UTF-8", to = "ASCII", sub = "")

##--------------------------------------------------------------------------------------------------------------------

# Cleaning corpus / pre_processing for the text-- The order of these processes matters, as results will change
# (stemming is disabled because of the poor aesthetics in the word cloud of stemmed words
# for better information, stemming should probably be in place, followed by good stem completion)
# NOTE: the first tm_mp will call 'corpus', all subsequent calls should be to 'cleaned_corpus'
#
# The phrase: "Thx for tweeting #MyDNAStory! Use this link for a special promotion to unbox your DNA Story with your family!"
# was repeated many times by the 23andMe twitter, it is effectively removed in the 2nd custom stopwords (placing it in the 
# first one did not work)
clean_corpus <- function(corpus){
  custom_stop_words <- c("amp", "can", "may", "via", "one", "httpstcodggafu", "youre", "nanopore", "nanoporeconf",
                         "illumina", "andme", "andmes", "sophia", "veritas", "veritasgenetics")
  cleaned_corpus <- tm_map(corpus, removeWords, custom_stop_words)
  cleaned_corpus <- tm_map(cleaned_corpus, removeWords, stopwords(kind="en"))
  cleaned_corpus <- tm_map(cleaned_corpus, removeWords, stopwords(kind="fr"))
  cleaned_corpus <- tm_map(cleaned_corpus, content_transformer(replace_abbreviation))
  cleaned_corpus <- tm_map(cleaned_corpus, content_transformer(replace_ordinal))
  cleaned_corpus <- tm_map(cleaned_corpus, content_transformer(tolower))
  cleaned_corpus <- tm_map(cleaned_corpus, removePunctuation)
  cleaned_corpus <- tm_map(cleaned_corpus, removeNumbers)
  custom_stop_words2 <- c("nanopore", "illumina", "andme", "sophia", "sophiagenetics", "veritas",
                          "ufc thx tweeting mydnastory use link special promotion unbox dna story family",
                          "ufc", "uff", "uufef", "ufd", "ufe", "ddm", "de", "la")
  cleaned_corpus <- tm_map(cleaned_corpus, removeWords, custom_stop_words2)
  # cleaned_corpus <- tm_map(cleaned_corpus, stemDocument)  # stemming looks bad in the wordcloud without stemCompletion
  # cleaned_corpus <- tm_map(cleaned_corpus, stemCompletion(dictionary = cleaned_corpus, type='prevalent'))
  cleaned_corpus <- tm_map(cleaned_corpus, removeWords, stopwords(kind="en"))  #repeated on purpose
  cleaned_corpus <- tm_map(cleaned_corpus, stripWhitespace)
  return(cleaned_corpus)
}


##--------------------------------------------------------------------------------------------------------------------


# this function processes and cleans the original csv of the tweets
clean_tweet_corpus <- function(gene_tweets){
  
  # Convert tweet text
  gene_tweets$text <- iconv(gene_tweets$text, from = "UTF-8", to = "ASCII", sub = "")
  
  # remove https addresses
  gene_tweets$text <- str_replace_all(string = gene_tweets$text, pattern =  "(https://(t.co/)[:alnum:]*)", replacement = "")
  
  # 
  # # Format the tweet time and remove web links and retain only alphanumeric characters
  # gene_tweets_clean <- gene_tweets %>%
  #   mutate(gene_tweets$created_at <- as.POSIXct(gene_tweets$created_at, format = "%a %b %d %H:%M:%S +0000 %Y"),
  #          gene_tweets$text <- gsub("\\s?(f|ht)(tp)(s?)(://)([^\\.]*)[\\.|/](\\S*)", "", gene_tweets$text),
  #          gene_tweets$text <- str_replace_all(gene_tweets$text, "[^A-Za-z0-9 _.,!"'/$@#:;%^&*_+-=\(\)]", " "))
  # 
  # # Create a corpus
  # tweet_corpus <- VCorpus(VectorSource(gene_tweets_clean$text))
  
  # Create a corpus
  tweet_corpus <- VCorpus(VectorSource(gene_tweets$text))
  
  cleaned_tweet_corpus <- clean_corpus(tweet_corpus)
  return(cleaned_tweet_corpus)
}

##--------------------------------------------------------------------------------------------------------------------


# This function uses cleans a text with the clean_tweet_corpus function, tokenizes according to the n-gram size
# creates a TDM and word frequency data.frame
# for best results, ngram should be 1, 2 or 3
# Available 'weight_type' weighting functions shipped with the tm package are weightTf, weightTfIdf, weightBin, and weightSMART.

make_wordcloud <- function(gene_tweets, ngram, weight_type='weightTf'){

  cleaned_tweet_corpus <- clean_tweet_corpus(gene_tweets)

  # function to tokenize the cleaned corpus according to the desired n-gram
  tokenizer <- function(x)
    NGramTokenizer(x,Weka_control(min=ngram, max=ngram))
  
  # create TDM from cleaned corpus of desired n-grams
  TDM_gene <- TermDocumentMatrix(cleaned_tweet_corpus,control = list(tokenize=tokenizer, weighting=weight_type))
  TDM_gene_m <- as.matrix(TDM_gene)
  
  # calculate term frequency and sort
  term_frequency <- rowSums(TDM_gene_m)
  term_frequency <- sort(term_frequency,dec=TRUE)
  
  # cleaned_tweet_corpus <- clean_corpus(gene_tweets)
  
  # term = names(term_frequency)
  # num = term_frequency
  # wordcloud(term, num, random.order = FALSE, vfont=c("serif","bold"), min.freq=15,
  #           max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)

  word_freqs <- data.frame(term = names(term_frequency), num = term_frequency)
  wordcloud(word_freqs$term, word_freqs$num, random.order = FALSE, vfont=c("serif","bold"), min.freq=15,
            max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
  
  # return(word_freqs)
}


##--------------------------------------------------------------------------------------------------------------------

# Similar to make_wordcloud, but with different pre-processing due to the multiple text sets
# This function uses cleans a text with the clean_tweet_corpus function, tokenizes according to the n-gram size
# creates a TDM and word frequency data.frame
# For best results, ngram should be 1, 2 or 3
# Available 'weight_type' weighting functions shipped with the tm package are weightTf, weightTfIdf, weightBin, and weightSMART.
# gene_tweets_list must be in the  form  list(company1, company2, company3)            
# list_names must be in the form  c('company1', 'company2', 'company3')

make_commonality_wordcloud <- function(gene_tweets_list, list_names, ngram, weight_type='weightTf'){

  # flatten the tweets first
  tweets <- c()

  for (i in c(1:length(gene_tweets_list))) {
    flat_tweets_list <- c(str_flatten(gene_tweets_list[[i]]$text, collapse = " "))
    tweets <-append(tweets, flat_tweets_list)}

  # remove https addresses
  tweets <- str_replace_all(string = tweets, pattern =  "(https://(t.co/)[:alnum:]*)", replacement = "")
  
  tweet_corpus <- VCorpus(VectorSource(tweets))
  
 # create a corpus and preprocess it
  cleaned_tweet_corpus <- clean_corpus(tweet_corpus)
  
  # function to tokenize the cleaned corpus according to the desired n-gram
  tokenizer <- function(x)
    NGramTokenizer(x,Weka_control(min=ngram, max=ngram))
  
  # create TDM from cleaned corpus of desired n-grams
  TDM_gene <- TermDocumentMatrix(cleaned_tweet_corpus,control = list(tokenize=tokenizer, weighting=weight_type))
  TDM_gene_m <- as.matrix(TDM_gene)
  
  # Create a commonality cloud
  commonality.cloud(TDM_gene_m, colors=brewer.pal(8, "Dark2"),max.words = 200, title.size = 2, scale = c(3, .3),
                    match.colors = TRUE,random.order=FALSE,rot.per=0, fixed.asp = FALSE)
  
  }

# gene_tweets_list <- list(tw23andMe, nanopore, VeritasGenetics)
# list_names <- c('tw23andMe', 'nanopore', 'VeritasGenetics')
# make_commonality_wordcloud(gene_tweets_list, list_names, 1)
##--------------------------------------------------------------------------------------------------------------------


# Similar to make_wordcloud, but with different pre-processing due to the multiple text sets
# This function uses cleans a text with the clean_tweet_corpus function, tokenizes according to the n-gram size
# creates a TDM and word frequency data.frame
# For best results, ngram should be 1, 2 or 3
# Available 'weight_type' weighting functions shipped with the tm package are weightTf, weightTfIdf, weightBin, and weightSMART.
# gene_tweets_list must be in the vector form  list(company1, company2, company3)            
# list_names must be in the vector form  c('company1', 'company2', 'company3')

make_comparison_wordcloud <- function(gene_tweets_list, list_names, ngram, weight_type='weightTf'){
  
  # flatten the tweets first
  tweets <- c()
  
  for (i in c(1:length(gene_tweets_list))) {
    flat_tweets_list <- c(str_flatten(gene_tweets_list[[i]]$text, collapse = " "))
    tweets <-append(tweets, flat_tweets_list)}
  
  # remove https addresses
  tweets <- str_replace_all(string = tweets, pattern =  "(https://(t.co/)[:alnum:]*)", replacement = "")

  tweet_corpus <- VCorpus(VectorSource(tweets))
  
  # create a corpus and preprocess it
  cleaned_tweet_corpus <- clean_corpus(tweet_corpus)
  
  # function to tokenize the cleaned corpus according to the desired n-gram
  tokenizer <- function(x)
    NGramTokenizer(x,Weka_control(min=ngram, max=ngram))
  
  # create TDM from cleaned corpus of desired n-grams
  TDM_gene <- TermDocumentMatrix(cleaned_tweet_corpus,control = list(tokenize=tokenizer, weighting=weight_type))
  colnames(TDM_gene) <- list_names
  TDM_gene_m <- as.matrix(TDM_gene)

  # create the comparison cloud
  comparison.cloud(TDM_gene_m, colors=brewer.pal(8, "Dark2"),max.words = 200, title.size = 2, scale = c(3, .3),
                   match.colors = TRUE,random.order=FALSE,rot.per=0, fixed.asp = FALSE)

}

##--------------------------------------------------------------------------------------------------------------------

### SENTIMENT ANALYSIS

# only the bing or loughran lexicon can be used with this function for meaningful results
make_sentiment_hist <- function(gene_tweets, lexicon="bing"){
  
  # # Format the tweet time and remove web links and retain only alphanumeric characters
  # gene_tweets_clean <- gene_tweets %>%
  #   mutate(gene_tweets$created_at <- as.POSIXct(gene_tweets$created_at, format = "%a %b %d %H:%M:%S +0000 %Y"),
  #          gene_tweets$text <- gsub("\\s?(f|ht)(tp)(s?)(://)([^\\.]*)[\\.|/](\\S*)", "", gene_tweets$text),
  #          gene_tweets$text <- str_replace_all(gene_tweets$text, "[^[:alnum:]]", " "))
  
  cleaned_tweet_corpus <- clean_tweet_corpus(gene_tweets)
  
  tidy_mytext <- tidy(TermDocumentMatrix(cleaned_tweet_corpus))
  sentiment_lex <- get_sentiments(lexicon)
  mytext_sentiment_lex <- inner_join(tidy_mytext, sentiment_lex, by = c("term" = "word"))
  mytext_sentiment_lex$sentiment_n <- ifelse(mytext_sentiment_lex$sentiment=="negative", -1, 1)
  mytext_sentiment_lex$sentiment_score <- mytext_sentiment_lex$count*mytext_sentiment_lex$sentiment_n
  
  sentiment_summary <- mytext_sentiment_lex %>%
    group_by(document) %>%
    summarize(tweet_sentiment = sum(sentiment_score)) %>%
    arrange(desc(tweet_sentiment))
  
  par(mar=c(5,4,4,2))
  hist(x = sentiment_summary$tweet_sentiment, breaks = 11, xlim = c(-5,5),
       main = c("Tweet Sentiment"), xlab = c("Sentiment Intensity"))

}

##--------------------------------------------------------------------------------------------------------------------
# Sentiment using sentimentr
# This produces an html page with the text highlighted as good, bad or neutral and gives a cumulative score for the document.

make_sentiment <- function(gene_tweets){
  
  # if umlauts and aceents are a problem then maybe stringi::stri_trans_general(fruits, 'latin-ascii')
  # SO FAR THIS ISN"T NEEDEDgene_tweets$created_at <- as.POSIXct(x= gene_tweets$created_at, format = "%a %b %d %H:%M:%S +0000 %Y")
  
  # remove https addresses
  gene_tweets$text <- str_replace_all(string = gene_tweets$text, pattern =  "(https://(t.co/)[:alnum:]*)", replacement = "")
  
  sentiment(gene_tweets$text, polarity_dt = lexicon::hash_sentiment_jockers_rinker,
            valence_shifters_dt = lexicon::hash_valence_shifters,
            amplifier.weight = 0.8, n.before = 5, n.after = 2,
            question.weight = 1, adversative.weight = 0.85)
  text_sentiment <- extract_sentiment_terms(gene_tweets$text)

  out <- sentiment_by(gene_tweets$text)
  highlight(out)

  return(text_sentiment)
}

##--------------------------------------------------------------------------------------------------------------------
# Create an emotional radar

make_emo_radar <- function(gene_tweets){
  
  cleaned_tweet_corpus_emo <-clean_tweet_corpus(gene_tweets)

  tidy_tweet_emo <- tidy(TermDocumentMatrix(cleaned_tweet_corpus_emo))
 
  nrc_lex <- get_sentiments("nrc")
  mytext_nrc <- inner_join(tidy_tweet_emo, nrc_lex, by = c("term" = "word"))
  # remove positive and negative first
  mytext_nrc_noposneg <- mytext_nrc[!(mytext_nrc$sentiment %in% c("positive","negative")),]
  emotion_summary <- mytext_nrc_noposneg %>%
    group_by(document,sentiment) %>%
    summarize(review_sentiment = sum(count)) %>%
    arrange(desc(review_sentiment))
  
  emotion_overall_summary <- mytext_nrc_noposneg %>%
    group_by(sentiment) %>%
    summarize(review_sentiment = sum(count)) %>%
    arrange(desc(review_sentiment))
  
  chartJSRadar(emotion_overall_summary)
}

##--------------------------------------------------------------------------------------------------------------------
# # Testing area
# gene_tweets_list <- list(tw23andMe, nanopore, VeritasGenetics)
# list_names <- c('tw23andMe', 'nanopore', 'VeritasGenetics')
# 
# make_emo_radar(tw23andMe)
# make_emo_radar(illumina)
# make_emo_radar(nanopore)
# make_emo_radar(SOPHiAGENETICS)
# make_emo_radar(VeritasGenetics)
# make_emo_radar(Ancestry)
# 
# make_wordcloud(tw23andMe, 1, weightTf)
# make_wordcloud(SOPHiAGENETICS, 2)
# make_wordcloud(make_wordcloud(tw23andMe, 1))
# make_wordcloud(Ancestry, 1, "weightTfIdf")
# make_wordcloud(Ancestry, 2, "weightTfIdf")
# 
# make_commonality_wordcloud(gene_tweets_list, list_names, 1)
# 
# make_comparison_wordcloud(gene_tweets_list, list_names, 1)
# 
# make_sentiment_hist(illumina)
# make_sentiment_hist(illumina, "loughran")
# 
# make_sentiment(illumina)
# make_sentiment(tw23andMe)
# 
# 
# # This might be need to create the right graphical parameters for the word cloud, default par(mar=c(5,4,4,2))
# # there doesn't seem to be much difference in the output plot
# 
# par(mar=c(1,1,1,1))
# 
# # An example word cloud. Default weight is term frequency, weightTf
# make_wordcloud(illumina,2)
# make_wordcloud(nanopore,1)
# make_wordcloud(nanopore,1,weightTfIdf)
# 
# 
# 
# ## This returns a dropdown list of all top terms with quoted area greyed out
# # global <- company <<- list("23andMe" = "word_freqs_tw23andMe",
# #                            "illumina" = cleaned_tweet_corpus_illumina,
# #                            "Oxford Nanopore" = cleaned_tweet_corpus_nanopore,
# #                            "SOPHiAGENETICS" = cleaned_tweet_corpus_SOPHiAGENETICS,
# #                            "Veritas Genetics" = cleaned_tweet_corpus_VeritasGenetics,
# #                            "Ancestry" = cleaned_tweet_corpus_Ancestry
# # )
# 
# company <<- list("23andMe" = tw23andMe,
#                            "illumina" = illumina,
#                            "Oxford Nanopore" = nanopore,
#                            "SOPHiAGENETICS" = SOPHiAGENETICS,
#                            "Veritas Genetics" = VeritasGenetics, 
#                            "Ancestry" = Ancestry
# )
# ############Word Cloud
# 
# # # Create word_freqs
# # word_freqs_tw23andMe <- data.frame(term = names(cleaned_tweet_corpus_tw23andMe), num = cleaned_tweet_corpus_tw23andMe)
# # word_freqs_illumina <- data.frame(term = names(cleaned_tweet_corpus_illumina), num = cleaned_tweet_corpus_illumina)
# # word_freqs_nanopore <- data.frame(term = names(cleaned_tweet_corpus_nanopore), num = cleaned_tweet_corpus_nanopore)
# # word_freqs_SOPHiAGENETICS <- data.frame(term = names(cleaned_tweet_corpus_SOPHiAGENETICS), num = cleaned_tweet_corpus_SOPHiAGENETICS)
# # word_freqs_VeritasGenetics <- data.frame(term = names(cleaned_tweet_corpus_VeritasGenetics), num = cleaned_tweet_corpus_VeritasGenetics)
# # word_freqs_Ancestry <- data.frame(term = names(cleaned_tweet_corpus_Ancestry), num = cleaned_tweet_corpus_Ancestry)
# # 
# # word_freqs_gene_tweets <- data.frame(term = names(cleaned_tweet_corpus_gene_tweets), num = cleaned_tweet_corpus_gene_tweets)
# # 
# # # Create a wordcloud for the values in word_freqs
# # par(mar=rep(0, 4))
# # # wordcloud(word_freqs$term, word_freqs$num, random.order = FALSE, vfont=c("serif","bold"), scale = c(7,1), min.freq=25, max.words=300,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# # wordcloud(word_freqs_gene_tweets$term, word_freqs_gene_tweets$num, random.order = FALSE, vfont=c("serif","bold"), scale = c(7,1), min.freq=25, max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# # 
# # wordcloud(word_freqs_tw23andMe$term, word_freqs_tw23andMe$num, random.order = FALSE, vfont=c("serif","bold"), min.freq=15, max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# # wordcloud(word_freqs_illumina$term, word_freqs_illumina$num, random.order = FALSE, vfont=c("serif","bold"), scale = c(1,1), min.freq=5, max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# # wordcloud(word_freqs_nanopore$term, word_freqs_nanopore$num, random.order = FALSE, vfont=c("serif","bold"), min.freq=25, max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# # wordcloud(word_freqs_SOPHiAGENETICS$term, word_freqs_SOPHiAGENETICS$num, random.order = FALSE, vfont=c("serif","bold"), scale = c(1,1), min.freq=25, max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# # wordcloud(word_freqs_VeritasGenetics$term, word_freqs_VeritasGenetics$num, random.order = FALSE, vfont=c("serif","bold"), scale = c(1,1), min.freq=25, max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# # wordcloud(word_freqs_Ancestry$term, word_freqs_Ancestry$num, random.order = FALSE, vfont=c("serif","bold"), scale = c(1,1), min.freq=25, max.words=100,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
# 
# 
# 
# illumina_sentiment <- make_sentiment(illumina)
# 
# # make_sentiment <- function(tweet_text){
# #   sentiment(tweet_text, polarity_dt = lexicon::hash_sentiment_jockers_rinker,
# #             valence_shifters_dt = lexicon::hash_valence_shifters,
# #             amplifier.weight = 0.8, n.before = 5, n.after = 2,
# #             question.weight = 1, adversative.weight = 0.85)
# #   text_sentiment <- extract_sentiment_terms(tweet_text)
# #   
# #   out <- sentiment_by(tweet_text)
# #   highlight(out)
# #   
# #   return(text_sentiment)
# # }
# # this text should be preprocessed first to improve the outcome ***************************************
# tw23andMe_sentiment <- make_sentiment(tw23andMe$text)
# illumina_sentiment <- make_sentiment(illumina$text)
# nanopore_sentiment <- make_sentiment(nanopore$text)
# SOPHiAGENETICS_sentiment <- make_sentiment(SOPHiAGENETICS$text)
# VeritasGenetics_sentiment <- make_sentiment(VeritasGenetics$text)
# Ancestry_sentiment <- make_sentiment(Ancestry$text)
# 
# check <- str_flatten(tw23andMe$text, collapse = " ")
# # check2 <- str_flatten(cleaned_tweet_corpus, collapse = " ")
# 
# make_sentiment(check)
# ###########################################################################################################
# ###########################################################################################################
# ###########################################################################################################
# ########## Emotions #############
# # afinn_lex <- get_sentiments("afinn")
# # table(afinn_lex$score)
# # 
# # nrc_lex <- get_sentiments("nrc")
# # table(nrc_lex$sentiment)
# # 
# # loughran_lex <- get_sentiments("loughran")
# # table(loughran_lex$sentiment)
# 
# # EMOTIONS
# # mytext <- read.csv("app_reviews.csv")
# # mytext$Review <- iconv(mytext$Review, from = "UTF-8", to = "ASCII", sub = "")
# 
# 
# # clean_corpus <- function(corpus){
# #   cleaned_corpus <- tm_map(corpus, content_transformer(replace_abbreviation))
# #   cleaned_corpus <- tm_map(corpus, content_transformer(replace_ordinal))
# #   cleaned_corpus <- tm_map(cleaned_corpus, content_transformer(tolower))
# #   cleaned_corpus <- tm_map(cleaned_corpus, removePunctuation)
# #   cleaned_corpus <- tm_map(cleaned_corpus, removeNumbers)
# #   # cleaned_corpus <- tm_map(cleaned_corpus, stemDocument)
# #   # cleaned_corpus <- tm_map(cleaned_corpus, stemCompletion(dictionary = cleaned_corpus, type='prevalent'))
# #   custom_stop_words <- c("amp", "can", "may", "via", "one", "httpstcodggafu", "youre", "nanopore", "nanoporeconf","illumina", "andme", "andmes", "sophia", "veritas", "veritasgenetics")
# #   cleaned_corpus <- tm_map(cleaned_corpus, removeWords, stopwords(kind="en"))
# #   cleaned_corpus <- tm_map(cleaned_corpus, removeWords, custom_stop_words)
# #   custom_stop_words2 <- c("nanopore", "illumina", "andme", "sophia", "sophiagenetics", "veritas")
# #   cleaned_corpus <- tm_map(cleaned_corpus, removeWords, custom_stop_words2)
# #   cleaned_corpus <- tm_map(cleaned_corpus, stripWhitespace)
# #   return(cleaned_corpus)
# # }
# ##############################################3
# 
# # tweet_corpus <- VCorpus(VectorSource(tw23andMe$text))
# # 
# # cleaned_tweet_corpus <- clean_corpus(tweet_corpus)
# ######
# # output of the function should be the word_freqs and the size of the ngram should determine some of the wordcloud attributes
# ############################################################################################################################################
# ############################################################################################################################################
# ############################################################################################################################################
# ############################################################################################################################################
# 
# # Define UI for application that draws a histogram
# ui <- fluidPage(
#   # Application title
#   titlePanel("Word Cloud"),
#   
#   sidebarLayout(
#     # Sidebar with a slider and selection inputs
#     sidebarPanel(
#       selectInput(inputId = "selection", label = "Choose a company:",
#                   choices = company),
#       actionButton("update", "Change"),
#       hr(),
#       sliderInput(inputId = "freq",
#                   label = "Minimum Frequency:",
#                   min = 1,  max = 50, value = 15),
#       sliderInput(inputId = "max",
#                   label = "Maximum Number of Words:",
#                   min = 1,  max = 300,  value = 100),
#       numericInput(inputId = "ngram_number", label = "Choose n for n-grams", value = 1,
#                    min = 1, max = 4, step = 1)
#     ),
#     
#     # Show Word Cloud
#     mainPanel(
#       plotOutput("plot")
#     )
#   )
# )
# 
# server <- function(input, output, session) {
#   # Define a reactive expression for the document term matrix
#   terms <- reactive({
#     # Change when the "update" button is pressed...
#     input$update
#     # ...but not for anything else
#     isolate({
#       withProgress({
#         setProgress(message = "Processing corpus...")
#         # return(make_wordcloud(input$selection, input$ngram_number))
#       })
#     })
#   })
# 
#   # Make the wordcloud drawing predictable during a session
#   # par(mar=rep(0, 4))
#   wordcloud_rep <- repeatable(wordcloud)
# 
#   output$plot <- renderPlot({
#     terms()})
#   #   # word_freqs <- data.frame(term = names(term_frequency), num = term_frequency)
#   #   wordcloud_rep(("word_freqs_tw"+input)$term, ("word_freqs_tw"+input)$num, random.order = FALSE,
#   #                 vfont=c("serif","bold"), scale = c(7,1), min.freq=25,
#   #                 max.words=300,colors=brewer.pal(10, "Paired"), rot.per=0, fixed.asp = FALSE)
#   #   # wordcloud_rep(names(v), v, scale=c(4,0.25),
#   #                 # min.freq = input$freq, max.words=input$max,
#   #                 # colors=brewer.pal(8, "Dark2"))
#   # })
# }
# 
# # Run the application 
# shinyApp(ui = ui, server = server)
# ############################################################################################################################################
# ############################################################################################################################################
# ############################################################################################################################################
# ############################################################################################################################################
# 
# 
# 
# ##########tf-idf weighting   taken care of above in the make_wordcloud function
# # cleaned_tweet_corpus <- clean_tweet_corpus(nanopore)
# # tfidf_tdm <- TermDocumentMatrix(cleaned_tweet_corpus,control=list(weighting=weightTfIdf))
# # tfidf_tdm_m <- as.matrix(tfidf_tdm)
# # 
# # # Term Frequency
# # term_frequency <- rowSums(tfidf_tdm_m)
# # # Sort term_frequency in descending order
# # term_frequency <- sort(term_frequency,dec=TRUE)
# # ############Word Cloud
# # # library(wordcloud)
# # # Create word_freqs
# # word_freqs <- data.frame(term = names(term_frequency), num = term_frequency)
# # # Create a wordcloud for the values in word_freqs
# # wordcloud(word_freqs$term, word_freqs$num, min.freq=5,random.order = FALSE, max.words=3000,colors=brewer.pal(8, "Paired"),rot.per=0, fixed.asp = FALSE)
# 
# 
# 
# 
# 
# 
# check <- str_flatten(tw23andMe$text, collapse = " ")
# 
# 
# # # 
# 
# 
# 
# # the text is an column of tweets, flatten it to make when text object 
# 
# # make sure to iconv the file$text to ASCII first
# tweets <- c(str_flatten(tw23andMe$text, collapse = " "), str_flatten(Ancestry$text, collapse = " "))
# 
# tweets <- c(str_flatten(tw23andMe$text, collapse = " "), str_flatten(illumina$text, collapse = " "),
#             str_flatten(nanopore$text, collapse = " "), str_flatten(SOPHiAGENETICS$text, collapse = " "),
#             str_flatten(VeritasGenetics$text, collapse = " "), str_flatten(Ancestry$text, collapse = " "))
# 
# # clean_tweet_corpus(tweets)
# #making a corpus of a vector source
# # Create a corpus
# tweet_corpus <- VCorpus(VectorSource(tweets))
# 
# 
# 
# #apply the clean_corpus function on your review_corpus
# cleaned_tweet_corpus <- clean_corpus(tweet_corpus)
# 
# ########### TDM########
# TDM_speech <- TermDocumentMatrix(cleaned_tweet_corpus)
# 
# TDM_speech_m <- as.matrix(TDM_speech)
# 
# #############################Commonality Cloud
# commonality.cloud(TDM_speech_m,colors=brewer.pal(8, "Dark2"),max.words = 200, random.order=FALSE,rot.per=0, fixed.asp = FALSE)
# 
# #################Comparison Cloud
# TDM_2tweeters <- TermDocumentMatrix(cleaned_tweet_corpus)
# colnames(TDM_2tweeters) <- c("23andMe","illumina","Oxford Nanopore", "SOPHiAGENETICS", "Veritas Genetics", "Ancestry")
# TDM_2tweeters_m <- as.matrix(TDM_2tweeters)
# comparison.cloud(TDM_2tweeters_m,colors=brewer.pal(8, "Dark2"),max.words = 200, title.size = 1, scale = c(3, .3), match.colors = TRUE,random.order=FALSE,rot.per=0, fixed.asp = FALSE)
# 
# ###########################################################################################################
# ###########################################################################################################
# ###########################################################################################################
# ###########################################################################################################
# ###########################################################################################################
# ###########################################################################################################
# ### SENTIMENT ANALYSIS
# tw23andMe$text <- iconv(tw23andMe$text, from = "UTF-8", to = "ASCII", sub = "")
# illumina$text <- iconv(illumina$text, from = "UTF-8", to = "ASCII", sub = "")
# nanopore$text <- iconv(nanopore$text, from = "UTF-8", to = "ASCII", sub = "")
# SOPHiAGENETICS$text <- iconv(SOPHiAGENETICS$text, from = "UTF-8", to = "ASCII", sub = "")
# VeritasGenetics$text <- iconv(VeritasGenetics$text, from = "UTF-8", to = "ASCII", sub = "")
# Ancestry$text <- iconv(Ancestry$text, from = "UTF-8", to = "ASCII", sub = "")
# 
# tweet_corpus <- VCorpus(VectorSource(tw23andMe$text))
# 
# cleaned_tweet_corpus <- clean_corpus(tweet_corpus)
# 
# 
# tidy_mytext <- tidy(TermDocumentMatrix(cleaned_tweet_corpus))
# bing_lex <- get_sentiments("bing")
# mytext_bing <- inner_join(tidy_mytext, bing_lex, by = c("term" = "word"))
# mytext_bing$sentiment_n <- ifelse(mytext_bing$sentiment=="negative", -1, 1)
# mytext_bing$sentiment_score <- mytext_bing$count*mytext_bing$sentiment_n
# 
# sentiment_summary <- mytext_bing %>%
#   group_by(document) %>%
#   summarize(review_sentiment = sum(sentiment_score)) %>%
#   arrange(desc(review_sentiment))
# 
# hist(sentiment_summary$review_sentiment)
# ########## Emotions #############
# 
# nrc_lex <- get_sentiments("nrc")
# table(nrc_lex$sentiment)
# 
# loughran_lex <- get_sentiments("loughran")
# table(loughran_lex$sentiment)
# 
# ###########################################################################################################
# # Sentiment using sentimentrsentimentr
# 
# 
# make_sentiment <- function(tweet_text){
#   sentiment(tweet_text, polarity_dt = lexicon::hash_sentiment_jockers_rinker,
#           valence_shifters_dt = lexicon::hash_valence_shifters,
#           amplifier.weight = 0.8, n.before = 5, n.after = 2,
#           question.weight = 1, adversative.weight = 0.85)
#   text_sentiment <- extract_sentiment_terms(tweet_text)
#   
#   out <- sentiment_by(tweet_text)
#   highlight(out)
# 
#   return(text_sentiment)
#   }
# 
# tw23andMe_sentiment <- make_sentiment(tw23andMe$text)
# illumina_sentiment <- make_sentiment(illumina$text)
# nanopore_sentiment <- make_sentiment(nanopore$text)
# SOPHiAGENETICS_sentiment <- make_sentiment(SOPHiAGENETICS$text)
# VeritasGenetics_sentiment <- make_sentiment(VeritasGenetics$text)
# Ancestry_sentiment <- make_sentiment(Ancestry$text)
# 
# check <- str_flatten(tw23andMe$text, collapse = " ")
# check2 <- str_flatten(cleaned_tweet_corpus, collapse = " ")
# 
# make_sentiment(check)
# ###########################################################################################################
# ###########################################################################################################
# ###########################################################################################################
# 
# # EMOTIONS
# # mytext <- read.csv("app_reviews.csv")
# # mytext$Review <- iconv(mytext$Review, from = "UTF-8", to = "ASCII", sub = "")
# 
# 
# clean_corpus <- function(corpus){
#   cleaned_corpus <- tm_map(corpus, content_transformer(replace_abbreviation))
#   cleaned_corpus <- tm_map(corpus, content_transformer(replace_ordinal))
#   cleaned_corpus <- tm_map(cleaned_corpus, content_transformer(tolower))
#   cleaned_corpus <- tm_map(cleaned_corpus, removePunctuation)
#   cleaned_corpus <- tm_map(cleaned_corpus, removeNumbers)
#   # cleaned_corpus <- tm_map(cleaned_corpus, stemDocument)
#   # cleaned_corpus <- tm_map(cleaned_corpus, stemCompletion(dictionary = cleaned_corpus, type='prevalent'))
#   custom_stop_words <- c("amp", "can", "may", "via", "one", "httpstcodggafu", "youre", "nanopore", "nanoporeconf","illumina", "andme", "andmes", "sophia", "veritas", "veritasgenetics")
#   cleaned_corpus <- tm_map(cleaned_corpus, removeWords, stopwords(kind="en"))
#   cleaned_corpus <- tm_map(cleaned_corpus, removeWords, custom_stop_words)
#   custom_stop_words2 <- c("nanopore", "illumina", "andme", "sophia", "sophiagenetics", "veritas")
#   cleaned_corpus <- tm_map(cleaned_corpus, removeWords, custom_stop_words2)
#   cleaned_corpus <- tm_map(cleaned_corpus, stripWhitespace)
#   return(cleaned_corpus)
# }
# ##############################################3
# 
# # tweet_corpus <- VCorpus(VectorSource(tw23andMe$text))
# # 
# # cleaned_tweet_corpus <- clean_corpus(tweet_corpus)
# 
# make_emo_radar <- function(tweets){
#   tweet_corpus_emo <- VCorpus(VectorSource(tweets$text))
#   cleaned_tweet_corpus_emo <- clean_corpus(tweet_corpus_emo)
# 
# 
#   tidy_tweet_emo <- tidy(TermDocumentMatrix(cleaned_tweet_corpus_emo))
#   bing_lex <- get_sentiments("bing")
#   mytext_bing <- inner_join(tidy_tweet_emo, bing_lex, by = c("term" = "word"))
#   mytext_bing$sentiment_n <- ifelse(mytext_bing$sentiment=="negative", -1, 1)
#   mytext_bing$sentiment_score <- mytext_bing$count*mytext_bing$sentiment_n
#   
#   sentiment_summary <- mytext_bing %>%
#     group_by(document) %>%
#     summarize(review_sentiment = sum(sentiment_score)) %>%
#     arrange(desc(review_sentiment))
#   
#   
#   nrc_lex <- get_sentiments("nrc")
#   mytext_nrc <- inner_join(tidy_tweet_emo, nrc_lex, by = c("term" = "word"))
#   # remove positive and negative first
#   mytext_nrc_noposneg <- mytext_nrc[!(mytext_nrc$sentiment %in% c("positive","negative")),]
#   emotion_summary <- mytext_nrc_noposneg %>%
#     group_by(document,sentiment) %>%
#     summarize(review_sentiment = sum(count)) %>%
#     arrange(desc(review_sentiment))
#   
#   emotion_overall_summary <- mytext_nrc_noposneg %>%
#     group_by(sentiment) %>%
#     summarize(review_sentiment = sum(count)) %>%
#     arrange(desc(review_sentiment))
#   
#   chartJSRadar(emotion_overall_summary)
# }
# 
# make_emo_radar(tw23andMe)
# make_emo_radar(illumina)
# make_emo_radar(nanopore)
# make_emo_radar(SOPHiAGENETICS)
# make_emo_radar(VeritasGenetics)
# make_emo_radar(Ancestry)
# 
# 
# # TODO
# # Add lexicon argument to emo function (instead of just bing)
# # sentiment function needs to be cleaned up, hopefully using the same as gene_tweets_clean (because I can't make it a corpus)
# # PowerBI or Shiny implementation
# 
# 
# 
# 
# 
# ##backup before making the function
# # tweet_corpus_emo <- VCorpus(VectorSource(tw23andMe$text))
# # cleaned_tweet_corpus_emo <- clean_corpus(tweet_corpus_emo)
# # 
# # 
# # tidy_tweet_emo <- tidy(TermDocumentMatrix(cleaned_tweet_corpus_emo))
# # bing_lex <- get_sentiments("bing")
# # mytext_bing <- inner_join(tidy_tweet_emo, bing_lex, by = c("term" = "word"))
# # mytext_bing$sentiment_n <- ifelse(mytext_bing$sentiment=="negative", -1, 1)
# # mytext_bing$sentiment_score <- mytext_bing$count*mytext_bing$sentiment_n
# # 
# # sentiment_summary <- mytext_bing %>%
# #   group_by(document) %>%
# #   summarize(review_sentiment = sum(sentiment_score)) %>%
# #   arrange(desc(review_sentiment))
# # 
# # 
# # nrc_lex <- get_sentiments("nrc")
# # mytext_nrc <- inner_join(tidy_tweet_emo, nrc_lex, by = c("term" = "word"))
# # # remove positive and negative first
# # mytext_nrc_noposneg <- mytext_nrc[!(mytext_nrc$sentiment %in% c("positive","negative")),]
# # emotion_summary <- mytext_nrc_noposneg %>%
# #   group_by(document,sentiment) %>%
# #   summarize(review_sentiment = sum(count)) %>%
# #   arrange(desc(review_sentiment))
# # 
# # emotion_overall_summary <- mytext_nrc_noposneg %>%
# #   group_by(sentiment) %>%
# #   summarize(review_sentiment = sum(count)) %>%
# #   arrange(desc(review_sentiment))
# # 
# # chartJSRadar(emotion_overall_summary)
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###########################################################################################################
# ###########################################################################################################
# ###########################################################################################################
# 
# #### Bigrams and trigrams will be too large a vector, so they would have to be done on 
# #### an individual company basis 
# 
# ##############bigrams and trigrams
# #create a structure for bi-gram 
# library(RWeka)
# tokenizer <- function(x)
#   NGramTokenizer(x,Weka_control(min=2,max=2))
# bigram_tdm <- TermDocumentMatrix(cleaned_tweet_corpus,control = list(tokenize=tokenizer))
# bigram_tdm_m <- as.matrix(bigram_tdm)
# 
# # Term Frequency
# term_frequency <- rowSums(bigram_tdm_m)
# # Sort term_frequency in descending order
# term_frequency <- sort(term_frequency,dec=TRUE)
# ############Word Cloud
# # Create word_freqs
# word_freqs <- data.frame(term = names(term_frequency), num = term_frequency)
# # Create a wordcloud for the values in word_freqs
# wordcloud(word_freqs$term, word_freqs$num,min.freq=5,max.words=500,colors=brewer.pal(8, "Paired"))
# 
# 
# 
# #create a structure for tri-gram 
# tokenizer <- function(x)
#   NGramTokenizer(x,Weka_control(min= 4,max=4))
# trigram_tdm <- TermDocumentMatrix(cleaned_tweet_corpus,control = list(tokenize=tokenizer))
# trigram_tdm_m <- as.matrix(trigram_tdm)
# 
# # Term Frequency
# term_frequency <- rowSums(trigram_tdm_m)
# # Sort term_frequency in descending order
# term_frequency <- sort(term_frequency,dec=TRUE)
# ############Word Cloud
# # Create word_freqs
# word_freqs <- data.frame(term = names(term_frequency), num = term_frequency)
# # Create a wordcloud for the values in word_freqs
# wordcloud(word_freqs$term, word_freqs$num,min.freq=5,max.words=500,colors=brewer.pal(8, "Paired"))
# 
# 
# 
# # columns <- c('created_at',	'screen_name',	'text',	'display_text_width',	'reply_to_screen_name',	'is_quote',	'is_retweet',
# #              'favorite_count',	'retweet_count',	'hashtags',	'urls_expanded_url',	'ext_media_expanded_url',	'ext_media_type',
# #              'mentions_screen_name',	'lang',	'quoted_text',	'quoted_created_at',	'quoted_favorite_count',	'quoted_retweet_count',
# #              'quoted_screen_name',	'quoted_name',	'quoted_followers_count',	'quoted_friends_count',	'quoted_statuses_count',
# #              'quoted_location',	'quoted_description',	'quoted_verified',	'retweet_text',	'retweet_created_at',	'retweet_favorite_count',
# #              'retweet_retweet_count',	'retweet_screen_name',	'retweet_name',	'retweet_followers_count',	'retweet_location',
# #              'retweet_description',	'retweet_verified',	'place_full_name',	'country',	'status_url')
# # 
# # # The csv has a lot of columns that aren't needed
# # gene_tweets <- gene_tweets %>%
# #   select(columns)
# 
# 
# # ### Count and plot tweet terms, minus stop words with at least 2 letters in the word, without company names
# # tweet_term_count <- freq_terms(gene_tweets_clean$text, stopwords=c(Top200Words, 'httpstco', 'amp', 'via', 'hi', 'sophia', 'nanopore',
# #                         'andme', 'ancestry', 'illumina', 'veritas'), at.least =2, top = 20)
# # plot(tweet_term_count)
# # 
# # ### Count and plot tweet terms, minus stop words with at least 2 letters in the word, including company names
# # tweet_term_count <- freq_terms(gene_tweets_clean$text, stopwords=c(Top200Words, 'httpstco', 'amp','via','hi' ), at.least =2, top = 20)
# # plot(tweet_term_count)
# 
#   #############################Commonality Cloud
# wordcloud(TDM_gene_m,colors=brewer.pal(8, "Dark2"),max.words = 100, random.order=FALSE)
# 
# #################Comparison Cloud
# TDM_speech <- TermDocumentMatrix(cleaned_tweet_corpus)
# TDM_speech_m <- as.matrix(TDM_speech)
# comparison.cloud(TDM_gene_m,colors=brewer.pal(8, "Dark2"),max.words = 200, random.order=FALSE)
# 
# #####
# #####
# #####
# #####
# ##############bigrams and trigrams
# #making a corpus of a vector source 
# obama_speech_corpus <- VCorpus(VectorSource(obama))
# 
# cleaned_obama_speech_corpus <- clean_corpus(obama_speech_corpus)
# 
# #create a structure for bi-gram 
# tokenizer <- function(x)
#   NGramTokenizer(x,Weka_control(min=2,max=2))
# bigram_tdm <- TermDocumentMatrix(cleaned_obama_speech_corpus,control = list(tokenize=tokenizer))
# bigram_tdm_m <- as.matrix(bigram_tdm)
# 
# # Term Frequency
# term_frequency <- rowSums(bigram_tdm_m)
# # Sort term_frequency in descending order
# term_frequency <- sort(term_frequency,dec=TRUE)
# ############Word Cloud
# # Create word_freqs
# word_freqs <- data.frame(term = names(term_frequency), num = term_frequency)
# # Create a wordcloud for the values in word_freqs
# wordcloud(word_freqs$term, word_freqs$num,min.freq=5,max.words=500,colors=brewer.pal(8, "Paired"))
# 
# #######
# ##############bigrams and trigrams
# #making a corpus of a vector source 
# trump_speech_corpus <- VCorpus(VectorSource(obama))
# 
# cleaned_trump_speech_corpus <- clean_corpus(obama_speech_corpus)
# 
# #create a structure for bi-gram 
# tokenizer <- function(x)
#   NGramTokenizer(x,Weka_control(min=2,max=2))
# bigram_tdm <- TermDocumentMatrix(cleaned_trump_speech_corpus,control = list(tokenize=tokenizer))
# bigram_tdm_m <- as.matrix(bigram_tdm)
# 
# # Term Frequency
# term_frequency <- rowSums(bigram_tdm_m)
# # Sort term_frequency in descending order
# term_frequency <- sort(term_frequency,dec=TRUE)
# ############Word Cloud
# # Create word_freqs
# word_freqs <- data.frame(term = names(term_frequency), num = term_frequency)
# # Create a wordcloud for the values in word_freqs
# wordcloud(word_freqs$term, word_freqs$num,min.freq=5,max.words=500,colors=brewer.pal(8, "Paired"))
# 
# #########################
