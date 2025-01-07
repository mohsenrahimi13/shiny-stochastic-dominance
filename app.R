library(shiny)
library(ggplot2)
library(dplyr)
library(boot)

# Define UI and Server in one file for ShinyApps.io compatibility
ui <- fluidPage(
  titlePanel("Stochastic Dominance Dashboard"),
  sidebarLayout(
    sidebarPanel(
      h4("Settings"),
      fileInput("dataFile", "Upload Data (CSV)", accept = ".csv"),
      selectInput("distType", "Distribution Type:", choices = c("Normal", "Uniform", "Custom")),
      numericInput("sampleSize", "Sample Size:", 100, min = 10, max = 1000),
      actionButton("runSim", "Run Simulation"),
      h5("Distortion Function"),
      selectInput("distortionType", "Choose Distortion:", choices = c("Linear", "Power", "Exponential")),
      sliderInput("distortionParam", "Distortion Intensity (k):", min = 0.1, max = 5, value = 1, step = 0.1),
      h5("Advanced Options"),
      selectInput("sampleFramework", "Sampling Framework:", choices = c("Independent", "Matched Pairs")),
      sliderInput("asdDelta", "ASD Delta Threshold:", min = 0, max = 0.1, value = 0.01, step = 0.005)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("SD Test Comparison",
                 plotOutput("ecdfPlot"),
                 plotOutput("ecdfPlotWithBands"),
                 tableOutput("testResults"),
                 verbatimTextOutput("summaryTextComparison")
        ),
        tabPanel("Distorted SD Visualization",
                 plotOutput("dsdPlot"),
                 tableOutput("comparisonResults"),
                 verbatimTextOutput("summaryTextDistorted")
        ),
        tabPanel("Almost SD Test",
                 verbatimTextOutput("asdTestResult")
        ),
        tabPanel("Quantile-Based Test",
                 verbatimTextOutput("quantileTest")
        )
      )
    )
  )
)

server <- function(input, output) {
  dataInput <- reactive({
    if (is.null(input$dataFile)) {
      if (input$distType == "Normal") {
        sample1 <- rnorm(input$sampleSize)
        sample2 <- rnorm(input$sampleSize, mean = 0.05)
      } else if (input$distType == "Uniform") {
        sample1 <- runif(input$sampleSize)
        sample2 <- runif(input$sampleSize, min = 0.45, max = 1.05)
      }
    } else {
      uploadedData <- read.csv(input$dataFile$datapath)
      sample1 <- uploadedData[, 1]
      sample2 <- if (ncol(uploadedData) >= 2) uploadedData[, 2] else sample1
    }
    if (length(sample1) == 0 || length(sample2) == 0) {
      showNotification("Error: No valid data available after loading. Please check the input.", type = "error")
      stop("Invalid data")
    }
    list(sample1 = sample1, sample2 = sample2)
  })
  
  distortionFunction <- reactive({
    k <- input$distortionParam
    switch(
      input$distortionType,
      "Linear" = function(x) x * k,
      "Power" = function(x) ifelse(x >= 0, x^k, NA),
      "Exponential" = function(x) exp(k * x)
    )
  })
  
  output$summaryTextComparison <- renderText({
    data <- dataInput()
    req(data)
    ksTest <- ks.test(data$sample1, data$sample2)
    wmwTest <- wilcox.test(data$sample1, data$sample2)
    
    paste0(
      "Kolmogorov-Smirnov Test: P-value = ", round(ksTest$p.value, 3), ". Result: ",
      ifelse(ksTest$p.value < 0.05, "Reject Null Hypothesis", "Fail to Reject Null Hypothesis"), "\n",
      "Wilcoxon-Mann-Whitney Test: P-value = ", round(wmwTest$p.value, 3), ". Result: ",
      ifelse(wmwTest$p.value < 0.05, "Reject Null Hypothesis", "Fail to Reject Null Hypothesis")
    )
  })
  
  output$summaryTextDistorted <- renderText({
    data <- dataInput()
    req(data)
    sample1 <- data$sample1
    sample2 <- data$sample2
    
    distortedSample2 <- distortionFunction()(sample2)
    distortedSample2 <- distortedSample2[!is.na(distortedSample2) & is.finite(distortedSample2)]
    
    if (length(distortedSample2) == 0) {
      showNotification("No valid distorted sample values for plotting.", type = "error")
      stop("No valid distorted values")
    }
    
    ksTestDistorted <- ks.test(sample1, distortedSample2)
    wmwTestDistorted <- wilcox.test(sample1, distortedSample2)
    
    paste0(
      "Distorted Kolmogorov-Smirnov Test: P-value = ", round(ksTestDistorted$p.value, 3), ". Result: ",
      ifelse(ksTestDistorted$p.value < 0.05, "Reject Null Hypothesis", "Fail to Reject Null Hypothesis"), "\n",
      "Distorted Wilcoxon-Mann-Whitney Test: P-value = ", round(wmwTestDistorted$p.value, 3), ". Result: ",
      ifelse(wmwTestDistorted$p.value < 0.05, "Reject Null Hypothesis", "Fail to Reject Null Hypothesis")
    )
  })
  
  output$testResults <- renderTable({
    data <- dataInput()
    req(data)
    sample1 <- data$sample1
    sample2 <- data$sample2
    
    ksTest <- ks.test(sample1, sample2)
    wmwTest <- wilcox.test(sample1, sample2)
    
    data.frame(
      Test = c("Kolmogorov-Smirnov", "Wilcoxon-Mann-Whitney"),
      Statistic = c(ksTest$statistic, wmwTest$statistic),
      `P-value` = c(ksTest$p.value, wmwTest$p.value),
      `Decision (α=0.05)` = ifelse(c(ksTest$p.value, wmwTest$p.value) < 0.05, "Reject H0", "Fail to Reject H0")
    )
  })
  
  output$ecdfPlot <- renderPlot({
    data <- dataInput()
    req(data)
    sample1 <- data$sample1
    sample2 <- data$sample2
    
    df <- data.frame(
      value = c(sample1, sample2),
      sample = rep(c("Sample 1", "Sample 2"), each = length(sample1))
    )
    
    ggplot(df, aes(x = value, color = sample)) +
      stat_ecdf(geom = "step") +
      ggtitle("Empirical CDF Comparison") +
      theme_minimal()
  })
  
  output$dsdPlot <- renderPlot({
    data <- dataInput()
    req(data)
    sample1 <- data$sample1
    sample2 <- data$sample2
    
    distortedSample2 <- distortionFunction()(sample2)
    distortedSample2 <- distortedSample2[!is.na(distortedSample2) & is.finite(distortedSample2)]
    
    if (length(distortedSample2) == 0) {
      showNotification("No valid distorted sample points to plot", type = "error")
      return()
    }
    
    df <- data.frame(
      value = c(sample1, distortedSample2),
      sample = rep(c("Sample 1", "Distorted Sample 2"), each = length(sample1))
    )
    
    ggplot(df, aes(x = value, color = sample)) +
      stat_ecdf(geom = "step") +
      ggtitle("Empirical CDF After Distortion") +
      theme_minimal()
  })
  
  output$ecdfPlotWithBands <- renderPlot({
    data <- dataInput()
    req(data)
    sample1 <- data$sample1
    ecdf1 <- ecdf(sample1)
    nBoot <- 1000
    
    bootBands <- replicate(nBoot, {
      bootSample <- sample(sample1, replace = TRUE)
      ecdf(bootSample)(sample1)
    })
    
    lowerBand <- apply(bootBands, 1, function(x) quantile(x, 0.025))
    upperBand <- apply(bootBands, 1, function(x) quantile(x, 0.975))
    
    ggplot() +
      geom_line(aes(x = sample1, y = ecdf1(sample1)), color = "blue") +
      geom_ribbon(aes(x = sample1, ymin = lowerBand, ymax = upperBand), alpha = 0.2, fill = "blue") +
      ggtitle("Empirical CDF with Bootstrap Confidence Bands") +
      theme_minimal()
  })
  
  output$asdTestResult <- renderText({
    data <- dataInput()
    req(data)
    delta <- input$asdDelta
    ecdf1 <- ecdf(data$sample1)
    ecdf2 <- ecdf(data$sample2)
    
    violations <- mean(ecdf1(data$sample2) - ecdf2(data$sample2) > delta)
    
    if (violations < 0.05) {
      "Almost stochastic dominance holds with small violations."
    } else {
      "Almost stochastic dominance does not hold."
    }
  })
  
  output$quantileTest <- renderText({
    data <- dataInput()
    req(data)
    q1 <- quantile(data$sample1, probs = seq(0.1, 0.9, 0.1))
    q2 <- quantile(data$sample2, probs = seq(0.1, 0.9, 0.1))
    
    if (all(q1 >= q2)) {
      "Sample 1 stochastically dominates Sample 2 in quantiles."
    } else {
      "No stochastic dominance in quantiles."
    }
  })
}

# Run App
shinyApp(ui = ui, server = server)
