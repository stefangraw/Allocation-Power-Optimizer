library(shiny)
library(shinyjs)
library(parallel)

server = function(input,output){
  # dynamic min for totalSampleSize
  # output$totalSampleSize <- renderUI({
  #   numericInput(inputId = "totalSampleSize", label = "Total sample size", value = NA, step = 1, min = input$numberOfArms*4)
  # })
  
  
  ## run program here
  source("apo.R")
  
  runProgram <- eventReactive(input$goButton, {
    # hide("goButton")
    randomSeed = sample(1:1000000, 1)
    out = list()
    runTimeStart = Sys.time()
    withProgress(message = 'Working...', detail = "", value = NULL, {
      sobj <- optimize.power(
        target_alpha = input$alpha,
        arm_mus = as.numeric(
          if(input$effectTrend == "oneBestArm"){
            c(rep(0,input$numberOfArms-1), ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom))
          } else if(input$effectTrend == "linear"){
            create.linear.means(num_arms = input$numberOfArms, best_arm_mu = as.numeric(ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom)))
          } else if(input$effectTrend == "custom"){
            c(0,as.numeric(unlist(strsplit(input$customMus,","))))
          }), 
        n_patients_TOTAL = input$totalSampleSize,
        post_sample_size = input$postSampleSize,
        simulation_runs_H0 = input$simsH0,
        simulation_runs_H1 = input$simsH1, 
        threads = input$threads,
        num_checked_sizes = input$checkSampleSize, 
        seed = ifelse(input$randomSeed == FALSE, input$seed, randomSeed),
        delta = input$thresholdFutility_delta,
        epsilon = input$thresholdFutility_epsilon
      )
    })
    runTimeStop = difftime(Sys.time(), runTimeStart, units = "auto")
    output$powerPlot <- renderPlot({plot.sim(sobj)})
    
    #Bonf-adjusted t-test (this needs to use all arms)
    if(input$effectTrend == "oneBestArm"){ 
      mus <- c(rep(0, (input$numberOfArms-1)), as.numeric(ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom)))
    } else if(input$effectTrend == "linear") {
      mus <- seq(0,ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom), length.out=input$numberOfArms)
    } else if(input$effectTrend == "custom") {
      mus <- c(0,as.numeric(unlist(strsplit(input$customMus,","))))
    }
    freqTest_single <- c()
    for(mus_indx in 2:length(mus))
    {
      freqTest_single[mus_indx-1] <- power.t.test(n = floor(input$totalSampleSize/length(mus)), delta = mus[mus_indx], alternative = "one.sided", sig.level = input$alpha/(length(mus)-1))[["power"]]
    }
    freqTest_bonf <- 1 - prod(1-freqTest_single)
    
    freqTest_bestArmKnow = power.t.test(n = floor(input$totalSampleSize/2), delta = max(mus), alternative = "one.sided", sig.level = input$alpha)
    
    table1 = data.frame(
      "Power" = c(round(sobj$results.best$power,3),  round(freqTest_bonf,3), round(freqTest_bestArmKnow$power,3)), 
      "N per arm (Stage 1)" = c(sobj$results.best$n.per.arm.I, " ", " "),
      "Additional N per arm (Stage 2)" = c(as.integer(sobj$results.best$n.per.arm.II), " ", " "), 
      "Total N" = as.integer(c(sobj$results.best$n.total, floor(input$totalSampleSize/input$numberOfArms)*input$numberOfArms, floor(input$totalSampleSize/2)*2)), 
      "Threshold" = c(round(sobj$results.best$threshold,3), " ", " "), 
      "Pr(Early stop: futility)" = c(round(sobj$results.best$pr.futility,3), " ", " "),
      "Pr(success, best arm wins)" = c(round(sobj$results.best$pr.best.tx.success,3)," ", " "),
      check.names = F)
    rownames(table1) = c("Adaptive design", "Bonferroni adjusted t", "If best arm known")
    
    logString = paste0(
      "Critical value alpha = ", input$alpha, "\n",
      "Total sample size = ", input$totalSampleSize, "\n",
      
      if(input$effectTrend == "oneBestArm"){
        paste0("Effect size (one best arm) = ", ifelse(input$effectSize == 0, input$effectSizeCustom, input$effectSize), "\n")
      } else if(input$effectTrend == "linear"){
        paste0("Effect size (linear trend) = ", ifelse(input$effectSize == 0, input$effectSizeCustom, input$effectSize), "\n")
      }, #else if(input$effectTrend == "custom"){
      #   paste0("Custom = ", "coming XXXXX", "\n")
      # },
      
      "H1 treatment means = ", paste(round(as.numeric(
        if(input$effectTrend == "oneBestArm"){
          c(rep(0,input$numberOfArms-1), ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom))
        } else if(input$effectTrend == "linear"){
          create.linear.means(num_arms = input$numberOfArms, best_arm_mu = as.numeric(ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom)))
        } else if(input$effectTrend == "custom"){
          c(0,as.numeric(unlist(strsplit(input$customMus,","))))
        }),2), collapse = ", "), "\n",
      
      "Post sample size = ", input$postSampleSize, "\n",
      "Simulation runs H0 = ", input$simsH0, "\n",
      "Simulation runs H1 = ", input$simsH1, "\n",
      "Number of checked sample sizes  = ", input$checkSampleSize, "\n",
      "Futility-Delta  = ", input$thresholdFutility_delta, "\n",
      "Futility-Epsilon  = ", input$thresholdFutility_epsilon, "\n",
      "Number of threads = ", input$threads, "\n",
      "Seed = ", ifelse(input$advancedSettings == 1 & input$randomSeed == FALSE, input$seed, randomSeed), "\n",
      "Run time = ", round(runTimeStop,1), " ", attr(runTimeStop, "units"),  
      sep = "")
    
    output$log <- renderText({logString})
    
    out[["sobj"]] = sobj
    out[["table1"]] = table1
    out
  }) # processbar done
  
  output$bestN <- renderTable({runProgram()[["table1"]]}, rownames = T, digits = 3)
  
  output$powerTable <- renderDataTable({data.frame(
    Power = round(runProgram()[["sobj"]]$results.all$power,4), 
    "N per arm (Stage 1)" = runProgram()[["sobj"]]$results.all$n.per.arm.I, 
    "Additional N per arm (Stage 2)" = runProgram()[["sobj"]]$results.all$n.per.arm.II, 
    "Total N" = runProgram()[["sobj"]]$results.all$n.total, 
    "Threshold" = round(runProgram()[["sobj"]]$results.all$threshold,5),
    "Pr(Early stop: futility)" = round(runProgram()[["sobj"]]$results.all$pr.futility,5),
    "Pr(success, best arm wins)" = round(runProgram()[["sobj"]]$results.all$pr.best.tx.success,5),
    check.names = F)
  })
  # hide("goButton")
}














ui = fluidPage(
  tags$head(tags$style(HTML(".shiny-notification {
              height: 100px;
              width: 600px;
              position:fixed;
              font-size: 500%;
              text-align: center;
              line-height: 75px;
              top: calc(50% - 50px);;
              left: calc(50% - 300px);;}"))),
  
  h1("Allocation & Power Optimizer"),
  h3("Adaptive Bayes Two-Stage Drop the Losers Design"),
  h5("Authors: Alex Karanevich, Richard Meier, Stefan Graw"),
  h5("Department of Biostatistics, University of Kansas School of Medicine"),
  h5("View the ",a("Instructions",target="_blank",href="instructions.pdf")),
  sidebarLayout( 
    sidebarPanel( # Box 2
      
      # show sample size input
      # Box 2 | sample size
      numericInput(inputId = "alpha", label = HTML("Critical value &alpha;"), value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput(inputId = "totalSampleSize", label = "Total sample size", value = 100, step = 1, min = 10), 
      # uiOutput("totalSampleSize"), # dynamic min calculated
      
      radioButtons(inputId = "effectTrend", label = "Effect size as", choices = list("Best arm" = "oneBestArm", "Linear trend" = "linear", "Custom" = "custom")),
      conditionalPanel(
        condition = "input.effectTrend == 'oneBestArm' || input.effectTrend == 'linear'",
        numericInput(inputId = "numberOfArms", label = "Number of arms (including control)", value = 2, step = 1, min = 2),
        radioButtons(inputId = "effectSize", label = "Standardized effect size", choices = list("Small (d=0.2)" = 0.2, "Medium (d=0.5)" = 0.5, "Large (d=0.8)" = 0.8, "Custom" = 0)),
        # CARE effect size has 2 variables now XXXXX
        conditionalPanel(
          condition = "input.effectSize == 0",
          numericInput(inputId = "effectSizeCustom", label = "Standardized effect size", value = 0.5, min = 0, max = 2, step = 0.1) # https://en.wikipedia.org/wiki/Effect_size#Cohen's_d
        )
      ),
      conditionalPanel(
        condition = "input.effectTrend == 'custom'",
        textInput(inputId = "customMus", label = HTML("H<sub>1</sub> treatment means </br> (comma separated; control excluded)") , value = "0.2, 0.5")
      ),
      
      numericInput(inputId = "threads", label = "Number of cores/threads", value = round(detectCores(all.tests = FALSE, logical = TRUE)/2), min = 1, max = detectCores(all.tests = FALSE, logical = TRUE), step = 1) # XXXXX set default to 1
      
    ),
    
    # power plot
    mainPanel(
      plotOutput(outputId = "powerPlot"),
      
      # "Best n1 n2 power:",
      tableOutput(outputId = "bestN")
    )
  ),
  
  
  # Box 3
  sidebarLayout(
    sidebarPanel(
      checkboxInput(inputId = "advancedSettings", label = "Advanced settings"),
      conditionalPanel( # Box 2 | power
        condition = "input.advancedSettings == 1",
        HTML("<b>Futility Stage 1: Stop trial early if<br /><div align='center'>P(&mu;<sub>Best</sub> - &mu;<sub>Ctrl</sub> &ge; <font color='red'>&delta;</font>) &le; <font color='red'>&epsilon;</font></div></b>"), 
        wellPanel(
          numericInput(inputId = "thresholdFutility_delta", label = HTML("<font color='red'>&delta;</font> = ?"), value = 0, min = 0, max = 1, step = 0.1),
          numericInput(inputId = "thresholdFutility_epsilon", label = HTML("<font color='red'>&epsilon;</font> = ?"), value = 0, min = 0, max = 1, step = 0.1)
        ),
        
        # numericInput(inputId = "superiority", label = HTML("Minimum treatment improvment:<br /> P(&mu;<sub>Best</sub> - &mu;<sub>Ctrl</sub> > <font color='red'>?</font>) > Th"), value = 0, min = 0), # max = 2*effect size
        #uiOutput("superiority"), 
        numericInput(inputId = "checkSampleSize", label = "Number of checked sample sizes", value = 20, min = 10, step = 1),
        numericInput(inputId = "postSampleSize", label = "Post sample size", value = 1000, min = 100, step = 100), 
        numericInput(inputId = "simsH0", label = HTML("Simulation runs H<sub>0</sub>"), value = 1000, min = 100, step = 100), 
        numericInput(inputId = "simsH1", label = HTML("Simulation runs H<sub>1</sub>"), value = 10000, min = 100, step = 1000),
        checkboxInput(inputId = "randomSeed", label = "Random seed", value = T),
        conditionalPanel(
          condition = "input.randomSeed == 0",
          numericInput(inputId = "seed", label = "Seed", value = 1234, min = 1, step = 1)
        )
      ), 
      actionButton("goButton", "Go!")
    ),
    
    mainPanel(
      # "Log:",
      verbatimTextOutput(outputId = "log"), 
      
      
      checkboxInput(inputId = "showTable", label = "Show table"),
      conditionalPanel(
        condition = "input.showTable == 1",
        # "tabular:",
        dataTableOutput(outputId = "powerTable")
      )
    )
  )  
)

shinyApp(ui = ui, server = server)
