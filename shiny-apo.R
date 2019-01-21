library(shiny)
library(shinyjs)
library(parallel)
library(doParallel)

server = function(input,output){
  
  source("apo.R")
  
  ## run program here
  
  runProgram1 <- eventReactive(input$goButton, {
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
        epsilon = input$thresholdFutility_epsilon,
		ci_sig_level = 1-input$simsCI
      )
      
      output$alphaCI <- renderText({HTML(paste0(
        "\nFor threshold of ",round(sobj$results.best$threshold,3)," the ","<b>"," simulated type I error rate (",input$simsCI * 100, "% CI)","</b>", 
        " is ", "&alpha;" ," = ", input$alpha, " (", round(sobj$results.best$alpha.low, 3),
        ", ", round(sobj$results.best$alpha.high, 3), ").\n", sep = ""
      ))})
      
      output$powerPlot <- renderPlot({plot.sim(sobj)})
      
      #Bonf-adjusted t-test (this needs to use all arms)
      if(input$effectTrend == "oneBestArm"){ 
        mus <- c(rep(0, (input$numberOfArms-1)), as.numeric(ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom)))
      } else if(input$effectTrend == "linear") {
        mus <- seq(0,ifelse(input$effectSize != 0, input$effectSize, input$effectSizeCustom), length.out=input$numberOfArms)
      } else if(input$effectTrend == "custom") {
        mus <- c(0,as.numeric(unlist(strsplit(input$customMus,","))))
      }
      
      freqTest_bonf <- bonf.power.parallel(n=floor(input$totalSampleSize/length(mus)),means=mus,alpha=input$alpha,simulations=input$simsFreq,threads=input$threads, seed=ifelse(input$randomSeed == FALSE, input$seed, randomSeed))
	  boniv = clopper.pearson.interval(size=input$simsFreq,prob=freqTest_bonf$pwr.raw,alpha=1-input$simsCI)
	  freqTest_bonf$power.low <- boniv[1]
	  freqTest_bonf$power.high <- boniv[2]
	  
	  freqTest_dunnet <- dunnett.power.parallel(n=floor(input$totalSampleSize/length(mus)),means=mus,alpha=input$alpha,simulations=input$simsFreq,threads=input$threads, seed=ifelse(input$randomSeed == FALSE, input$seed, randomSeed))
	  duniv = clopper.pearson.interval(size=input$simsFreq,prob=freqTest_dunnet$pwr.raw,alpha=1-input$simsCI)
	  freqTest_dunnet$power.low <- duniv[1]
	  freqTest_dunnet$power.high <- duniv[2]
	  
	  freqTest_bestArmKnow = list(pwr.raw=(power.t.test(n = floor(input$totalSampleSize/2), delta = max(mus), alternative = "one.sided", sig.level = input$alpha))$power)
	  
	  runTimeStop = difftime(Sys.time(), runTimeStart, units = "auto")
	  
    })
    
    table1 = data.frame(
      "thisWillBeChanged" = c( 
		paste(round(sobj$results.best$power,3), " (", round(sobj$results.best$power.low,3), ", ", round(sobj$results.best$power.high,3), ")", sep = ""), 
		paste(round(freqTest_bonf$pwr.raw,3), " (", round(freqTest_bonf$power.low,3), ", ", round(freqTest_bonf$power.high,3), ")", sep = ""), 
		paste(round(freqTest_dunnet$pwr.raw,3), " (", round(freqTest_dunnet$power.low,3), ", ", round(freqTest_dunnet$power.high,3), ")", sep = ""),
		paste(round(freqTest_bestArmKnow$pwr.raw,3), " (exact)", sep = "")
	  ), 
      "N per arm (Stage 1)" = c(sobj$results.best$n.per.arm.I, " ", " ", " "),
      "Additional N per arm (Stage 2)" = c(as.integer(sobj$results.best$n.per.arm.II), " ", " ", " "), 
      "Total N" = as.integer(c(sobj$results.best$n.total, floor(input$totalSampleSize/input$numberOfArms)*input$numberOfArms, floor(input$totalSampleSize/input$numberOfArms)*input$numberOfArms, floor(input$totalSampleSize/2)*2)), 
      "Threshold" = c(round(sobj$results.best$threshold,3), " ", " ", " "), 
      "Pr(Early stop: futility)" = c(round(sobj$results.best$pr.futility,3), " ", "", " "),
      "Pr(success, best arm wins)" = c(round(sobj$results.best$pr.best.tx.success,3)," ", " ", " "),
      check.names = F)
    rownames(table1) = c("Adaptive DTL design", "Bonferroni adjusted t", "Dunnett adjusted t", "If best arm known")
    colnames(table1)[1] = paste("Power (", input$simsCI * 100, "% CI)", sep = "")
    
    
    logString = paste0(
      "Critical value alpha = ", input$alpha, "\n",
      "Total sample size = ", input$totalSampleSize, "\n",
      
      if(input$effectTrend == "oneBestArm"){
        paste0("Effect size (one best arm) = ", ifelse(input$effectSize == 0, input$effectSizeCustom, input$effectSize), "\n")
      } else if(input$effectTrend == "linear"){
        paste0("Effect size (linear trend) = ", ifelse(input$effectSize == 0, input$effectSizeCustom, input$effectSize), "\n")
      }, 
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
	  "Confidence level for simulated alpha and power = ", input$simsCI, "\n",
	  "Simulation runs other = ", input$simsFreq, "\n",
      "Number of n1 and n2 combinations  = ", input$checkSampleSize, "\n",
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
  
  runProgram2 <- eventReactive(input$goButtonB, {
    # hide("goButton")
    randomSeedB = sample(1:1000000, 1)
    out = list()
    runTimeStartB = Sys.time()
    withProgress(message = 'Working...', detail = "", value = NULL, {
      sobj <- optimize.N(
        target_alpha = input$alphaB,
        arm_mus = as.numeric(
          if(input$effectTrendB == "oneBestArm"){
            c(rep(0,input$numberOfArmsB-1), ifelse(input$effectSizeB != 0, input$effectSizeB, input$effectSizeCustomB))
          } else if(input$effectTrendB == "linear"){
            create.linear.means(num_arms = input$numberOfArmsB, best_arm_mu = as.numeric(ifelse(input$effectSizeB != 0, input$effectSizeB, input$effectSizeCustomB)))
          } else if(input$effectTrendB == "custom"){
            c(0,as.numeric(unlist(strsplit(input$customMusB,","))))
          }), 
        
        n_patients_start = input$nStart,
        n_patients_stop = input$nEnd,
        n_patients_step = input$nStep,
        post_sample_size = input$postSampleSizeB,
        simulation_runs_H0 = input$simsH0B,
        simulation_runs_H1 = input$simsH1B, 
        threads = input$threadsB,
        num_checked_sizes = input$checkSampleSizeB, 
        seed = ifelse(input$randomSeedB == FALSE, input$seedB, randomSeedB),
        delta = input$thresholdFutility_deltaB,
        epsilon = input$thresholdFutility_epsilonB
      )
      
      runTimeStopB = difftime(Sys.time(), runTimeStartB, units = "auto")
      
      output$powerPlotB <- renderPlot({plot.sim2(sobj)})
    })
    
    #out[["blub"]] = sobj$results.smooth
    tout = sobj$results.smooth
    tout$power = round(tout$power,digits=4)
    
    logStringB = paste0(
      
      if(max(sobj$results.smooth[,"power"]) > 0.8) {paste0("Minimal total sample size for 80% power = ", sobj$results.smooth[sobj$results.smooth[,"power"] > 0.8,][1,1], "\n", sep = "")},
      if(max(sobj$results.smooth[,"power"]) > 0.9) {paste0("Minimal total sample size for 90% power = ", sobj$results.smooth[sobj$results.smooth[,"power"] > 0.9,][1,1], "\n", sep = "")},
      "Critical value alpha = ", input$alphaB, "\n",
      "Total sample size = ", paste( seq(from = input$nStart, to = input$nEnd, by = input$nStep), collapse=", "), "\n",
      
      if(input$effectTrendB == "oneBestArm"){
        paste0("Effect size (one best arm) = ", ifelse(input$effectSizeB == 0, input$effectSizeCustomB, input$effectSizeB), "\n")
      } else if(input$effectTrendB == "linear"){
        paste0("Effect size (linear trend) = ", ifelse(input$effectSizeB == 0, input$effectSizeCustomB, input$effectSizeB), "\n")
      },
      "H1 treatment means = ", paste(round(as.numeric(
        if(input$effectTrendB == "oneBestArm"){
          c(rep(0,input$numberOfArmsB-1), ifelse(input$effectSizeB != 0, input$effectSizeB, input$effectSizeCustomB))
        } else if(input$effectTrendB == "linear"){
          create.linear.means(num_arms = input$numberOfArmsB, best_arm_mu = as.numeric(ifelse(input$effectSizeB != 0, input$effectSizeB, input$effectSizeCustomB)))
        } else if(input$effectTrendB == "custom"){
          c(0,as.numeric(unlist(strsplit(input$customMusB,","))))
        }),2), collapse = ", "), "\n",
      
      "Post sample size = ", input$postSampleSizeB, "\n",
      "Simulation runs H0 = ", input$simsH0B, "\n",
      "Simulation runs H1 = ", input$simsH1B, "\n",
      # "Confidence level for simulated alpha and power = ", input$simsCIB, "\n",
      # "Simulation runs other = ", input$simsFreqB, "\n",
      "Number of n1 and n2 combinations per N = ", input$checkSampleSizeB, "\n",
      "Futility-Delta  = ", input$thresholdFutility_deltaB, "\n",
      "Futility-Epsilon  = ", input$thresholdFutility_epsilonB, "\n",
      "Number of threads = ", input$threadsB, "\n",
      "Seed = ", ifelse(input$advancedSettingsB == 1 & input$randomSeedB == FALSE, input$seedB, randomSeedB), "\n",
      "Run time = ", round(runTimeStopB,1), " ", attr(runTimeStopB, "units"),  
      sep = "")
    
    output$logB <- renderText({logStringB})
    
    tout
  }) # processbar done
  
  output$bestN <- renderTable({runProgram1()[["table1"]]}, rownames = T, digits = 3)
  
#  cistring = paste0("\"Power (",round(input$simsCI,3),"% CI)\"")
#  eval(parse(text=cistring))
  
  output$powerTable <- renderDataTable({table2 = data.frame(
    "thisWillBeChanged" = paste0(round(runProgram1()[["sobj"]]$results.all$power,4), " (", round(runProgram1()[["sobj"]]$results.all$power.low,4), ", ", round(runProgram1()[["sobj"]]$results.all$power.high,4), ")", sep = ""),
    "N per arm (Stage 1)" = runProgram1()[["sobj"]]$results.all$n.per.arm.I, 
    "Additional N per arm (Stage 2)" = runProgram1()[["sobj"]]$results.all$n.per.arm.II, 
    "Total N" = runProgram1()[["sobj"]]$results.all$n.total, 
    "Threshold" = round(runProgram1()[["sobj"]]$results.all$threshold,5),
    "Pr(Early stop: futility)" = round(runProgram1()[["sobj"]]$results.all$pr.futility,5),
    "Pr(success, best arm wins)" = round(runProgram1()[["sobj"]]$results.all$pr.best.tx.success,5),
    check.names = F)
    colnames(table2)[1] = paste("Power (", input$simsCI * 100, "% CI)", sep = "")
    table2
  })
  
  output$overall <- renderDataTable({table3 = runProgram2()
  colnames(table3) = c("Total N", "Power")
  table3}, options = list(pageLength = 10) )
  
  # hide("goButton")
}



ui = fluidPage(
  
  tags$head(tags$style(HTML(
    ".shiny-notification {
                            height: 100px;
                            width: 600px;
                            position:fixed;
                            font-size: 500%;
                            text-align: center;
                            line-height: 75px;
                            top: calc(50% - 50px);;
                            left: calc(50% - 300px);;}"
  ))),
  
  h2("Allocation & Power Optimizer"),
  h3("Adaptive Bayes Two-Stage Drop the Losers Design"),
  h5("Authors: Alex Karanevich, Richard Meier, Stefan Graw"),
  h5("Department of Biostatistics, University of Kansas School of Medicine"),
  h5("View the ",a("Instructions",target="_blank",href="instructions.pdf")),
  
  tabsetPanel( ########################################################################################################################################### FIRST TAB
    type = "tabs", 
    tabPanel(
      "Total sample size", 
      
      sidebarLayout( 
        sidebarPanel( # Box 2
          
          # show sample size input
          # Box 2 | sample size
          numericInput(inputId = "alphaB", label = HTML("Critical value &alpha; (one-sided)"), value = 0.05, min = 0, max = 1, step = 0.01),
          # uiOutput("totalSampleSize"), # dynamic min calculated
          
          numericInput(inputId = "nStart", label = "Start point of total sample size", value = 20, min = 5, step = 1),
          numericInput(inputId = "nEnd", label = "End point of total sample size", value = 100, min = 10, step = 1),
          numericInput(inputId = "nStep", label = "Step size", value = 10, min = 1, step = 1),
          
          radioButtons(inputId = "effectTrendB", label = "Effect size as", choices = list("Best arm" = "oneBestArm", "Linear trend" = "linear", "Custom" = "custom")),
          conditionalPanel(
            condition = "input.effectTrendB == 'oneBestArm' || input.effectTrendB == 'linear'",
            numericInput(inputId = "numberOfArmsB", label = "Number of arms (including control)", value = 3, step = 1, min = 2),
            radioButtons(inputId = "effectSizeB", label = "Standardized effect size", choices = list("Small (d=0.2)" = 0.2, "Medium (d=0.5)" = 0.5, "Large (d=0.8)" = 0.8, "Custom" = 0)),
            # CARE effect size has 2 variables now XXXXX
            conditionalPanel(
              condition = "input.effectSizeB == 0",
              numericInput(inputId = "effectSizeCustomB", label = "Standardized effect size", value = 0.5, min = 0, max = 2, step = 0.1) # https://en.wikipedia.org/wiki/Effect_size#Cohen's_d
            )
          ),
          conditionalPanel(
            condition = "input.effectTrendB == 'custom'",
            textInput(inputId = "customMusB", label = HTML("H<sub>1</sub> treatment means </br> (comma separated; control excluded)") , value = "0.2, 0.5")
          ),
          
          numericInput(inputId = "threadsB", label = "Number of cores/threads", value = round(detectCores(all.tests = FALSE, logical = TRUE)/2), min = 1, max = detectCores(all.tests = FALSE, logical = TRUE), step = 1) # XXXXX set default to 1
          , width = 3      
        ),
        
        # power plot
        mainPanel(
          # htmlOutput(outputId = "alphaCIB"),
          fluidRow(
            plotOutput(outputId = "powerPlotB"),
            verbatimTextOutput(outputId = "logB")
            
            # column(5,dataTableOutput(outputId = "overall")),
            # column(7,plotOutput(outputId = "powerPlotB"))
          )
          
        )
      ),
      
      # Box 3
      sidebarLayout(
        sidebarPanel(
          checkboxInput(inputId = "advancedSettingsB", label = "Advanced settings"),
          conditionalPanel( # Box 2 | power
            condition = "input.advancedSettingsB == 1",
            HTML("<b>Futility Stage 1: Stop trial early if<br /><div align='center'>P(&mu;<sub>Best</sub> - &mu;<sub>Ctrl</sub> &ge; <font color='red'>&delta;</font>) &le; <font color='red'>&epsilon;</font></div></b>"), 
            wellPanel(
              numericInput(inputId = "thresholdFutility_deltaB", label = HTML("<font color='red'>&delta;</font> = ?"), value = 0, step = 0.1),
              numericInput(inputId = "thresholdFutility_epsilonB", label = HTML("<font color='red'>&epsilon;</font> = ?"), value = 0, min = 0, max = 1, step = 0.1)
            ),
            
            numericInput(inputId = "checkSampleSizeB", label = HTML("Number of n<sub>1</sub> and n<sub>2</sub> combinations per N"), value = 15, min = 10, step = 1),
            
            numericInput(inputId = "postSampleSizeB", label = "Post sample size", value = 300, min = 100, step = 100), 
            numericInput(inputId = "simsH0B", label = HTML("Simulation runs H<sub>0</sub>"), value = 500, min = 100, step = 100), 
            numericInput(inputId = "simsH1B", label = HTML("Simulation runs H<sub>1</sub>"), value = 700, min = 100, step = 100),
            checkboxInput(inputId = "randomSeedB", label = "Random seed", value = T),
            conditionalPanel(
              condition = "input.randomSeedB == 0",
              numericInput(inputId = "seedB", label = "Seed", value = 1234, min = 1, step = 1)
            )
          ), 
          actionButton("goButtonB", "Go!"),
          width = 3),
        
        mainPanel(
           # "Log:",
           # verbatimTextOutput(outputId = "logB"),
          
          dataTableOutput(outputId = "overall"),
          #  
          #  checkboxInput(inputId = "showTableB", label = "Show table"),
          #  conditionalPanel(
          #    condition = "input.showTableB == 1",
          #    # "tabular:",
          #    dataTableOutput(outputId = "powerTableB")
          #  )
          ""
        )
      )
      
    ), 
    tabPanel( ########################################################################################################################################### SECOND TAB
      "Stage allocation optimization", 
      
      sidebarLayout( 
        sidebarPanel( # Box 2
          
          # show sample size input
          # Box 2 | sample size
          numericInput(inputId = "alpha", label = HTML("Critical value &alpha; (one-sided)"), value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(inputId = "totalSampleSize", label = "Total sample size", value = 100, step = 1, min = 10), 
          # uiOutput("totalSampleSize"), # dynamic min calculated
          
          radioButtons(inputId = "effectTrend", label = "Effect size as", choices = list("Best arm" = "oneBestArm", "Linear trend" = "linear", "Custom" = "custom")),
          conditionalPanel(
            condition = "input.effectTrend == 'oneBestArm' || input.effectTrend == 'linear'",
            numericInput(inputId = "numberOfArms", label = "Number of arms (including control)", value = 3, step = 1, min = 2),
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
          , width = 3      
        ),
        
        # power plot
        mainPanel(
          htmlOutput(outputId = "alphaCI"),
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
              numericInput(inputId = "thresholdFutility_delta", label = HTML("<font color='red'>&delta;</font> = ?"), value = 0, step = 0.1),
              numericInput(inputId = "thresholdFutility_epsilon", label = HTML("<font color='red'>&epsilon;</font> = ?"), value = 0, min = 0, max = 1, step = 0.1)
            ),
            
            # numericInput(inputId = "superiority", label = HTML("Minimum treatment improvment:<br /> P(&mu;<sub>Best</sub> - &mu;<sub>Ctrl</sub> > <font color='red'>?</font>) > Th"), value = 0, min = 0), # max = 2*effect size
            #uiOutput("superiority"), 
            numericInput(inputId = "checkSampleSize", label = HTML("Number of n<sub>1</sub> and n<sub>2</sub> combinations"), value = 20, min = 10, step = 1),
            numericInput(inputId = "postSampleSize", label = "Post sample size", value = 1000, min = 100, step = 100), 
            numericInput(inputId = "simsH0", label = HTML("Simulation runs H<sub>0</sub>"), value = 1000, min = 100, step = 100), 
            numericInput(inputId = "simsH1", label = HTML("Simulation runs H<sub>1</sub>"), value = 10000, min = 100, step = 1000),
			numericInput(inputId = "simsCI", label = HTML("Confidence level for simulated &alpha; and power"), value = 0.95, min = 0.01, max = 0.999, step = 0.01),
			numericInput(inputId = "simsFreq", label = HTML("Simulation runs Dunnett/Bonferroni"), value = 10000, min = 100, step = 1000),
            checkboxInput(inputId = "randomSeed", label = "Random seed", value = T),
            conditionalPanel(
              condition = "input.randomSeed == 0",
              numericInput(inputId = "seed", label = "Seed", value = 1234, min = 1, step = 1)
            )
          ), 
          actionButton("goButton", "Go!"),
          width = 3),
        
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
      
    ) #################################################################################################################################### END OF TABS
    
  )
  
)

shinyApp(ui = ui, server = server)
