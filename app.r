library(gridExtra)
library(tidyverse)
library(deSolve)
library(ggplot2)
library(plotly)
library(shinythemes)
library(shinyjs)
#library(covdata)

## Configuration
theme_set(theme_minimal(base_size = 12))

## CONSTANTS

# Population size 
N <- 1394000

# Rate at which person stays in the infectious compartment (disease specific and tracing specific)
gamma <- 0.14

# Initial number of infected people
n_init <- 1

x_max = 365
max_time <- 365
# # Grid where to evaluate
# max_time <- 365 # 150
# times <- seq(0, max_time, by = 0.1)

# R0 for the beta and gamma values
# R0 <- beta*N/gamma

# calculate beta
# Infectious contact rate - beta = R0/N*gamma and when R0  ~2.25 then  2.25/N*gamma
# beta <- 4.5e-07 

# Function to compute the derivative of the ODE system
# -----------------------------------------------------------
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system (beta, gamma)

diffsdp<-function(sdp)
{
  # Compute the square of integer `n`
  dif <- as.numeric(as.Date(sdp) - as.Date("2020-03-03"))
  return(dif)
}


sir <- function(t, y, parms, 
                social_dist_period, 
                reduction) {
  
  beta0 <- parms[1]
  gamma <- parms[2]

  #print(social_dist_period)
  
  if (t <= social_dist_period[1]) {
    beta_t = beta0;
    #cat('1. This is loop number',t)
  }
  else if (t >=social_dist_period[1] && t <= social_dist_period[2]) {
    beta_t =  beta0 * reduction[1];
    #cat('2. This is loop number',t)
  }
  else if (t >=social_dist_period[2] && t <= social_dist_period[3]) {
    beta_t =  beta0 * reduction[2];
    #cat('3. This is loop number',t)
  }
  else if (t >=social_dist_period[3] && t <= social_dist_period[4]) {
    beta_t =  beta0 * reduction[3];
    #cat('3. This is loop number',t)
  }
  else {
    beta_t =  beta0 * reduction[4];
    #cat('4. This is loop number',t)
  }
  
  # Reduce contact rate 
  #beta_t <- if_else(t <= social_dist_period[1], 
  #                  beta0,
  #                  if_else(t <= social_dist_period[2], 
  #                          beta0 * reduction[1],
  #                          beta0 * reduction[2]
  #                  )
  #)
  
  S <- y[1]
  I <- y[2]
  R <- y[3]
  return(list(c(S = -beta_t * S * I, 
                I =  beta_t * S * I - gamma * I,
                R =  gamma * I)))
}

## we assume that some globals exist...
solve_ode <- function(sdp, red, typ, beta, max_time) {
  
  # Grid where to evaluate
  times <- reactive({ seq(0, max_time, by = 0.1) })
  
  ode_solution <- lsoda(y = c(N - n_init, n_init, 0), 
                        times = times(), 
                        func  = sir, 
                        parms = c(beta, gamma), 
                        social_dist_period = sdp,
                        reduction = red) %>%
    as.data.frame() %>%
    setNames(c("t", "S", "I", "R")) %>%
    mutate(beta = beta, 
           gama = gamma,
           R0 = N * beta / gamma, 
           s  = S, 
           i  = I, 
           type = typ)
  
  daily <- ode_solution %>%
    filter(t %in% seq(0, max_time, by = 1)) %>%
    mutate(C = if_else(row_number() == 1, 0, lag(S, k = 1) - S), 
           c = C)
  
  daily
}


# add results with intervention and plot
run <- function(sdp, red, r0, max_time) {
  #print(sdp)
  #print("sdp")
  #print(sdp2)
  #print("sdp2")
  print(red)
  shinyjs::logjs(red)
  #print("red")
  
  beta <- r0 / N * gamma
  
  ode_solution_daily <- solve_ode(
    sdp = c(0, max_time,max_time,max_time),  # social_dist_period
    red = c(1, 1, 1, 1),         # reduction
    typ = "If there was no Lockdown", 
    beta = beta, 
    max_time
  )
  
  # solve with interventions
  ode_solution2_daily <- solve_ode(
    sdp = sdp,
    red = red,
    typ = "with Lockdown", 
    beta = beta, 
    max_time
  )
  
  # solve with interventions
  #ode_solution3_daily <- solve_ode(
  #  sdp = sdp2,
  #  red = c(0.66, 0.8),
  #  typ = "with Lockdown2", 
  #  beta = beta, 
  #  max_time
  #)  
  
  # Combine the two solutions into one dataset
  ode_df <- rbind(ode_solution_daily, ode_solution2_daily)
  
}

plot_lockdown <- function(ode_df, sdp) { 
  #print("Plot Result")
  #print(sdp)
  #print(ode_df)
  #write.csv(x=ode_df, file="ode_df")
  max_time = 365
  # The final size in the two cases:
  final_sizes <- ode_df %>%
    group_by(type) %>%
    filter(row_number() == n()) %>%
    mutate("Proportion Infected" = scales::percent(1 - s, accuracy = 1)) %>%
    select("Proportion Infected", interventions = type) %>% 
    arrange(desc(interventions))
  
  # Plot
  y_axis_fixed = FALSE
  if (y_axis_fixed) {
    y_max <- 0.09
  } else {
    y_max <- max(ode_df$c, na.rm = TRUE) * 1.05
  }
  y_arrow <- y_max * 0.975
  y_text  <- y_arrow + y_max * 0.01 
  col_sdp <- "deeppink4"
  
  x_labs <- sort(c(0, 150, 200, 300, 365, sdp))
  start = as.Date("2020-03-03")
  ode_df$date <- as.Date(ode_df$t+start)
  
  pp <-ggplotly(
    ggplot(ode_df, 
           aes(x = t, 
               y = 0, 
               xend = t, 
               yend = c, 
               color = type,
               name = date)) + 
      geom_segment(alpha = 0.7) + 
      geom_line(aes(x = t, y = c)) + 
      labs(
        x = "Days from first known case", 
        y = "Forecasts for Infected Daily", 
        subtitle = "Daily new cases", 
        caption  = "SIQ: Separate, Isolate and Quarantine period") +
      scale_x_continuous(labels = x_labs, 
                         breaks = x_labs) +
      scale_y_continuous(limits = c(0, y_max), labels = scales::comma) + 
      #scale_color_brewer(name = "Key", 
      #                   type = "qual", 
      #                   palette = 6, 
      #                   guide = guide_legend(reverse = TRUE)) +
      # sdp 
      #geom_hline(yintercept = 10900, color = "deeppink4") +
      geom_vline(xintercept = sdp, lty = 2, color = "gray30") +
      geom_text(aes(x = 10 + sdp[1] + (sdp[2] - sdp[1])/2, 
                    y = y_text, 
                    label = "With No Lockdown"),
                vjust = 0,
                color = col_sdp) +   
      geom_text(aes(x = sdp[2] + 100, 
                    y = y_text, 
                    label = "With Lockdowns"),
                vjust = 0,
                color = "springgreen3") +   
      geom_segment(aes(
        x = sdp[1],
        y = y_arrow, 
        xend = sdp[2] * 0.99, # shorten
        yend = y_arrow
      ),
      size = 0.3, 
      color = col_sdp,
      arrow = arrow(length = unit(2, "mm"))) +
      geom_segment(aes(
        x = sdp[2]*1.01, # shorten
        y = y_arrow, 
        xend = max_time,
        yend = y_arrow
      ),
      size = 0.3, 
      color = col_sdp,
      arrow = arrow(length = unit(2, "mm")))  +
      # Add final size as table
      annotation_custom(tableGrob(final_sizes, 
                                  rows = NULL,
                                  theme = ttheme_minimal(
                                    core    = list(fg_params = list(hjust = 0, x = 0.1)),
                                    rowhead = list(fg_params = list(hjust = 0, x = 0))
                                  )),
                        xmin = max_time * 0.4,
                        xmax = max_time,
                        ymin = y_max * 0.8,
                        ymax = y_max * 0.8
      )
  )
  ggplotly(pp, tooltip="c")
}

plot_predict <- function(ode_df, sdp) {    
  print(sdp)
  max_time = 365
  
  x_labs <- sort(c(0, 200, 300, 365, sdp))
  start = as.Date("2020-03-03")
  ends = as.Date("2020-03-03")
  ode_df$date <- as.Date(ode_df$t+start)
  #filtered <- subset(ode_df, type == "with SIQ")
  result <- filter(ode_df, type == "with Lockdown")
  result$Dead <-result$C * 0.005
  total_deaths = sum(result$Dead)
  total_cases = sum(result$C)
  total_detected_cases = sum(result$C) * 0.2
  y_max <- max(result$c, na.rm = TRUE) * 1.05
  y_arrow <- y_max * 0.975
  y_text  <- y_arrow + y_max * 0.01 
  col_sdp <- "deeppink4"
  
  result = head(result, 365)
  names(result)[names(result) == "C"] <- "daily_infections"
  result$CumulativeCases <-cumsum(result$daily_infections)
  #print(result)
  
  pp1 <-ggplotly(
    ggplot(result, aes(date, daily_infections)) +
      geom_bar(stat="identity", na.rm = TRUE, fill = "darkturquoise") + 
      labs(
        x = "Date", 
        y = "Infected Daily") + #ggtitle("Predicted Daily Infections with Lockdown Interventions") + 
        #title = "Predicted Daily Infections with Lockdown Interventions") + 
      geom_vline(xintercept=as.numeric(ode_df$date[sdp]), colour="gray30",linetype="dotted") +
      scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%B") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
      theme(axis.text=element_text(size=8),
            axis.title=element_text(size=8)) +
      annotate(geom="text",x=start + 45,
               y=y_max, label="Full Lockdown",fontface="bold", color = "deeppink4") + 
      annotate(geom="text",x=ends + sdp[2] + 4,
               y=y_max,label="Ph 1",fontface="bold",color = "deeppink4") + 
      annotate(geom="text",x=ends + sdp[3] + 4,
               y=y_max,label="Ph 2",fontface="bold",color = "deeppink4") +
      annotate(geom="text",x=ends + sdp[4] + 4,
               y=y_max,label="Ph 3",fontface="bold",color = "deeppink4") +
      annotate(geom="text",x=ends + 300,
               y=y_max, label=sprintf("Total Cases = %s",toString(ceiling(total_cases))),fontface="bold",color = "red") +
    annotate(geom="text",x=ends + 300,
             y=y_max * .9, label=sprintf("Total Detected Cases = %s",toString(ceiling(total_detected_cases))),fontface="bold",color = "red") +
    annotate(geom="text",x=ends + 300,
             y=y_max * .8, label=sprintf("Total Dead = %s",toString(ceiling(total_deaths))),fontface="bold",color = "red")
      + scale_y_continuous(labels = scales::comma)
    
  )
  ggplotly(pp1, tooltip="daily_infections")
}

make_summary_table <- function(ode_df, sdp) {
  
  start = as.Date("2020-03-03")
  ends = as.Date("2020-03-03") + sdp[1]
  ode_df$date <- as.Date(ode_df$t+start)
  prediction <- filter(ode_df, type == "with Lockdown")
  
  prediction$Dead <-prediction$C * 0.005
  prediction$Detected_Cases <-prediction$C * 0.1
  
  summary <- prediction %>%mutate(Year_Month = format(date, "%Y-%m")) %>%
    group_by(Year_Month) %>%
    summarise(total_predicted_cases = ceiling(sum(C)), total_detected_cases = ceiling(sum(Detected_Cases)), total_predicted_deaths = ceiling(sum(Dead)))
  
  summary <- summary[order(summary$Year_Month),]
  
  return(summary)
  
}

plot_real <- function(ode_df, sdp) {    
  guy  <- read.csv(file = 'Trinidad.csv')
  guy$date <- as.Date(guy$date)
  pp3 <-ggplotly(
    ggplot(guy, aes(x=date)) +
      geom_line(aes(y = cu_cases, color = "Cumulative Cases")) + 
      geom_line(aes(y = cu_deaths, color = "Cumulative Deaths")) +
      geom_point(aes(y = cases, color = "Cases", alpha=1/3), shape=1) +
      geom_point(aes(y = deaths, color = "Deaths",alpha=1/3), shape = 4) + 
      scale_colour_manual("", 
                          breaks = c("Cumulative Cases", "Cumulative Deaths", "Cases","Deaths"),
                          values = c("turquoise3", "red", "turquoise3","red")) #+
    #scale_y_continuous(labels = scales::comma) + geom_vline(xintercept=as.numeric(ode_df$date[sdp]), colour="gray30") +
    #annotate(geom="text",x=start + 50,
    #         y=top,label="Lockdown starts",fontface="bold", color = "deeppink4") + 
    #annotate(geom="text",x=ends + sdp[1] + 30,
    #         y=top,label="Phase 1",fontface="bold",color = "deeppink4") + 
    #annotate(geom="text",x=ends + sdp[2] + 30,
    #         y=top,label="Phase 2",fontface="bold",color = "deeppink4") 
  )
  
  ggplotly(pp3)
}

plot_cum <- function(ode_df, sdp) {    
  
  max_time = 365
  
  x_labs <- sort(c(0, 200, 300, 365, sdp))
  start = as.Date("2020-03-03")
  ends = as.Date("2020-03-03") + sdp[1]
  ode_df$date <- as.Date(ode_df$t+start)
  #filtered <- subset(ode_df, type == "with SIQ")
  result <- filter(ode_df, type == "with Lockdown")
  y_max <- max(result$c, na.rm = TRUE) * 1.05
  y_arrow <- y_max * 0.975
  y_text  <- y_arrow + y_max * 0.01 
  
  col_sdp <- "deeppink4"
  
  result = head(result, 365)
  names(result)[names(result) == "C"] <- "daily_infections"
  result$CumulativeCases <-cumsum(result$daily_infections)
  result$Dead <-result$daily_infections * 0.005
  result$CumulativeSus <-result$S - result$Dead
  result$CumulativeDead <-cumsum(result$Dead)
  result$CumulativeRec <- result$S - result$CumulativeCases
  top = head(result$CumulativeSus, n=1) * 1.1
  #top = result["CumulativeCases"].iloc[-1]
  #print(result)
  
  pp3 <-ggplotly(
    ggplot(result, aes(x=date)) +
      geom_line(aes(y = CumulativeCases, color = "Cases")) + 
      geom_line(aes(y = CumulativeDead, color = "Dead")) +
      geom_line(aes(y = CumulativeSus, color = "Susceptible")) +
      geom_line(aes(y = R, color = "Recovered")) + 
      scale_colour_manual("", 
                          breaks = c("Cases", "Dead", "Susceptible","Recovered"),
                          values = c("steelblue", "darkred", "goldenrod", "green")) +
      #scale_color_discrete(name = "Descriptions", labels = c("CumulativeCases", "CumulativeDead", "CumulativeSus", "R")) +
      scale_y_continuous(labels = scales::comma) + geom_vline(xintercept=as.numeric(ode_df$date[sdp]), colour="gray30") +
      annotate(geom="text",x=start + 40,
               y=top,label="Lockdown",fontface="bold", color = "deeppink4") + 
      annotate(geom="text",x=ends + sdp[2] - 26,
               y=top,label="Ph 1",fontface="bold",color = "deeppink4") + 
      annotate(geom="text",x=ends + sdp[3] - 26,
               y=top,label="Ph 2",fontface="bold",color = "deeppink4") +
      annotate(geom="text",x=ends + sdp[4] - 26,
               y=top,label="Ph 3",fontface="bold",color = "deeppink4") 
  )
  
  ggplotly(pp3, tooltip="CumulativeCases")
}


ui <- fluidPage(theme = shinytheme("flatly"),
  titlePanel("What Happens When T&T Re-Opens? An Interactive Covid-19 SIR Policy Simulator"),
  h4("Charts update in real-time, please allow ~5 seconds for charts to load"),
  hr(),
  sidebarLayout(
    sidebarPanel(
      h3("Phase 1 of Lockdown Removal"),
      br(),
      dateInput("end",
                "Set start date:",
                value = "2020-06-08",
                format = "dd/mm/yyyy"),
      br(),
      h4("Select Policy Choices"),
      checkboxInput("foo", label = "Food Establistments and related stores re-opened", value = TRUE),
      checkboxInput("non", label = "Manufacturing & Public Section re-started", value = TRUE),
      checkboxInput("bus", label = "Medium to low risk businesses and malls re-opened", value = TRUE),
      checkboxInput("res", label = "Medium to high risk businesses re-openes (resturants, bars, cinemas, salons, offices & gyms)", value = FALSE),
      checkboxInput("gat", label = "Mass Gatherings of 25+ allowed", value = FALSE),
      checkboxInput("sch", label = "Schools re-opened", value = FALSE),
      sliderInput("sda",
                  "Control Social Distance/Mask Wearing/Hand Washing Adherance (100% is everyone follows)",
                  min   = 1,
                  max   = 100,
                  value = 60, 
                  step  = 1), 
      hr(),
      h3("Phase 2 of Lockdown Removal"),
      br(),
      dateInput("end2",
                "Set start date:",
                value = "2020-06-22",
                format = "dd/mm/yyyy"),
      br(),
      h4("Select Policy Choices"),
      checkboxInput("foo2", label = "Food Establistments and related stores re-opened", value = TRUE),
      checkboxInput("non2", label = "Manufacturing & Public Section re-started", value = TRUE),
      checkboxInput("bus2", label = "Medium to low risk businesses and malls re-opened", value = TRUE),
      checkboxInput("res2", label = "Medium to high risk businesses re-openes (resturants, bars, cinemas, offices, salons & gyms)", value = TRUE),
      checkboxInput("gat2", label = "Mass Gatherings of 25+ allowed", value = TRUE),
      checkboxInput("sch2", label = "Schools re-opened", value = TRUE),
      sliderInput("sda2",
                  "Control Social Distance/Mask Wearing/Hand Washing Adherance (100% is everyone follows)",
                  min   = 1,
                  max   = 100,
                  value = 50, 
                  step  = 1), 
      hr(),
      h3("Phase 3 of Lockdown Removal"),
      br(),
      dateInput("end3",
                "Set start date:",
                value = "2020-08-17",
                format = "dd/mm/yyyy"),
      br(),
      h4("Select Policy Choices"),
      checkboxInput("foo3", label = "Food Establistments and related stores re-opened", value = TRUE),
      checkboxInput("non3", label = "Manufacturing & Public Section re-started", value = TRUE),
      checkboxInput("bus3", label = "Medium to low risk businesses and malls re-opened", value = TRUE),
      checkboxInput("res3", label = "Medium to high risk businesses re-openes (resturants, bars, cinemas, offices, salons & gyms)", value = FALSE),
      checkboxInput("gat3", label = "Mass Gatherings of 25+ allowed", value = FALSE),
      checkboxInput("sch3", label = "Schools re-opened", value = FALSE),
      sliderInput("sda3",
                  "Control Social Distance/Mask Wearing/Hand Washing Adherance (100% is everyone follows)",
                  min   = 1,
                  max   = 100,
                  value = 70, 
                  step  = 1), 
      #radioButtons("sda", label = h4("Control Social Distance Adherance (Mask Wearing etc)"),
      #             choices = list("0 to 30% People wear masks in public" = 7,
      #                            "30% to 70% People wear masks in public" = 13,
      #                            "70% to 100% People wear masks in public" = 16),
      #             selected = 13),
      #radioButtons("sdo", label = h4("Social Distancing Obeyed"),
      #             choices = list("0 to 30% People object 6ft spacing" = 7,
      #                            "30% to 70% People object 6ft spacing" = 13,
      #                            "70% to 100% People object 6ft spacing" = 15),
      #             selected = 13),
      #radioButtons("hwh", label = h4("Hand Washing and Basic Hygenine"),
      #             choices = list("0 to 30% People follow" = 7,
      #                            "30% to 70% People follow" = 13,
      #                           "70% to 100% People follow" = 15),
      #             selected = 13),
      br(),
      hr(),         
      tags$div("Developed by Rajev Ratan"), 
      tags$div("Version 2.1, 16-May-2020"), 
      tags$a(href="https://caribbeandatafam.com/","Caribbean Data Fam"), 
      
      br(),
      width = 3,
    ),
    
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Forecasted Cases", plotlyOutput("chart1",  height = "auto", width = "auto"), 
                           DT::dataTableOutput('table', width = "55%", height = "auto")),
                  tabPanel("Forecasted Cumulative Data", plotlyOutput("chart3",  height = "auto",width = "auto")),
                  tabPanel("Actual Data", plotlyOutput("chart2",  height = "auto",width = "100%")),
                  tabPanel("What if there were no Lockdowns?", plotlyOutput("chart4",  height = "auto",width = "auto"))
                  #tabPanel("Table", tableOutput("table"))
      ),
      
      
      
      #br(),
      #br(),
      #hr(),
      div(HTML("Parameters Used:")), 
      div(HTML("R<sub>0</sub> of 2.25")), 
      div(HTML("Gamma = 0.2, beta <- 4.5e-07 ")), 
      div(HTML("beta <- 4.5e-07 ")),
      div(HTML("Population = 1,394,000")), 
      div(HTML("Initial Infect People = 2")),
      div(HTML("Day 0 is 3rd March 202")),
      br(),
      width = 9
    )
    
  )
)


server <- function(input, output) {
  
  res <- reactive({
    
    run(sdp = c(c(30, diffsdp(input$end)), diffsdp(input$end2),  diffsdp(input$end3)),
        #sdp2 = c(c(diffsdp(input$end)),c(diffsdp(input$end2))),
        #red = c(0.4, (100-(input$red_two)/1.5)/100), 
        #test =  input$mask,
        red = c(0.3, (((100 - ((-13*as.integer(input$sch)) + (-13*as.integer(input$res)) + (-13*as.integer(input$gat)) +
                                 (-10*as.integer(input$bus))+ (-10*as.integer(input$foo)) + 
                                 (-10*as.integer(input$non))  ))/100)- 0.7) * (1-(as.integer(input$sda)/300)), 
                (((100 - ((-13*as.integer(input$sch2)) + (-13*as.integer(input$res2)) + (-13*as.integer(input$gat2)) +
                            (-10*as.integer(input$bus2))+ (-10*as.integer(input$foo2)) + 
                            (-10*as.integer(input$non2))  ))/100)- 0.7) * (1-(as.integer(input$sda2)/300)),
                (((100 - ((-13*as.integer(input$sch3)) + (-13*as.integer(input$res3)) + (-13*as.integer(input$gat3)) +
                            (-10*as.integer(input$bus3))+ (-10*as.integer(input$foo3)) + 
                            (-10*as.integer(input$non3))  ))/100)- 0.7) * (1-(as.integer(input$sda3)/300))),
        
        #(((100 - ((-5*as.integer(input$sch2)) + (-5*as.integer(input$res2)) + (-5*as.integer(input$gat2)) +
        #           (-5*as.integer(input$bus2))+ (-5*as.integer(input$foo2))+ (-5*as.integer(input$non2)) + (-1*(10-(as.integer(input$sda2)/10) )) ))/100)- 0.4)), 
        #red2 = c(0(((100 - ((-5*as.integer(input$sch)) + (-5*as.integer(input$res)) + (-4*as.integer(input$gat)) +
        #                      (-3*as.integer(input$bus))+ (-3*as.integer(input$foo))+ (-3*as.integer(input$non)) + (-1*(7-(as.integer(input$sda)/10) )) ))/100)- 0.4), 
        #         (((100 - ((-5*as.integer(input$sch2)) + (-5*as.integer(input$res2)) + (-4*as.integer(input$gat2)) +
        #                     (-3*as.integer(input$bus2))+ (-3*as.integer(input$foo2))+ (-3*as.integer(input$non2)) + (-1*(7-(as.integer(input$sda2)/10) )) ))/100)- 0.4)), 
        r0  = 2.25,
        max_time = 365)
    #start  = input$start,
    # 
    #)
  })
  
  output$table <- DT::renderDataTable(make_summary_table(res(), c(c(30,diffsdp(input$end),diffsdp(input$end2),diffsdp(input$end3))) ), options=list(scrollX=T))
  
  output$chart1 <- renderPlotly({
    plot_predict(res(), c(c(30,diffsdp(input$end),diffsdp(input$end2)), diffsdp(input$end3)) )
  })
  
  output$chart2 <- renderPlotly({
    plot_real(res(), c(c(30,diffsdp(input$end),diffsdp(input$end2)), diffsdp(input$end3)) )
  })
  
  
  output$chart3 <- renderPlotly({
    plot_cum(res(), c(c(30,diffsdp(input$end),diffsdp(input$end2)), diffsdp(input$end3)) )
  })
  
  output$chart4 <- renderPlotly({
    plot_lockdown(res(), c(c(30,diffsdp(input$end),diffsdp(input$end2)), diffsdp(input$end3)) )
  })  
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)