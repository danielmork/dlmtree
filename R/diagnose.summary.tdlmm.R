#' @method diagnose summary.tdlmm
#' @rdname diagnose
diagnose.summary.tdlmm <- function(x, ...) {
  `%notin%` <- Negate(`%in%`)

  if (is.null(x$mcmc.samples)){
    stop("MCMC samples are missing. Make sure to set `mcmc = T` when running `summary()` function.")
  }
  mcmc.samples <- x$mcmc.samples
  
  # Shiny app components
  shinyApp(
    ui = fluidPage(
      theme = shinytheme("flatly"),
      tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"),
      navbarPage(
        "Diagnostics",
        tabPanel(
          "Model",
          sidebarPanel(
            selectInput("exp.choice", "Select exposure:", choices = x$exp.names, )
          ),
          mainPanel(fluidRow(
            column(6, plotOutput("dlm.trace")), column(6, plotOutput("ce.trace.plot"))
          ), fluidRow(
            column(6, plotOutput("dlm.density")), column(6, plotOutput("ce.density"))
          ))
        ),
        tabPanel(
          "Lag effects",
          sidebarPanel(
            selectInput("exp.lag.choice", "Select exposure:", choices = x$exp.names, ),
            selectInput("lag.choice", "Select lag:", choices = 1:nrow(x$DLM[[1]]$mcmc))
          ),
          mainPanel(fluidRow(
            column(6, plotOutput("lag.trace")),
            column(6, plotOutput("lag.density"))  
          ), fluidRow(
            column(6, plotOutput("lag.acf")),
            column(6, wellPanel(
              tableOutput("lag.table"),
              uiOutput("rule.of.thumb1")
            ))
          ))
        ),
        tabPanel("Tree", fluidRow(
          fluidRow(column(7, plotOutput("tree1.size.plot")), 
                   column(5, plotOutput("accept1.plot"))), 
          fluidRow(column(7, plotOutput("tree2.size.plot")), 
                   column(5, plotOutput("accept2.plot")))
        )),
        tabPanel(
          "Hyperparameters",
          sidebarPanel(
            selectInput(
              "hyper.choice",
              "Select parameter:",
              choices = colnames(mcmc.samples$hyper)
            )
          ),
          mainPanel(fluidRow(
            column(6, plotOutput("hyper.trace")), 
            column(6, plotOutput("hyper.density"))  
          ), fluidRow(
            column(6, plotOutput("hyper.acf")),
            column(6, wellPanel(
              tableOutput("hyper.table"), 
              uiOutput("rule.of.thumb2") 
            ))
          ))
        ),
        tabPanel("Exposure selection", fluidRow(column(
          10, plotOutput("cate.trace")
        )))
        
      )
    ),
    
    server = function(input, output, session) {
      # *** Tab 1: Model diagnostics ***
      output$dlm.trace <- renderPlot({
        exp.choice <- input$exp.choice
        
        # Extract data
        dlm.fit <- t(x$DLM[[exp.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        iterations <- 1:nrow(dlm.fit)
        fit.mean <- apply(dlm.fit, 2, mean)
        fit.ci <- apply(dlm.fit, 2, function(x)
          quantile(x, probs = c(0.025, 0.975)))
        
        # Prepare data frame for ggplot
        fit.df <- data.frame(
          iteration = rep(iterations, ncol(dlm.fit)),
          time = rep(1:ncol(dlm.fit), each = nrow(dlm.fit)),
          fit = as.vector(dlm.fit)
        )
        
        fit.summary <- data.frame(
          time = 1:ncol(dlm.fit),
          mean = fit.mean,
          lower = fit.ci[1, ],
          upper = fit.ci[2, ]
        )
        
        ggplot() +
          geom_line(
            data = fit.df,
            aes(x = time, y = fit, group = iteration),
            color = "grey",
            alpha = 0.5
          ) +
          geom_line(
            data = fit.summary,
            aes(x = time, y = mean),
            color = "#e74c3c",
            linewidth = 1
          ) +
          geom_line(
            data = fit.summary,
            aes(x = time, y = lower),
            color = "#0073b7",
            linetype = "dashed"
          ) +
          geom_line(
            data = fit.summary,
            aes(x = time, y = upper),
            color = "#0073b7",
            linetype = "dashed"
          ) +
          labs(
            title = "Posterior samples of distributed lag effects",
            x = "Lag",
            y = "Effect",
            caption = "Grey line: 1 MCMC iteration"
          ) +
          theme_bw(base_size = 16)
      })
      
      
      output$dlm.density <- renderPlot({
        exp.choice <- input$exp.choice
        
        # Extract data
        dlm.fit <- t(x$DLM[[exp.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        iterations <- 1:nrow(dlm.fit)
        fit.mean <- apply(dlm.fit, 2, mean)
        fit.ci <- apply(dlm.fit, 2, function(x)
          quantile(x, probs = c(0.025, 0.975)))
        
        # Prepare data frame for ggplot
        fit.df <- data.frame(
          iteration = rep(iterations, ncol(dlm.fit)),
          time = rep(1:ncol(dlm.fit), each = nrow(dlm.fit)),
          fit = as.vector(dlm.fit)
        )
        
        fit.summary <- data.frame(
          time = 1:ncol(dlm.fit),
          mean = fit.mean,
          lower = fit.ci[1, ],
          upper = fit.ci[2, ]
        )
        
        # Reshape to long format
        fit_long <- dlm.fit %>%
          as.data.frame() %>%
          pivot_longer(cols = everything(),
                       names_to = "lag",
                       values_to = "effect") %>%
          mutate(lag = factor(lag, levels = 1:ncol(dlm.fit)))
        
        # Plot with ggridges
        ggplot(fit_long, aes(x = effect, y = lag)) +
          geom_density_ridges_gradient(
            scale = 1.5,
            rel_min_height = 0.01,
            fill = "#00a65a",
            color = "black"
          ) +
          labs(title = "Density plot of effects across lags", x = "Effect", y = "Lag") +
          theme_bw(base_size = 16) +  # Rotate x-axis labels for better readability
          coord_flip()
      })
      
      
      # Cumulative effects
      output$ce.trace.plot <- renderPlot({
        exp.choice <- input$exp.choice
        
        # Extract data
        dlm.fit <- t(x$DLM[[exp.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        iterations <- 1:nrow(dlm.fit)
        fit.mean <- apply(dlm.fit, 2, mean)
        fit.ci <- apply(dlm.fit, 2, function(x)
          quantile(x, probs = c(0.025, 0.975)))
        
        # Prepare data frame for ggplot
        fit.df <- data.frame(
          iteration = rep(iterations, ncol(dlm.fit)),
          time = rep(1:ncol(dlm.fit), each = nrow(dlm.fit)),
          fit = as.vector(dlm.fit)
        )
        
        fit.summary <- data.frame(
          time = 1:ncol(dlm.fit),
          mean = fit.mean,
          lower = fit.ci[1, ],
          upper = fit.ci[2, ]
        )
        
        # Use coda for cumulative effects
        ce <- rowSums(dlm.fit)
        mcmc.obj <- coda::mcmc(ce)
        cum.mean <- mean(ce)
        cum.ci <- quantile(ce, probs = c(0.025, 0.975))
        
        trace.df <- data.frame(iteration = 1:length(ce), value = ce)
        
        ggplot(trace.df, aes(x = iteration, y = value)) +
          geom_line(color = "grey") +
          geom_hline(
            yintercept = cum.mean,
            color = "#e74c3c",
            linewidth = 1
          ) +
          geom_hline(
            yintercept = cum.ci[1],
            color = "#0073b7",
            linetype = "dashed"
          ) +
          geom_hline(
            yintercept = cum.ci[2],
            color = "#0073b7",
            linetype = "dashed"
          ) +
          labs(title = "Traceplot of cumulative effects with posterior mean and 95% CrI", 
               x = "MCMC iteration", y = "Cumulative effect") +
          theme_bw(base_size = 16)
      })
      
      
      output$ce.density <- renderPlot({
        exp.choice <- input$exp.choice
        
        # Extract data 
        dlm.fit <- t(x$DLM[[exp.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        iterations <- 1:nrow(dlm.fit)
        fit.mean <- apply(dlm.fit, 2, mean)
        fit.ci <- apply(dlm.fit, 2, function(x)
          quantile(x, probs = c(0.025, 0.975)))
        
        # Prepare data frame for ggplot
        fit.df <- data.frame(
          iteration = rep(iterations, ncol(dlm.fit)),
          time = rep(1:ncol(dlm.fit), each = nrow(dlm.fit)),
          fit = as.vector(dlm.fit)
        )
        
        fit.summary <- data.frame(
          time = 1:ncol(dlm.fit),
          mean = fit.mean,
          lower = fit.ci[1, ],
          upper = fit.ci[2, ]
        )
        
        # Use coda for cumulative effects
        ce <- rowSums(dlm.fit)
        mcmc.obj <- coda::mcmc(ce)
        cum.mean <- mean(ce)
        cum.ci <- quantile(ce, probs = c(0.025, 0.975))
        
        trace.df <- data.frame(iteration = 1:length(ce), value = ce)
        
        # Calculate 95% CI
        ce.ci <- quantile(trace.df$value, probs = c(0.025, 0.975))
        
        ggplot(trace.df, aes(x = value)) +
          geom_density(fill = "#00a65a", alpha = 0.5) +
          geom_vline(
            xintercept = ce.ci,
            linetype = "dashed",
            color = "#0073b7",
            linewidth = 1
          ) +
          labs(
            title = "Density plot for cumulative effect",
            x = "Cumulative effects",
            caption = paste(
              "95% Credible interval: [",
              round(ce.ci[1], 2),
              ", ",
              round(ce.ci[2], 2),
              "]"
            )
          ) +
          theme_bw(base_size = 16)
      })
      
      
      # *** Tab 2: Convergence Diagnostics ***
      output$lag.trace <- renderPlot({
        dlm.fit <- t(x$DLM[[input$exp.lag.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        selected.lag <- input$lag.choice
        lag.mcmc <- dlm.fit[, selected.lag]
        
        # Calculate posterior mean and 95% CI
        posterior_mean <- mean(lag.mcmc)
        ci <- quantile(lag.mcmc, probs = c(0.025, 0.975))
        
        trace.df <- data.frame(iteration = 1:nrow(dlm.fit), value = lag.mcmc)
        
        ggplot(trace.df, aes(x = iteration, y = value)) +
          geom_line(color = "grey") +  # Trace lines
          geom_hline(
            yintercept = posterior_mean,
            color = "#e74c3c",
            linetype = "solid",
            linewidth = 1
          ) +  
          geom_hline(
            yintercept = ci[1],
            color = "#0073b7",
            linetype = "dashed",
            linewidth = 1
          ) +  
          geom_hline(
            yintercept = ci[2],
            color = "#0073b7",
            linetype = "dashed",
            linewidth = 1
          ) +  
          labs(
            title = paste("Traceplot for lag", selected.lag),
            x = "MCMC Iteration",
            y = paste("Effect sampled at lag", selected.lag),
            caption = paste("Posterior mean: ", round(posterior_mean, 2))
          ) +
          theme_bw(base_size = 16)
      })
      
      output$lag.density <- renderPlot({
        dlm.fit <- t(x$DLM[[input$exp.lag.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        selected.lag <- input$lag.choice
        lag.mcmc <- dlm.fit[, selected.lag]
        
        # Calculate 95% CI
        lag.ci <- quantile(lag.mcmc, probs = c(0.025, 0.975))
        
        density.df <- data.frame(value = lag.mcmc)
        
        ggplot(density.df, aes(x = value)) +
          geom_density(fill = "#00a65a", alpha = 0.5) +
          geom_vline(
            xintercept = lag.ci,
            linetype = "dashed",
            color = "#0073b7",
            linewidth = 1
          ) +
          labs(
            title = paste("Density plot for lag", selected.lag),
            x = paste("Effect sampled at lag", selected.lag),
            caption = paste("95% CI: [", round(lag.ci[1], 2), ", ", round(lag.ci[2], 2), "]")
          ) +
          theme_bw(base_size = 16)
      })
      
      output$lag.acf <- renderPlot({
        dlm.fit <- t(x$DLM[[input$exp.lag.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        selected.lag <- input$lag.choice
        lag.mcmc <- dlm.fit[, selected.lag]
        
        acf.data <- acf(lag.mcmc, plot = FALSE)
        acf.df <- data.frame(lag = acf.data$lag, acf = acf.data$acf)
        ggplot(acf.df, aes(x = lag, y = acf)) +
          geom_bar(stat = "identity",
                   fill = "#00a65a",
                   alpha = 0.5) +
          labs(
            title = paste("Autocorrelation plot for lag", selected.lag),
            x = "Lag",
            y = "Autocorrelation"
          ) +
          theme_bw(base_size = 16)
      })
      
      output$lag.table <- renderTable({
        dlm.fit <- t(x$DLM[[input$exp.lag.choice]]$mcmc)
        colnames(dlm.fit) <- 1:ncol(dlm.fit)
        selected.lag <- input$lag.choice
        lag.mcmc <- dlm.fit[, selected.lag]
        
        ess <- effectiveSize(lag.mcmc)
        geweke <- geweke.diag(lag.mcmc)$z
        heidelberg <- heidel.diag(lag.mcmc)[, 3]
        acf.data <- acf(lag.mcmc, plot = FALSE)
        acf.value <- acf.data$acf[2]  # lag 1 autocorrelation
        
        data.frame(
          ESS = ess,
          Geweke = geweke,
          Heidelberger.Welch = heidelberg,
          Autocorrelation = acf.value
        )
      })
      
      output$rule.of.thumb1 <- renderUI({
        tagList(
          h4("Rule of Thumb for Convergence Diagnostics:"),
          p(
            HTML(
              "1. <strong>ESS (Effective Sample Size):</strong> ESS > 200 is generally good."
            )
          ),
          p(
            HTML(
              "2. <strong>Geweke Diagnostic:</strong> Values within [-1.96, 1.96] indicate convergence."
            )
          ),
          p(
            HTML(
              "3. <strong>Heidelberger-Welch Test:</strong> p-value > 0.05 suggests convergence."
            )
          ),
          p(
            HTML(
              "4. <strong>Autocorrelation:</strong> Low autocorrelation (close to 0) suggests good mixing."
            )
          )
        )
      })
      
      
      # *** Tab 3: Tree size ***
      tree1.df <- as.data.frame(mcmc.samples$tree1.size)
      tree2.df <- as.data.frame(mcmc.samples$tree2.size)
      
      # Reshape the data for ggplot (long format)
      tree1.df_long <- reshape(
        tree1.df,
        varying = names(tree1.df)[1:x$ctr$n.trees],
        v.names = "size",
        timevar = "tree",
        times = 1:x$ctr$n.trees,
        direction = "long"
      )
      tree1.df_long$tree <- as.factor(tree1.df_long$tree)
      
      # tree size bar plot
      output$tree1.size.plot <- renderPlot({
        ggplot(tree1.df_long, aes(x = tree, y = 1, fill = size)) + 
          geom_bar(stat = "identity", show.legend = TRUE) +  
          scale_fill_gradient(low = "green",
                              high = "#e74c3c",
                              name = "Tree Size") +  
          labs(title = "Number of terminal nodes (leaves) of DLM tree 1 across MCMC iterations", 
               x = "Tree pair", y = "MCMC iterations") +
          theme_bw(base_size = 16) +
          coord_flip() +
          scale_y_continuous(expand = c(0, 0))
      })
      
      # Reshape the data for ggplot (long format)
      tree2.df_long <- reshape(
        tree2.df,
        varying = names(tree2.df)[1:x$ctr$n.trees],
        v.names = "size",
        timevar = "tree",
        times = 1:x$ctr$n.trees,
        direction = "long"
      )
      tree2.df_long$tree <- as.factor(tree2.df_long$tree)
      
      # tree size bar plot
      output$tree2.size.plot <- renderPlot({
        ggplot(tree2.df_long, aes(x = tree, y = 1, fill = size)) +  
          geom_bar(stat = "identity", show.legend = TRUE) +
          scale_fill_gradient(low = "green",
                              high = "#e74c3c",
                              name = "Tree Size") +  
          labs(title = "Number of terminal nodes (leaves) of DLM tree 2 across MCMC iterations", 
               x = "Tree pair", y = "MCMC iterations",
               caption = "More red indicates a larger tree, which may require caution, especially for trees with more than 8 leaves.") +
          theme_bw(base_size = 16) + 
          coord_flip() +
          scale_y_continuous(expand = c(0, 0))
      })
      

      # MH accept ratio
      data.accept <- as.data.frame(mcmc.samples$accept)  
      
      tree1.accept <- subset(data.accept, tree == 1)
      tree2.accept <- subset(data.accept, tree == 2)
      
      output$accept1.plot <- renderPlot({
        tree1.accept <- tree1.accept %>% as.data.frame() %>%
          mutate(
            step = recode(
              step,
              "0" = "grow",
              "1" = "prune",
              "2" = "change",
              "3" = "switch-exposure"
            ),
            decision = ifelse(success < 2, "reject", "accept")
          ) %>%
          mutate(step = factor(
            step,
            levels = c("grow", "prune", "change", "switch-exposure"),
            ordered = TRUE
          ))
        
        # Create a plot
        ggplot(tree1.accept, aes(x = factor(step), fill = decision)) +
          geom_bar(position = "stack", stat = "count") +
          scale_fill_manual(values = c(
            "reject" = "#e74c3c",
            "accept" = "#00a65a"
          )) +
          labs(
            title = "Metropolis-Hastings acceptance",
            x = "Tree steps",
            y = "Decision counts",
            fill = "Decision"
          ) +
          theme_bw(base_size = 16)
      })
      
      output$accept2.plot <- renderPlot({
        tree2.accept <- tree2.accept %>% as.data.frame() %>%
          mutate(
            step = recode(
              step,
              "0" = "grow",
              "1" = "prune",
              "2" = "change",
              "3" = "switch-exposure"
            ),
            decision = ifelse(success < 2, "reject", "accept")
          ) %>%
          mutate(step = factor(
            step,
            levels = c("grow", "prune", "change", "switch-exposure"),
            ordered = TRUE
          ))
        
        # Create a plot
        ggplot(tree2.accept, aes(x = factor(step), fill = decision)) +
          geom_bar(position = "stack", stat = "count") +
          scale_fill_manual(values = c(
            "reject" = "#e74c3c",
            "accept" = "#00a65a"
          )) +
          labs(
            title = "Metropolis-Hastings acceptance",
            x = "Tree steps",
            y = "Decision counts",
            fill = "Decision",
            caption = "y-axis: MCMC iteraction x Number of trees"
          ) +
          theme_bw(base_size = 16)
      })
      

      
      # *** Tab 4: Hyperparameters ***
      hyper.data <- mcmc.samples[["hyper"]]
      
      output$hyper.trace <- renderPlot({
        hyper.selected <- input$hyper.choice
        hyper.mcmc <- hyper.data[, hyper.selected]
        
        # Calculate posterior mean and 95% CI
        posterior_mean <- mean(hyper.mcmc)
        ci <- quantile(hyper.mcmc, probs = c(0.025, 0.975))
        
        trace.df <- data.frame(iteration = 1:length(hyper.mcmc),
                               value = hyper.mcmc)
        
        ggplot(trace.df, aes(x = iteration, y = value)) +
          geom_line(color = "grey") +  
          geom_hline(
            yintercept = posterior_mean,
            color = "#e74c3c",
            linetype = "solid",
            linewidth = 1
          ) + 
          geom_hline(
            yintercept = ci[1],
            color = "#0073b7",
            linetype = "dashed",
            linewidth = 1
          ) + 
          geom_hline(
            yintercept = ci[2],
            color = "#0073b7",
            linetype = "dashed",
            linewidth = 1
          ) +
          labs(
            title = paste("Traceplot for", hyper.selected),
            x = "MCMC Iteration",
            y = "Posterior sample",
            caption = paste("Posterior mean: ", round(posterior_mean, 2))
          ) +
          theme_bw(base_size = 16)
      })
      
      output$hyper.density <- renderPlot({
        hyper.selected <- input$hyper.choice
        hyper.mcmc <- hyper.data[, hyper.selected]
        
        # Calculate 95% CI
        hyper.ci <- quantile(hyper.mcmc, probs = c(0.025, 0.975))
        density.df <- data.frame(value = hyper.mcmc)
        
        ggplot(density.df, aes(x = value)) +
          geom_density(fill = "#00a65a", alpha = 0.5) +
          geom_vline(
            xintercept = hyper.ci,
            linetype = "dashed",
            color = "#0073b7",
            linewidth = 1
          ) +
          labs(
            title = paste("Density Plot for", hyper.selected),
            x = "Posterior sample",
            caption = paste(
              "95% CI: [",
              round(hyper.ci[1], 2),
              ", ",
              round(hyper.ci[2], 2),
              "]"
            )
          ) +
          theme_bw(base_size = 16)
      })
      
      output$hyper.acf <- renderPlot({
        hyper.selected <- input$hyper.choice
        hyper.mcmc <- hyper.data[, hyper.selected]
        
        acf.data <- acf(hyper.mcmc, plot = FALSE)
        acf.df <- data.frame(lag = acf.data$lag, acf = acf.data$acf)
        ggplot(acf.df, aes(x = lag, y = acf)) +
          geom_bar(stat = "identity",
                   fill = "#00a65a",
                   alpha = 0.5) +
          labs(
            title = paste("Autocorrelation Plot for", hyper.selected),
            x = "Lag",
            y = "Autocorrelation"
          ) +
          theme_bw(base_size = 16)
      })
      
      output$hyper.table <- renderTable({
        hyper.selected <- input$hyper.choice
        hyper.mcmc <- hyper.data[, hyper.selected]
        
        ess <- effectiveSize(hyper.mcmc)
        geweke <- geweke.diag(hyper.mcmc)$z
        heidelberg <- heidel.diag(hyper.mcmc)[, 3]
        acf.data <- acf(hyper.mcmc, plot = FALSE)
        acf.value <- acf.data$acf[2]  # lag 1 autocorrelation
        
        data.frame(
          ESS = ess,
          Geweke = geweke,
          Heidelberger.Welch = heidelberg,
          Autocorrelation = acf.value
        )
      })
      
      output$rule.of.thumb2 <- renderUI({
        tagList(
          h4("Rule of Thumb for Convergence Diagnostics:"),
          p(
            HTML(
              "1. <strong>ESS (Effective Sample Size):</strong> ESS > 200 is generally good."
            )
          ),
          p(
            HTML(
              "2. <strong>Geweke Diagnostic:</strong> Values within [-1.96, 1.96] indicate convergence."
            )
          ),
          p(
            HTML(
              "3. <strong>Heidelberger-Welch Test:</strong> p-value > 0.05 suggests convergence."
            )
          ),
          p(
            HTML(
              "4. <strong>Autocorrelation:</strong> Low autocorrelation (close to 0) suggests good mixing."
            )
          )
        )
      })
      
      
      # *** Tab 5: Exposure selection ***
      output$cate.trace <- renderPlot({
        exp.count.df <- mcmc.samples$exp.count
        
        # Reshape data to long format and add row numbers
        df_long <- exp.count.df %>% as.data.frame() %>%
          pivot_longer(cols = everything(),
                       names_to = "Exposure",
                       values_to = "Count") %>%
          mutate(Row.num = rep(1:nrow(exp.count.df), each = ncol(exp.count.df)))
        
        # Create the plot
        ggplot(df_long, aes(x = Row.num, y = Count, fill = Exposure)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_brewer(palette = "Set2") +
          labs(
            title = "Exposures assigned to DLM trees ",
            x = "MCMC iteration",
            y = "Number of trees",
            fill = "Exposure"
          ) +
          theme_bw(base_size = 16) +
          scale_x_continuous(expand = c(0, 0))
      })
    }
  )
}

