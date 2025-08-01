#' @method diagnose summary.hdlm
#' @rdname diagnose
diagnose.summary.hdlm <- function(x, ...) {
  
  if (is.null(x$mcmc.samples)){
    stop("MCMC samples are missing. Make sure to set `mcmc = T` when running `summary()` function.")
  }
  mcmc.samples <- x$mcmc.samples
  
  # Function to generate selectInput elements based on column names
  generateSelectInputs <- function(mod_names) {
    mod_data <- x$data[, mod_names]
    
    inputs <- lapply(mod_names, function(mod) {
      mod_col <- mod_data[, mod]
      
      if (x$modIsNum[mod]) {
        if (all(mod_col %in% c(0, 1))) {
          radioButtons(
            paste0("ind_", mod),
            paste0(mod),
            choices = c(0, 1),
            selected = 0,
            inline = TRUE
          )
        } else if (all(mod_col > 0)) {
          # mod_scale
          sliderInput(
            paste0("ind_", mod),
            paste0(mod),
            min = floor(min(mod_col)),
            max = ceiling(max(mod_col)),
            value = round(apply(
              mod_col %>% as.data.frame(), 2, median
            ), 2)
          )
        } else {
          # mod_cont
          sliderInput(
            paste0("ind_", mod),
            paste0(mod),
            min = floor(min(mod_col)),
            max = ceiling(max(mod_col)),
            value = 0
          )
        }
      } else {
        # For non-continuous columns, add buttons
        selectInput(paste0("ind_", mod), paste0(mod), choices = unique(mod_col))
      }
    })
    
    do.call(tagList, inputs)
  }
  
  

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
            helpText(
              "Select the level of modifiers and click the button below to collect MCMC samples."
            ),
            tags$h4("Modifier input: "),
            generateSelectInputs(x$mod.names),
            actionButton("generate", "Collect MCMC samples"),
          ),
          mainPanel(fluidRow(
            column(6, plotOutput("dlm.trace")), 
            column(6, plotOutput("ce.trace.plot"))
          ), fluidRow(
            column(6, plotOutput("dlm.density")), 
            column(6, plotOutput("ce.density"))
          ))
        ),
        tabPanel(
          "Lag effects",
          sidebarPanel(
            helpText(
              "This tab displays the MCMC samples for each lag based on the modifier level specified in the 'Model' tab."
            ),
            selectInput(
            "lag.choice", "Select lag:", choices = 1:x$n.lag
          )),
          mainPanel(
            fluidRow(
              column(6, plotOutput("lag.trace")),
              column(6, plotOutput("lag.density")) 
            ),
            fluidRow(column(6, plotOutput("lag.acf")),
                     column(6, wellPanel(
                       tableOutput("lag.table"),
                       uiOutput("rule.of.thumb1")
                     )))
          )
        ),
        tabPanel(
          "Tree",
          fluidRow(column(7, plotOutput(
            "modtree.size.plot"
          )), column(5, plotOutput("modaccept.plot"))),
          fluidRow(column(7, plotOutput(
            "dlmtree.size.plot"
          )), column(5, plotOutput("dlmaccept.plot")))
        ),
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
        )
      )
    ),
    
    server = function(input, output, session) {
      
      # Compute heterogeneous effect estimation
      het.data <- eventReactive(input$generate, {
        TreeStructs <- mcmc.samples$dlm.mcmc
        
        # Modifier data frame
        mod <- data.frame(setNames(lapply(x$mod.names, function(col)
          input[[paste0("ind_", col)]]), x$mod.names))
        n <- nrow(mod)
        
        # Build 'draws' list for each exposure
        draws <- lapply(1:max(TreeStructs$Iter), function(i)
          matrix(0.0, n, x$n.lag))
        
        withProgress(message = 'Calculating...', value = 0, {
          for (i in 1:nrow(TreeStructs)) {
            if ((i %% floor(nrow(TreeStructs) / 100)) == 0) {
              incProgress(1 / 100, detail = " Collecting MCMC samples")
            }
            
            # Extract mcmc iteration count and the rule
            Iter <- TreeStructs$Iter[i]
            Rule <- TreeStructs$Rule[i]
            if (Rule == "") {
              idx <- 1:n
            } else {
              idx <- which(eval(parse(text = Rule)))
            }
            
            t   <- TreeStructs$tmin[i]:TreeStructs$tmax[i]
            est <- TreeStructs$est[i]
            draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
          }
        })
        
        # posterior mcmc calculation
        draws <- array(do.call(c, draws), c(n, x$n.lag, max(TreeStructs$Iter)))
        
        t(draws[1, , ])
      })
      
      
      # *** Tab 1: Model diagnostics ***
      output$dlm.trace <- renderPlot({
        req(het.data())
        
        # Extract data
        dlm.fit <- het.data()
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
        dlm.fit <- het.data()
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
          theme_bw(base_size = 16) +
          coord_flip()
      })
      
      
      # Trace plot using ggplot2
      output$ce.trace.plot <- renderPlot({
        dlm.fit <- het.data()
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
        dlm.fit <- het.data()
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
        dlm.fit <- het.data()
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
        dlm.fit <- het.data()
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
        dlm.fit <- het.data()
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
        dlm.fit <- het.data()
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
      tree1.df <- as.data.frame(mcmc.samples$dlmtree.size)
      tree2.df <- as.data.frame(mcmc.samples$modtree.size)
      
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
      output$dlmtree.size.plot <- renderPlot({
        ggplot(tree1.df_long, aes(x = tree, y = 1, fill = size)) + 
          geom_bar(stat = "identity", show.legend = TRUE) +
          scale_fill_gradient(low = "green",
                              high = "#e74c3c",
                              name = "Tree Size") + 
          labs(title = "Number of terminal nodes (leaves) of DLM tree across MCMC iterations", 
               x = "Tree", y = "MCMC iterations") +
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
      output$modtree.size.plot <- renderPlot({
        ggplot(tree2.df_long, aes(x = tree, y = 1, fill = size)) + 
          geom_bar(stat = "identity", show.legend = TRUE) +
          scale_fill_gradient(low = "green",
                              high = "#e74c3c",
                              name = "Tree Size") + 
          labs(title = "Number of terminal nodes (leaves) of modifier tree across MCMC iterations", 
               x = "Tree", y = "MCMC iterations",
               caption = "More red indicates a larger tree, which may require caution, especially for trees with more than 8 leaves.") +
          theme_bw(base_size = 16) +
          coord_flip() +
          scale_y_continuous(expand = c(0, 0))
      })
      
      
      # MH acceptance rate
      dlmtree.accept <- mcmc.samples$dlm.accept
      modtree.accept <- mcmc.samples$mod.accept
      
      output$dlmaccept.plot <- renderPlot({
        dlmtree.accept <- dlmtree.accept %>% as.data.frame() %>%
          mutate(
            step = recode(
              step,
              "0" = "grow",
              "1" = "prune",
              "2" = "change"
            ),
            decision = ifelse(success < 2, "reject", "accept")
          ) %>%
          mutate(step = factor(
            step,
            levels = c("grow", "prune", "change"),
            ordered = TRUE
          ))
        
        # Create a plot
        ggplot(dlmtree.accept, aes(x = factor(step), fill = decision)) +
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
            caption = "y-axis: MCMC iteraction x Number of DLM trees"
          ) +
          theme_bw(base_size = 16)
      })
      
      output$modaccept.plot <- renderPlot({
        modtree.accept <- modtree.accept %>% as.data.frame() %>%
          mutate(
            step = recode(
              step,
              "0" = "grow",
              "1" = "prune",
              "2" = "change",
              "3" = "swap"
            ),
            decision = ifelse(success < 2, "reject", "accept")
          ) %>%
          mutate(step = factor(
            step,
            levels = c("grow", "prune", "change", "swap"),
            ordered = TRUE
          ))
        
        # Create a plot
        ggplot(modtree.accept, aes(x = factor(step), fill = decision)) +
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
            caption = "y-axis: MCMC iteraction x Number of modifier trees"
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
    }
  )
}
