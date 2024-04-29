#' shiny.hdlm
#'
#' @description A function to execute a shinyApp to provide comprehensive analysis with HDLM. The shinyApp includes PIP, split points, individualized & subgroup-specific effect.
#'
#' @param fit an object of class 'hdlm'
#'
#' @examples
#' D <- sim.hdlmm(sim = "B", n = 1000)
#' fit <- dlmtree(y ~ ., 
#'                data = D$dat,
#'                exposure.data = D$exposures,
#'                dlm.type = "linear",
#'                family = "gaussian",
#'                het = TRUE)
#' if (interactive()) {
#'  shiny(fit)
#' }
#'
#' @returns A shinyapp interface 
#' @export
shiny.hdlm <- function(fit)
{
  `%notin%` <- Negate(`%in%`)
  
  if(!(class(fit) %in% c("hdlm", "dlmtree"))){
    stop("The class of the model fit must be 'hdlm' or 'dlmtree'")
  }

  # Data preparation
  factor_cols <- sapply(fit$data, is.factor)  # Identify factor columns
  fit$data[factor_cols] <- lapply(fit$data[factor_cols], as.character) 
  fit$data <- fit$data %>% as_tibble()

  
  # Functions for weighted subgroup DLM effect
  plotDLM <- function(estDLM_data, groups = 1) {
    if (groups == 2) {
      estDLM_data$plotData$grp1 <- sapply(strsplit(estDLM_data$plotDat$group, " & "), function(g) g[1])
      estDLM_data$plotData$grp2 <- sapply(strsplit(estDLM_data$plotDat$group, " & "), function(g) g[2])
      
      # Order for nice facet_wrap
      estDLM_data$plotData$grp1 <- factor(estDLM_data$plotData$grp1, levels = unique(estDLM_data$plotData$grp1))
      estDLM_data$plotData$grp2 <- factor(estDLM_data$plotData$grp2, levels = unique(estDLM_data$plotData$grp2))
    }
    
    # Order for nice facet_wrap
    estDLM_data$plotData$group  <- factor(estDLM_data$plotData$group, levels = unique(estDLM_data$plotData$group))
    
    # Plot
    p <- ggplot(estDLM_data$plotData) +
          geom_hline(yintercept = 0, color = "red") +
          geom_ribbon(aes(x = as.numeric(time - 1), ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
          geom_line(aes(x = as.numeric(time - 1), y = est)) +
          facet_wrap(~group)
    if (groups == 2) {
      p <- p + facet_grid(grp1 ~ grp2)
    } else {
      p <- p + facet_wrap(~group)
    }
    p <- p + theme_bw(base_size = 24) +
          scale_x_continuous(expand = c(0, 0)) +
          labs(x = "Lags", y = "Effect")
    p
  }

  createGrpIdx <- function(model, data, mod1) {
    m1    <- which(model$modNames == mod1)
    m1sp  <- model$modSplitValRef[[m1]]
    m1num <- model$modIsNum[m1]
    if (model$modIsNum[m1]) {
      grp1_idx <- lapply(1:(length(m1sp) + 1), function(i) {
        if (i == 1) {
          which(data[[mod1]] < m1sp[i])
        } else if (i <= length(m1sp)) {
          which(data[[mod1]] >= m1sp[i-1] & data[[mod1]] < m1sp[i])
        } else {
          which(data[[mod1]] >= m1sp[i-1])
        }
      })
      names(grp1_idx) <- sapply(1:(length(m1sp) + 1), function(i) {
        if (i == 1) {
          paste0(mod1, " < ", round(m1sp[i], 2))
        } else if (i <= length(m1sp)) {
          paste0(mod1, " in [", round(m1sp[i-1], 2), ", ", round(m1sp[i], 2),")")
        } else {
          paste0(mod1, " >= ", round(m1sp[i-1], 2))
        }
      })
    } else {
      grp1_idx <- lapply(m1sp, function(s) {
        which(data[[mod1]] == s)
      })
      names(grp1_idx) <- sapply(m1sp, function(s) {
        paste0(mod1, " = ", s)
      })
    }
    return(grp1_idx)
  }

  create2GrpIdx <- function(model, data, mod1, mod2) {
    grp1_idx <- createGrpIdx(model, data, mod1)
    grp2_idx <- createGrpIdx(model, data, mod2)
    comb_idx <- list()
    for (i in 1:length(grp1_idx)) {
      for (j in 1:length(grp2_idx)) {
        comb_idx[[paste0(names(grp1_idx)[i], " & ", names(grp2_idx)[j])]] <-
          intersect(grp1_idx[[i]], grp2_idx[[j]])
      }
    }
    return(comb_idx)
  }

  # Function to generate selectInput elements based on column names
  generateSelectInputs <- function(mod_names) {
    mod_data <- fit$data[, mod_names]
    
    inputs <- lapply(mod_names, function(mod) {
      mod_col <- mod_data[, mod]
      
      if (fit$modIsNum[mod]) {
        if (all(mod_col %in% c(0, 1))){
          radioButtons(paste0("ind_", mod), paste0(mod), choices = c(0, 1), selected = 0, inline = TRUE)
        } else if (all(mod_col > 0)){
          # mod_scale
          sliderInput(paste0("ind_", mod), paste0(mod), min = floor(min(mod_col)), max = ceiling(max(mod_col)), value = round(apply(mod_col %>% as.data.frame(), 2, median), 2))
        } else {
          # mod_cont
          sliderInput(paste0("ind_", mod), paste0(mod), min = floor(min(mod_col)), max = ceiling(max(mod_col)), value = 0)
        }
      } else {
        # For non-continuous columns, add buttons
        selectInput(paste0("ind_", mod), paste0(mod), choices = unique(mod_col))
      }
    })
    
    do.call(tagList, inputs)
  }


  ui <- fluidPage(theme = shinytheme("flatly"),
                  tags$style(type = "text/css", ".irs-grid-pol.small {height: 0px;}"),
                  navbarPage("Heterogeneous DLM",
                            tabPanel("Modifier",
                                      sidebarPanel(
                                        selectInput("mod", "Select a numeric modifier", 
                                                    choices = fit$modNames[fit$modIsNum]),
                                        #h3("Help text"),
                                        helpText("Note: The PIP plot is not subject to changes in the selection."),
                                      ), 
                                      mainPanel(
                                        h3("Posterior inclusion probability (PIP) & Split points"),
                                        
                                        fluidRow(
                                          column(6, plotOutput(outputId = "piphist", height = "600px")),
                                          column(6, plotOutput(outputId = "splithist", height = "600px"))
                                        ),
                                      )
                                      
                            ),
                            tabPanel("Individual",
                                      sidebarPanel(
                                        tags$h4("Modifier input: "),
                                        generateSelectInputs(fit$modNames),
                                        actionButton("ind_btn", "Submit")
                                      ), 
                                      mainPanel(
                                        h3("Individualized distributed lag effect"),
                                        
                                        plotOutput(outputId = "individualDLM", height = "600px") # Individualized DLM effect
                                        
                                      )
                                      
                            ),
                            tabPanel("Subgroup",
                                      sidebarPanel(
                                        tags$h4("Individual effect per subgroup input: "),
                                        sliderInput("sub_n", "Sample size per subgroup", min = 0, max = 200, value = 50),
                                        selectInput("subgroup1", "Select a first modifier", choices = fit$modNames),
                                        selectInput("subgroup2", "Select a second modifier", choices = NULL),
                                        helpText("Note: 'No selection' for the second modifier will result in a plot featuring only the first modifier."),
                                        actionButton("sub_btn", "Submit"),
                                        
                                        tags$hr(),
                                        
                                        tags$h4("Weighted subgroup specific effect input: "),
                                        selectInput("subgroup1_w", "Select a first modifier", choices = fit$modNames),
                                        selectInput("subgroup2_w", "Select a second modifier", choices = NULL),
                                        helpText("Caution: Choosing two continuous modifiers will lead to an excessive number of plots"),
                                        actionButton("sub_btn_w", "Submit"),
                                        
                                        
                                      ),
                                      mainPanel(
                                        h3("Subgroup specific effect"),
                                        
                                        plotOutput(outputId = "subgroupDLM", height = "600px"), # Individualized DLM effect
                                        tags$hr(),
                                        plotOutput(outputId = "subgroupWeight", height = "1000px") # Subgroup effect
                                      )
                            )
                  )
  )




  server <- function(input, output, session){
    # *** Panel 1: Split point plot ***
    output$piphist <- renderPlot({
      pip_df <- data.frame("Modifier" = names(pip(fit)), "PIP" = pip(fit))
      
      # Your ggplot code for plot2 here
      ggplot(pip_df, aes(x = reorder(Modifier, PIP), y = PIP)) + 
        geom_col(fill = "aquamarine3", color = "black") + 
        theme_bw(base_size = 20) + 
        labs(x = "Modifier", y = "Posterior inclusion probability") + 
        coord_flip()
    })
    
    output$splithist <- renderPlot({
      # Modifier split points
      spl_df <- splitpoints(fit, var = input$mod, round = 2)
      
      # Bar plot
      split_p <- ggplot(spl_df, aes(x = as.factor(location), y = proportion)) +
        geom_bar(stat = "identity", fill = "skyblue", color = "black") +
        theme_minimal(base_size = 20) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(title = paste0(input$mod), x = "Location", y = "Proportion")
      
      if(nrow(spl_df) != 1){ # Adjust the angle for numeric modifiers
        split_p <- split_p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      }
      
      split_p
    })
    
    
    # *** Panel 2: Individualized DLM ***
    # Build a new data frame (or a vector) called "mod"
    output$individualDLM <- renderPlot({
      input$ind_btn
      
      # Modifier data frame
      mod <- isolate(data.frame(setNames(lapply(fit$modNames, function(col) input[[paste0("ind_", col)]]), fit$modNames)))
      n <- nrow(mod) 
      
      if (input$ind_btn > 0) {
        # Build 'draws' list for each exposure
        draws <- lapply(1:fit$mcmcIter, function(i) matrix(0.0, n, fit$pExp))
        
        withProgress(message = 'Calculating...', value = 0, {
          # Iterate through DLM matrix and add to draws
          for(i in 1:nrow(fit$TreeStructs)){
            if((i %% floor(nrow(fit$TreeStructs)/5)) == 0){
              incProgress(1/5, detail = " individualized effect")
            }
            
            # Extract mcmc iteration count and the rule
            Iter <- fit$TreeStructs$Iter[i]
            Rule <- fit$TreeStructs$Rule[i]
            if (Rule == ""){
              idx <- 1:n
            } else {
              idx <- which(eval(parse(text = Rule)))
            }
            
            t   <- fit$TreeStructs$tmin[i]:fit$TreeStructs$tmax[i]
            est <- fit$TreeStructs$est[i]
            draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
          }
        })
        
        # posterior mcmc calculation
        draws <- array(do.call(c, draws), c(n, fit$pExp, fit$mcmcIter))
        
        # Exposure effect plot
        dlmest        <- sapply(1:(fit$pExp), function(t) {rowMeans(draws[, t, , drop=F])}) # All of n / t = 1 / all mcmc
        dlmest.lower  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = (1 - 0.95)/2)})
        dlmest.upper  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = 1 - (1 - 0.95)/2)})
        
        # Prepare a data.frame for plotting
        dlmest_df <- data.frame("time" = 1:(fit$pExp), 
                                "est" = dlmest, 
                                "est.lower" = dlmest.lower, 
                                "est.upper" = dlmest.upper)
        
        # Generate an exposure effect plot
        ggplot(data = dlmest_df, aes(x = time, y = est)) +
          geom_line(col = "coral", linewidth = 1) +
          geom_ribbon(aes(ymin = est.lower, ymax = est.upper), linetype = "dotted", alpha = 0.1, fill = "coral") +
          theme_bw(base_size = 24) + 
          scale_x_continuous(expand = c(0, 0)) +
          #scale_y_continuous(expand = c(0, 0)) +
          labs(x = "Lag", y = "Effect")
      } else {
        return(NULL) 
      }
    })

    
    # *** Panel 3: subgroup DLM ***
    # Input change so that modifier inputs are not the same
    observeEvent(input$subgroup1, {
      # Update choices for second input based on the selection in the first input
      remaining_options <- setdiff(fit$modNames, input$subgroup1)
      updateSelectInput(session, "subgroup2", choices = c("No selection", remaining_options))
    })
    
    observeEvent(input$subgroup1_w, {
      # Update choices for second input based on the selection in the first input
      remaining_options <- setdiff(fit$modNames, input$subgroup1_w)
      updateSelectInput(session, "subgroup2_w", choices = c("No selection", remaining_options))
    })
    
    # Build a new data frame (or a vector) called "mod"
    output$subgroupDLM <- renderPlot({
      input$sub_btn
      
      # Inputs
      sub_mod1    <- isolate(input$subgroup1)
      sub_mod2    <- isolate(input$subgroup2)
      sample_size <- isolate(input$sub_n)
      
      
      if (input$sub_btn > 0) {
        # Sort out the cases
        if(sub_mod2 == "No selection"){ # Single modifier case
          if (fit$modIsNum[sub_mod1]) { # Single numerical modifier
            # print("num")

            # Find the middle splitpoint
            spl       <- splitpoints(fit, var = sub_mod1, round = 2)
            spl$cumu  <- cumsum(spl$proportion)
            mid_loc   <- spl$location[length(which(spl$cumu < 0.5)) + 1]
            
            mod_num   <- c(paste0(sub_mod1," < ", mid_loc), paste0(mid_loc, " =< ", sub_mod1))
            list_num  <- vector("list", length = length(mod_num))
            
            # DLM calculation
            dlm_list        <- vector("list", length = length(mod_num))
            dlm_lower_list  <- vector("list", length = length(mod_num))
            dlm_upper_list  <- vector("list", length = length(mod_num))
            
            names(list_num) <- names(dlm_list) <- names(dlm_lower_list) <- names(dlm_upper_list) <- mod_num
            
            list_num[[paste0(sub_mod1," < ", mid_loc)]]   <- fit$data[fit$data[, sub_mod1] < mid_loc, ] %>% sample_n(sample_size)
            list_num[[paste0(mid_loc, " =< ", sub_mod1)]] <- fit$data[fit$data[, sub_mod1] >= mid_loc, ] %>% sample_n(sample_size)
            
            withProgress(message = 'Calculating...', value = 0, {
              # DLM for each subgroup
              for(cluster in mod_num){
                incProgress(1/length(mod_num), detail = paste(" for a subgroup: ", cluster))
                
                mod <- list_num[[cluster]]
                n <- nrow(mod)
                
                # Build 'draws' list for each exposure
                draws <- lapply(1:fit$mcmcIter, function(i) matrix(0.0, n, fit$pExp))
                
                # Iterate 
                for(i in 1:nrow(fit$TreeStructs)){
                  
                  # Extract mcmc iteration count and the rule
                  Iter <- fit$TreeStructs$Iter[i]
                  Rule <- fit$TreeStructs$Rule[i]
                  if (Rule == ""){
                    idx <- 1:n
                  } else {
                    idx <- which(eval(parse(text = Rule)))
                  }
                  
                  t   <- fit$TreeStructs$tmin[i]:fit$TreeStructs$tmax[i]
                  est <- fit$TreeStructs$est[i]
                  draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
                }
                
                # posterior mcmc calculation
                draws <- array(do.call(c, draws), c(n, fit$pExp, fit$mcmcIter))
                
                # Exposure effect plot
                dlmest        <- sapply(1:(fit$pExp), function(t) {rowMeans(draws[, t, , drop=F])}) # All of n / t = 1 / all mcmc
                dlmest.lower  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = (1 - 0.95)/2)})
                dlmest.upper  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = 1 - (1 - 0.95)/2)})
                
                dlm_list[[cluster]]       <- dlmest
                dlm_lower_list[[cluster]] <- dlmest.lower
                dlm_upper_list[[cluster]] <- dlmest.upper
              }
            })
            
            # Combine DLM to a longer data frame with the cluster assignment
            dlm_df <- as.data.frame(do.call(rbind, dlm_list))
            colnames(dlm_df) <- 1:fit$pExp
            
            # list_num
            list_num[[paste0(sub_mod1," < ", mid_loc)]]   <- list_num[[paste0(sub_mod1," < ", mid_loc)]] %>% mutate(Var1 := paste0(sub_mod1," < ", mid_loc))
            list_num[[paste0(mid_loc, " =< ", sub_mod1)]] <- list_num[[paste0(mid_loc, " =< ", sub_mod1)]] %>% mutate(Var1 := paste0(mid_loc, " =< ", sub_mod1))
            
            dlm_df      <- cbind(dlm_df, as.data.frame(do.call(rbind, list_num))) # Combine modifier information for colors
            dlm_df$Obs  <- rep(1:sample_size, length(mod_num))
            dlm_df      <- dlm_df %>% pivot_longer(1:fit$pExp, names_to = "Week", values_to = "Effect") %>% mutate(Week = as.numeric(Week))
            
            # Plot
            ggplot(dlm_df, aes(x = Week, y = Effect)) +
              geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
              geom_line(aes(group = Obs), alpha = 0.3, col = "dodgerblue") +
              theme_bw(base_size = 24) +
              theme(plot.title = element_text(hjust = 0.5)) +
              labs(x = "Lag", y = "Effect", title = paste0("Individual DLMs per subgroups (n = ", sample_size, " per subgroup)")) +
              scale_x_continuous(expand = c(0, 0)) +
              facet_wrap(as.formula(paste("~ Var1"))) 
            
            
          } else { # categorical modifier            
            mod_cat  <- pull(unique(fit$data[, sub_mod1]))       
            list_cat <- vector("list", length = length(mod_cat))
            
            # DLM calculation
            dlm_list        <- vector("list", length = length(mod_cat))
            dlm_lower_list  <- vector("list", length = length(mod_cat))
            dlm_upper_list  <- vector("list", length = length(mod_cat))
            
            names(list_cat) <- names(dlm_list) <- names(dlm_lower_list) <- names(dlm_upper_list) <- mod_cat
            
            for(cluster in mod_cat){
              list_cat[[cluster]] <- fit$data[fit$data[, sub_mod1] == mod_cat[which(mod_cat == cluster)], ] %>% sample_n(sample_size)
            }
            
            withProgress(message = 'Calculating...', value = 0, {
              # DLM for each subgroup
              for(cluster in mod_cat){
                incProgress(1/length(mod_cat), detail = paste("effect for individuals in a subgroup: ", cluster))
                
                mod <- list_cat[[cluster]]
                n <- nrow(mod)
                
                # Build 'draws' list for each exposure
                draws <- lapply(1:fit$mcmcIter, function(i) matrix(0.0, n, fit$pExp))
                
                # Iterate 
                for(i in 1:nrow(fit$TreeStructs)){
                  
                  # Extract mcmc iteration count and the rule
                  Iter <- fit$TreeStructs$Iter[i]
                  Rule <- fit$TreeStructs$Rule[i]
                  if (Rule == ""){
                    idx <- 1:n
                  } else {
                    idx <- which(eval(parse(text = Rule)))
                  }
                  
                  t   <- fit$TreeStructs$tmin[i]:fit$TreeStructs$tmax[i]
                  est <- fit$TreeStructs$est[i]
                  draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
                }
                
                # posterior mcmc calculation
                draws <- array(do.call(c, draws), c(n, fit$pExp, fit$mcmcIter))
                
                # Exposure effect plot
                dlmest        <- sapply(1:(fit$pExp), function(t) {rowMeans(draws[, t, , drop=F])}) # All of n / t = 1 / all mcmc
                dlmest.lower  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = (1 - 0.95)/2)})
                dlmest.upper  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = 1 - (1 - 0.95)/2)})
                
                dlm_list[[cluster]]       <- dlmest
                dlm_lower_list[[cluster]] <- dlmest.lower
                dlm_upper_list[[cluster]] <- dlmest.upper
              }
            })
            
            # Combine DLM to a longer data frame with the cluster assignment
            dlm_df <- as.data.frame(do.call(rbind, dlm_list))
            colnames(dlm_df) <- 1:fit$pExp

            dlm_df      <- cbind(dlm_df, as.data.frame(do.call(rbind, list_cat))) # Combine modifier information for colors
            dlm_df$Obs  <- rep(1:sample_size, length(mod_cat))
            dlm_df      <- dlm_df %>% pivot_longer(1:fit$pExp, names_to = "Week", values_to = "Effect") %>% mutate(Week = as.numeric(Week))
            
            # Plot
            ggplot(dlm_df, aes(x = Week, y = Effect)) +
              geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
              geom_line(aes(group = Obs), alpha = 0.3, col = "dodgerblue") +
              theme_bw(base_size = 24) +
              theme(plot.title = element_text(hjust = 0.5)) +
              labs(x = "Lag", y = "Effect", title = paste0("Individual DLMs per subgroups (n = ", sample_size, " per subgroup)")) +
              scale_x_continuous(expand = c(0, 0)) +
              facet_wrap(as.formula(paste("~", sub_mod1))) 
          }
        } else {
          if (fit$modIsNum[sub_mod1]) { # numerical modifier
            if (fit$modIsNum[sub_mod2]){      # num x num              
              # Find the middle splitpoint
              # modifier 1
              spl       <- splitpoints(fit, var = sub_mod1, round = 2)
              spl$cumu  <- cumsum(spl$proportion)
              mid_loc1  <- spl$location[length(which(spl$cumu < 0.5)) + 1]
              
              # modifier 1
              spl <- splitpoints(fit, var = sub_mod2, round = 2)
              spl$cumu <- cumsum(spl$proportion)
              mid_loc2 <- spl$location[length(which(spl$cumu < 0.5)) + 1]
              
              mod1_num <- c(paste0(sub_mod1, " < ", mid_loc1), paste0(mid_loc1, " =< ", sub_mod1))
              mod2_num <- c(paste0(sub_mod2, " < ", mid_loc2), paste0(mid_loc2, " =< ", sub_mod2))
              comb_num <- expand.grid(mod1_num, mod2_num) %>% mutate(comb = paste0(Var1, " & ", Var2))
              names(comb_num) <- c(sub_mod1, sub_mod2, "comb")
              
              # DLM calculation
              list_num        <- vector("list", length = nrow(comb_num))
              dlm_list        <- vector("list", length = nrow(comb_num))
              dlm_lower_list  <- vector("list", length = nrow(comb_num))
              dlm_upper_list  <- vector("list", length = nrow(comb_num))
              
              names(list_num) <- names(dlm_list) <- names(dlm_lower_list) <- names(dlm_upper_list) <- comb_num$comb
              
              list_num[[paste0(sub_mod1, " < ", mid_loc1, " & ", sub_mod2, " < ", mid_loc2)]] <- fit$data[fit$data[, sub_mod1] < mid_loc1 & fit$data[, sub_mod2] < mid_loc2, ] %>% 
                sample_n(sample_size) %>% mutate(Var1 = paste0(mod1_num[1])) %>% mutate(Var2 = mod2_num[1])
              list_num[[paste0(mid_loc1, " =< ", sub_mod1, " & ", sub_mod2, " < ", mid_loc2)]] <- fit$data[fit$data[, sub_mod1] >= mid_loc1 & fit$data[, sub_mod2] < mid_loc2, ] %>% 
                sample_n(sample_size) %>% mutate(Var1 = mod1_num[2]) %>% mutate(Var2 = mod2_num[1])
              list_num[[paste0(sub_mod1, " < ", mid_loc1, " & ", mid_loc2, " =< ", sub_mod2)]] <- fit$data[fit$data[, sub_mod1] < mid_loc1 & fit$data[, sub_mod2] >= mid_loc2, ] %>% 
                sample_n(sample_size) %>% mutate(Var1 = mod1_num[1]) %>% mutate(Var2 = mod2_num[2])
              list_num[[paste0(mid_loc1, " =< ", sub_mod1, " & ", mid_loc2, " =< ", sub_mod2)]] <- fit$data[fit$data[, sub_mod1] >= mid_loc1 & fit$data[, sub_mod2] >= mid_loc2, ] %>% 
                sample_n(sample_size) %>% mutate(Var1 = mod1_num[2]) %>% mutate(Var2 = mod2_num[2])

              withProgress(message = 'Calculating...', value = 0, {
                # DLM for each subgroup
                for(cluster in comb_num$comb){
                  incProgress(1/length(comb_num$comb), detail = paste("effect for individuals in for a subgroup: ", cluster))
                  
                  mod <- list_num[[cluster]]
                  n   <- nrow(mod)
                  
                  # Build 'draws' list for each exposure
                  draws <- lapply(1:fit$mcmcIter, function(i) matrix(0.0, n, fit$pExp))
                  
                  # Iterate 
                  for(i in 1:nrow(fit$TreeStructs)){
                    # Extract mcmc iteration count and the rule
                    Iter <- fit$TreeStructs$Iter[i]
                    Rule <- fit$TreeStructs$Rule[i]
                    if (Rule == ""){
                      idx <- 1:n
                    } else {
                      idx <- which(eval(parse(text = Rule)))
                    }
                    
                    t   <- fit$TreeStructs$tmin[i]:fit$TreeStructs$tmax[i]
                    est <- fit$TreeStructs$est[i]
                    draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
                  }
                  
                  # posterior mcmc calculation
                  draws <- array(do.call(c, draws), c(n, fit$pExp, fit$mcmcIter))
                  
                  # Exposure effect plot
                  dlmest <- sapply(1:(fit$pExp), function(t) {rowMeans(draws[, t, , drop=F])}) # All of n / t = 1 / all mcmc
                  dlmest.lower <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = (1 - 0.95)/2)})
                  dlmest.upper <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = 1 - (1 - 0.95)/2)})
                  
                  dlm_list[[cluster]]       <- dlmest
                  dlm_lower_list[[cluster]] <- dlmest.lower
                  dlm_upper_list[[cluster]] <- dlmest.upper
                }
              })
              
              # Combine DLM to a longer data frame with the cluster assignment
              dlm_df <- as.data.frame(do.call(rbind, dlm_list))
              colnames(dlm_df) <- 1:fit$pExp

              dlm_df      <- cbind(dlm_df, as.data.frame(do.call(rbind, list_num))) # Combine modifier information for colors
              dlm_df$Obs  <- rep(1:sample_size, nrow(comb_num))
              dlm_df      <- dlm_df %>% pivot_longer(1:fit$pExp, names_to = "Week", values_to = "Effect") %>% mutate(Week = as.numeric(Week))
              
              # Plot
              ggplot(dlm_df, aes(x = Week, y = Effect)) +
                geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
                geom_line(aes(group = Obs), alpha = 0.3, col = "dodgerblue") +
                theme_bw(base_size = 24) +
                theme(plot.title = element_text(hjust = 0.5)) +
                labs(x = "Lag", y = "Effect", title = paste0("Individual DLMs per subgroups (n = ", sample_size, " per subgroup)")) +
                scale_x_continuous(expand = c(0, 0)) +
                facet_grid(Var2 ~ Var1)
              
              
            } else { # num x cat
              # Find the middle splitpoint
              # modifier 1
              spl       <- splitpoints(fit, var = sub_mod1, round = 2)
              spl$cumu  <- cumsum(spl$proportion)
              mid_loc1  <- spl$location[length(which(spl$cumu < 0.5)) + 1]
              
              mod1_num    <- c(paste0(sub_mod1, " < ", mid_loc1), paste0(mid_loc1, " =< ", sub_mod1))
              mod2_cat    <- pull(unique(fit$data[, sub_mod2]))
              comb        <- expand.grid(mod1_num, mod2_cat) %>% mutate(comb = paste0(Var1, " & ", Var2))
              names(comb) <- c(sub_mod1, sub_mod2, "comb")
              
              # DLM calculation
              list_num        <- vector("list", length = nrow(comb))
              dlm_list        <- vector("list", length = nrow(comb))
              dlm_lower_list  <- vector("list", length = nrow(comb))
              dlm_upper_list  <- vector("list", length = nrow(comb))
              
              names(list_num) <- names(dlm_list) <- names(dlm_lower_list) <- names(dlm_upper_list) <- comb$comb
              
              for(cat in mod2_cat){
                list_num[[paste0(sub_mod1, " < ", mid_loc1, " & ", cat)]] <- fit$data[fit$data[, sub_mod1] < mid_loc1 & fit$data[, sub_mod2] == cat, ] %>% 
                  sample_n(sample_size) %>% mutate(Var1 = mod1_num[1])
                list_num[[paste0(mid_loc1, " =< ", sub_mod1, " & ", cat)]] <- fit$data[fit$data[, sub_mod1] >= mid_loc1 & fit$data[, sub_mod2]  == cat, ] %>% 
                  sample_n(sample_size) %>% mutate(Var1 = mod1_num[2])
              }
              
              withProgress(message = 'Calculating...', value = 0, {
                # DLM for each subgroup
                for(cluster in comb$comb){
                  incProgress(1/length(comb$comb), detail = paste("effect for individuals in a subgroup: ", cluster))
                  
                  mod <- list_num[[cluster]]
                  n   <- nrow(mod)
                  
                  # Build 'draws' list for each exposure
                  draws <- lapply(1:fit$mcmcIter, function(i) matrix(0.0, n, fit$pExp))
                  
                  # Iterate 
                  for(i in 1:nrow(fit$TreeStructs)){
                    
                    # Extract mcmc iteration count and the rule
                    Iter <- fit$TreeStructs$Iter[i]
                    Rule <- fit$TreeStructs$Rule[i]
                    if (Rule == ""){
                      idx <- 1:n
                    } else {
                      idx <- which(eval(parse(text = Rule)))
                    }
                    
                    t   <- fit$TreeStructs$tmin[i]:fit$TreeStructs$tmax[i]
                    est <- fit$TreeStructs$est[i]
                    draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
                  }
                  
                  # posterior mcmc calculation
                  draws <- array(do.call(c, draws), c(n, fit$pExp, fit$mcmcIter))
                  
                  # Exposure effect plot
                  dlmest        <- sapply(1:(fit$pExp), function(t) {rowMeans(draws[, t, , drop=F])}) # All of n / t = 1 / all mcmc
                  dlmest.lower  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = (1 - 0.95)/2)})
                  dlmest.upper  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = 1 - (1 - 0.95)/2)})
                  
                  dlm_list[[cluster]]       <- dlmest
                  dlm_lower_list[[cluster]] <- dlmest.lower
                  dlm_upper_list[[cluster]] <- dlmest.upper
                }
              })
              
              # Combine DLM to a longer data frame with the cluster assignment
              dlm_df <- as.data.frame(do.call(rbind, dlm_list))
              colnames(dlm_df) <- 1:fit$pExp

              dlm_df      <- cbind(dlm_df, as.data.frame(do.call(rbind, list_num))) # Combine modifier information for colors
              dlm_df$Obs  <- rep(1:sample_size, nrow(comb))
              dlm_df      <- dlm_df %>% pivot_longer(1:fit$pExp, names_to = "Week", values_to = "Effect") %>% mutate(Week = as.numeric(Week))
              
              # Plot
              ggplot(dlm_df, aes(x = Week, y = Effect)) +
                geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
                geom_line(aes(group = Obs), alpha = 0.3, col = "dodgerblue") +
                theme_bw(base_size = 24) +
                theme(plot.title = element_text(hjust = 0.5)) +
                labs(x = "Lag", y = "Effect", title = paste0("Individual DLMs per subgroups (n = ", sample_size, " per subgroup)")) +
                scale_x_continuous(expand = c(0, 0)) +
                facet_grid(as.formula(paste(sub_mod2, "~ Var1")))
              
            }
          } else { # categorical modifier
            if (fit$modIsNum[sub_mod2]){      # cat x num
              # Find the middle splitpoint
              # modifier 1
              mod1_cat <- pull(unique(fit$data[, sub_mod1]))
              
              # modifier 2
              spl       <- splitpoints(fit, var = sub_mod2, round = 2)
              spl$cumu  <- cumsum(spl$proportion)
              mid_loc2  <- spl$location[length(which(spl$cumu < 0.5)) + 1]
              mod2_num  <- c(paste0(sub_mod2, " < ", mid_loc2), paste0(mid_loc2, " =< ", sub_mod2))
              
              comb        <- expand.grid(mod1_cat, mod2_num) %>% mutate(comb = paste0(Var1, " & ", Var2))
              names(comb) <- c(sub_mod1, sub_mod2, "comb")
              
              # DLM calculation
              list_num        <- vector("list", length = nrow(comb))
              dlm_list        <- vector("list", length = nrow(comb))
              dlm_lower_list  <- vector("list", length = nrow(comb))
              dlm_upper_list  <- vector("list", length = nrow(comb))
              
              names(list_num) <- names(dlm_list) <- names(dlm_lower_list) <- names(dlm_upper_list) <- comb$comb
              
              for(cat in mod1_cat){
                list_num[[paste0(cat, " & ", sub_mod2, " < ", mid_loc2)]] <- fit$data[fit$data[, sub_mod2] < mid_loc2 & fit$data[, sub_mod1] == cat, ] %>% 
                  sample_n(sample_size) %>% mutate(Var1 = mod2_num[1])
                list_num[[paste0(cat, " & ", mid_loc2, " =< ", sub_mod2)]] <- fit$data[fit$data[, sub_mod2] >= mid_loc2 & fit$data[, sub_mod1]  == cat, ] %>% 
                  sample_n(sample_size) %>% mutate(Var1 = mod2_num[2])
              }
              
              withProgress(message = 'Calculating...', value = 0, {
                # DLM for each subgroup
                for(cluster in comb$comb){
                  incProgress(1/length(comb$comb), detail = paste("effect for individuals in a subgroup: ", cluster))
                  
                  mod <- list_num[[cluster]]
                  n <- nrow(mod)
                  
                  # Build 'draws' list for each exposure
                  draws <- lapply(1:fit$mcmcIter, function(i) matrix(0.0, n, fit$pExp))
                  
                  # Iterate 
                  for(i in 1:nrow(fit$TreeStructs)){
                    
                    # Extract mcmc iteration count and the rule
                    Iter <- fit$TreeStructs$Iter[i]
                    Rule <- fit$TreeStructs$Rule[i]
                    if (Rule == ""){
                      idx <- 1:n
                    } else {
                      idx <- which(eval(parse(text = Rule)))
                    }
                    
                    t   <- fit$TreeStructs$tmin[i]:fit$TreeStructs$tmax[i]
                    est <- fit$TreeStructs$est[i]
                    draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
                  }
                  
                  # posterior mcmc calculation
                  draws <- array(do.call(c, draws), c(n, fit$pExp, fit$mcmcIter))
                  
                  # Exposure effect plot
                  dlmest        <- sapply(1:(fit$pExp), function(t) {rowMeans(draws[, t, , drop=F])}) # All of n / t = 1 / all mcmc
                  dlmest.lower  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = (1 - 0.95)/2)})
                  dlmest.upper  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = 1 - (1 - 0.95)/2)})
                  
                  dlm_list[[cluster]]       <- dlmest
                  dlm_lower_list[[cluster]] <- dlmest.lower
                  dlm_upper_list[[cluster]] <- dlmest.upper
                }
              })
              
              # Combine DLM to a longer data frame with the cluster assignment
              dlm_df <- as.data.frame(do.call(rbind, dlm_list))
              colnames(dlm_df) <- 1:fit$pExp

              dlm_df      <- cbind(dlm_df, as.data.frame(do.call(rbind, list_num))) # Combine modifier information for colors
              dlm_df$Obs  <- rep(1:sample_size, nrow(comb))
              dlm_df      <- dlm_df %>% pivot_longer(1:fit$pExp, names_to = "Week", values_to = "Effect") %>% mutate(Week = as.numeric(Week))
              
              # Plot
              ggplot(dlm_df, aes(x = Week, y = Effect)) +
                geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
                geom_line(aes(group = Obs), alpha = 0.3, col = "dodgerblue") +
                theme_bw(base_size = 24) +
                theme(plot.title = element_text(hjust = 0.5)) +
                labs(x = "Lag", y = "Effect", title = paste0("Individual DLMs per subgroups (n = ", sample_size, " per subgroup)")) +
                scale_x_continuous(expand = c(0, 0)) +
                facet_grid(as.formula(paste("Var1 ~ ", sub_mod1)))
              
            } else {      # cat x cat
              mod1_cat        <- pull(unique(fit$data[, sub_mod1]))
              mod2_cat        <- pull(unique(fit$data[, sub_mod2])) # Extract as a vector
              comb_cat        <- expand.grid(mod1_cat, mod2_cat) %>% mutate(comb = paste0(Var1, " & ", Var2))
              names(comb_cat) <- c(sub_mod1, sub_mod2, "comb")
              list_cat        <- vector("list", length = nrow(comb_cat))
              
              # DLM calculation
              dlm_list        <- vector("list", length = nrow(comb_cat))
              dlm_lower_list  <- vector("list", length = nrow(comb_cat))
              dlm_upper_list  <- vector("list", length = nrow(comb_cat))

              names(list_cat) <- names(dlm_list) <- names(dlm_lower_list) <- names(dlm_upper_list) <- comb_cat$comb
              
              for(cluster in 1:nrow(comb_cat)){
                cluter_name <- comb_cat$comb[cluster]
                list_cat[[cluter_name]] <- fit$data[pull(fit$data[, sub_mod1]) == comb_cat[cluster, 1] & 
                                                    pull(fit$data[, sub_mod2]) == comb_cat[cluster, 2], ] %>% sample_n(sample_size)
              }
              
              withProgress(message = 'Calculating...', value = 0, {
                # DLM for each subgroup
                for(cluster in comb_cat$comb){
                  incProgress(1/length(comb_cat$comb), detail = paste("effect for individuals in for a subgroup: ", cluster))
                  
                  mod <- list_cat[[cluster]]
                  n <- nrow(mod)
                  
                  # Build 'draws' list for each exposure
                  draws <- lapply(1:fit$mcmcIter, function(i) matrix(0.0, n, fit$pExp))
                  
                  # Iterate 
                  for(i in 1:nrow(fit$TreeStructs)){
                    
                    # Extract mcmc iteration count and the rule
                    Iter <- fit$TreeStructs$Iter[i]
                    Rule <- fit$TreeStructs$Rule[i]
                    if (Rule == ""){
                      idx <- 1:n
                    } else {
                      idx <- which(eval(parse(text = Rule)))
                    }
                    
                    t   <- fit$TreeStructs$tmin[i]:fit$TreeStructs$tmax[i]
                    est <- fit$TreeStructs$est[i]
                    draws[[Iter]][idx, t] <- draws[[Iter]][idx, t] + est
                  }
                  
                  # posterior mcmc calculation
                  draws <- array(do.call(c, draws), c(n, fit$pExp, fit$mcmcIter))
                  
                  # Exposure effect plot
                  dlmest        <- sapply(1:(fit$pExp), function(t) {rowMeans(draws[, t, , drop=F])}) # All of n / t = 1 / all mcmc
                  dlmest.lower  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = (1 - 0.95)/2)})
                  dlmest.upper  <- sapply(1:(fit$pExp), function(t) {apply(draws[, t, , drop=F], 1, quantile, probs = 1 - (1 - 0.95)/2)})
                  
                  dlm_list[[cluster]]       <- dlmest
                  dlm_lower_list[[cluster]] <- dlmest.lower
                  dlm_upper_list[[cluster]] <- dlmest.upper
                }
              })
              
              # Combine DLM to a longer data frame with the cluster assignment
              dlm_df <- as.data.frame(do.call(rbind, dlm_list))
              colnames(dlm_df) <- 1:fit$pExp

              dlm_df      <- cbind(dlm_df, as.data.frame(do.call(rbind, list_cat))) # Combine modifier information for colors
              dlm_df$Obs  <- rep(1:sample_size, nrow(comb_cat))
              dlm_df      <- dlm_df %>% pivot_longer(1:fit$pExp, names_to = "Week", values_to = "Effect") %>% mutate(Week = as.numeric(Week))
              
              # Plot
              ggplot(dlm_df, aes(x = Week, y = Effect)) +
                geom_hline(yintercept = 0, linetype = "dashed", col = "black") +
                geom_line(aes(group = Obs), alpha = 0.3, col = "dodgerblue") +
                theme_bw(base_size = 24) +
                theme(plot.title = element_text(hjust = 0.5)) +
                labs(x = "Lag", y = "Effect", title = paste0("Individual DLMs per subgroups (n = ", sample_size, " per subgroup)")) +
                scale_x_continuous(expand = c(0, 0)) +
                facet_grid(as.formula(paste(sub_mod2, "~", sub_mod1)))
              
            }
          }
        }
      } else {
        return(NULL)
      }
    })
    
    output$subgroupWeight <- renderPlot({
      input$sub_btn_w
      
      # Inputs
      sub_mod1    <- isolate(input$subgroup1_w)
      sub_mod2    <- isolate(input$subgroup2_w)
      sample_size <- isolate(input$sub_n)
      
      if (input$sub_btn_w > 0) {
        withProgress(message = 'Calculating...', value = 0, {
          if (sub_mod2 == "No selection"){
            # 1 modifier
            incProgress(1/2, detail = paste(" weighted subgroup effect: ", sub_mod1))
            grpDLM <- estDLM(fit, fit$data[complete.cases(fit$data),],
                            createGrpIdx(fit, fit$data[complete.cases(fit$data),], sub_mod1), verbose = F)
            incProgress(1/2, detail = "")
            
            plotDLM(grpDLM) + 
              labs(title = "Weighted subgroup average effect") +
              theme(plot.title = element_text(hjust = 0.5))
          } else {
            # 2 modifiers
            incProgress(1/2, detail = paste(" weighted subgroup effect: ", sub_mod1, " & ", sub_mod2))
            grpDLM <- estDLM(fit, fit$data[complete.cases(fit$data),],
                            create2GrpIdx(fit, fit$data[complete.cases(fit$data),], sub_mod1, sub_mod2), verbose = F)
            incProgress(1/2, detail = "")
            
            plotDLM(grpDLM, groups = 2) + 
              labs(title = "Weighted subgroup average effect") +
              theme(plot.title = element_text(hjust = 0.5))
          }
        })
      } else {
        return(NULL)
      }
    })
  }

  shinyApp(ui = ui, server = server)
}
