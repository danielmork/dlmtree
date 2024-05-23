\donttest{
  
  # The first three examples are for one lagged exposure
  
  
  # treed distributed lag model (TDLM)
  # binary outcome with logit link
  
  D <- sim.tdlmm(sim = "A", mean.p = 0.5, n = 1000)
  fit_tdlm <- dlmtree(y ~ .,
                      data = D$dat,
                      exposure.data = D$exposures[[1]],
                      dlm.type = "linear",
                      family = "logit",
                      binomial.size = 1)
  
  # summarize results
  s_fit_tdlm <- summary(fit_tdlm)
  s_fit_tdlm
  
  # plot results
  plot(s_fit_tdlm)
  
  
  
  # Treed distributed lag nonlinear model (TDLNM)
  # Gaussian regression model
  D <- sim.tdlnm(sim = "A", error.to.signal = 1)
  fit_tdlnm <- dlmtree(formula = y ~ .,
                       data = D$dat,
                       exposure.data = D$exposures,
                       dlm.type = "nonlinear",
                       family = "gaussian")
  
  # summarize results
  s_fit_tdlnm <- summary(fit_tdlnm)
  s_fit_tdlnm
  
  # plot results
  plot(s_fit_tdlnm)
  
  
  
  # Heterogenious TDLM (HDLM), similar to first example but with heterogenious exposure response
  D <- sim.hdlmm(sim = "B", n = 1000)
  fit_hdlm <- dlmtree(y ~ .,
                      data = D$dat,
                      exposure.data = D$exposures,
                      dlm.type = "linear",
                      family = "gaussian",
                      het = TRUE)
  
  # summarize results
  s_fit_hdlm <- summary(fit_hdlm)
  s_fit_hdlm
  
  # shiny app for HDLM
  if (interactive()) {
    shiny(fit_hdlm)
  }
  
  
  
  # The next two examples are for a mixture (or multivariate) exposure 
  
  
  # Treed distributed lag mixture model (TDLMM)
  # Model for mixutre (or multivariate) lagged exposures
  # with a homogenious exposure-time-response function
  D <- sim.tdlmm(sim = "B", error = 25, n = 1000)
  fit_tdlmm <- dlmtree(y ~ .,
                       data = D$dat, exposure.data = D$exposures,
                       mixture.interactions = "noself",
                       dlm.type = "linear", family = "gaussian",
                       mixture = TRUE)
  
  # summarize results
  s_fit_tdlmm <- summary(fit_tdlmm)
  
  # plot the marginal exposure-response for one exposure
  plot(s_fit_tdlmm, exposure1 = "e1")
  
  # plot exposure-response surface
  plot(s_fit_tdlmm, exposure1 = "e1", exposure2 = "e2")
  
  
  
  # heterogenious version of TDLMM
  D <- sim.hdlmm(sim = "D", n = 1000)
  fit_hdlmm <- dlmtree(y ~ .,
                       data = D$dat,
                       exposure.data = D$exposures,
                       dlm.type = "linear",
                       family = "gaussian",
                       mixture = TRUE,
                       het = TRUE)
  
  # summarize results
  s_fit_hdlmm <- summary(fit_hdlmm)
  s_fit_hdlmm
  
  # summarize results
  if (interactive()) {
    shiny(s_fit_hdlmm)
  }
  
  
}
