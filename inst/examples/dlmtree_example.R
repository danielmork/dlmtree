\donttest{

  # The first three examples are for one lagged exposure


  # treed distributed lag model (TDLM)
  # binary outcome with logit link

  D <- sim.tdlmm(sim = "A", mean.p = 0.5, n = 1000)
  tdlm.fit <- dlmtree(y ~ .,
                      data = D$dat,
                      exposure.data = D$exposures[[1]],
                      dlm.type = "linear",
                      family = "logit",
                      binomial.size = 1)

  # summarize results
  tdlm.sum <- summary(tdlm.fit)
  tdlm.sum

  # plot results
  plot(tdlm.sum)



  # Treed distributed lag nonlinear model (TDLNM)
  # Gaussian regression model
  D <- sim.tdlnm(sim = "A", error.to.signal = 1)
  tdlnm.fit <- dlmtree(formula = y ~ .,
                       data = D$dat,
                       exposure.data = D$exposures,
                       dlm.type = "nonlinear",
                       family = "gaussian")

  # summarize results
  tdlnm.sum <- summary(tdlnm.fit)
  tdlnm.sum

  # plot results
  plot(tdlnm.sum)



  # Heterogeneous TDLM (HDLM), similar to first example but with heterogeneous exposure response
  D <- sim.hdlmm(sim = "B", n = 1000)
  hdlm.fit <- dlmtree(y ~ .,
                      data = D$dat,
                      exposure.data = D$exposures,
                      dlm.type = "linear",
                      family = "gaussian",
                      het = TRUE)

  # summarize results
  hdlm.sum <- summary(hdlm.fit)
  hdlm.sum

  # shiny app for HDLM
  if (interactive()) {
    shiny(hdlm.fit)
  }



  # The next two examples are for a mixture (or multivariate) exposure


  # Treed distributed lag mixture model (TDLMM)
  # Model for mixutre (or multivariate) lagged exposures
  # with a homogenious exposure-time-response function
  D <- sim.tdlmm(sim = "B", error = 25, n = 1000)
  tdlmm.fit <- dlmtree(y ~ .,
                       data = D$dat, exposure.data = D$exposures,
                       mixture.interactions = "noself",
                       dlm.type = "linear", family = "gaussian",
                       mixture = TRUE)

  # summarize results
  tdlmm.sum <- summary(tdlmm.fit)

  # plot the marginal exposure-response for one exposure
  plot(tdlmm.sum, exposure1 = "e1")

  # plot exposure-response surface
  plot(tdlmm.sum, exposure1 = "e1", exposure2 = "e2")



  # heterogeneous version of TDLMM
  D <- sim.hdlmm(sim = "D", n = 1000)
  hdlmm.fit <- dlmtree(y ~ .,
                       data = D$dat,
                       exposure.data = D$exposures,
                       dlm.type = "linear",
                       family = "gaussian",
                       mixture = TRUE,
                       het = TRUE)

  # summarize results
  hdlmm.sum <- summary(hdlmm.fit)
  hdlmm.sum

  # summarize results
  if (interactive()) {
    shiny(hdlmm.fit)
  }


}
