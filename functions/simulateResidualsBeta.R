simulateResidualsBeta <- function(model){

# function to simulate residuals beta distribution  
beta.sim <- function(gam.object, rep){
    
  
  
    mu <- gam.object$fitted.values
    phi <- model$coefficients$precision
    
    replicate(rep, rbeta(length(mu), 
                         shape1 = (mu*phi), 
                         shape2 = (phi - (mu*phi))))
    
  }
  
gam.sim = beta.sim(model, 250)

gam.res <- createDHARMa(simulatedResponse = gam.sim,
                        observedResponse = model$model[,1],
                        fittedPredictedResponse = predict(model, 
                                                          type = "response"),
                        integerResponse = TRUE)
gam.res$refit=FALSE
gam.res$simulatedResponse = gam.sim

return(gam.res)
}
