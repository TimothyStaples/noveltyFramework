gamPredInt <- function(model, quan, xvar, yvar){
  
SSE <- sum(model$residuals^2)
MSE <- SSE / nrow(model$model)
SSsd <- sd(model$residuals)
 
x <- model$model[,xvar]
y <- model$model[,yvar]

t.quantile <- qt(quan, nrow(model$model))

prediction <- predict(model)

return(prediction + (t.quantile * sqrt(MSE * (1 + (1/length(y) + ((x - mean(x)^2) / sum((x - mean(x)^2))))))))
}