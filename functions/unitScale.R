# convert numeric data to scale (default 0-1). 
#Primarily for use in colour ramps for plotting.
unitScale <- function(x, custMin = NULL, custMax=NULL){
  if(is.null(custMax)){custMax = max(x, na.rm=TRUE)}
  if(is.null(custMin)){custMin = min(x, na.rm=TRUE)}
  (x - custMin) / (custMax - custMin)
}