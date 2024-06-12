# function provides numeric centers of numeric -> factor conversion with
# the "cut" function
cutCenter <- function(x, ...){
  
  cutX = cut(x, ...)
  
  leftNum = as.numeric(substr(cutX, 
                              regexpr("\\[|\\(", cutX)+1, 
                              regexpr(",", cutX)-1))
  
  rightNum = as.numeric(substr(cutX, 
                               regexpr(",", cutX)+1, 
                               regexpr("\\]", cutX)-1))
  
  return(rowMeans(cbind(leftNum, rightNum)))
}