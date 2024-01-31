colourBlend <- function(c1, c2){
  
  rgb1 <- col2rgb(c1)
  rgb2 <- col2rgb(c2)
  
  # Calculate composite color using Multiply Blend Mode
  compCol <- sapply(1:3, function(n){
    (rgb1[n] * rgb2[n]) / 255
  }) / 255
 
  return(rgb(compCol[1], compCol[2], compCol[3]))
  
}