#' Not my function

plot.progress <- function(...)	{
  vectOfBar <- c(...)*100
  numOfBar <- length(vectOfBar)
  plot(c(0,100), c(0,numOfBar), type='n', xlab='', ylab='', yaxt='n', mar=c(3,3,3,3))
  for(i in 1:numOfBar) {
    rect(0, 0.1+i-1, vectOfBar[i], 0.9+i-1, col=rainbow(numOfBar)[i])
    text(0.5, 0.5+i-1, paste('Running!!) ', i, ': ', round(vectOfBar[i],2), '%', sep=''), adj=0)
  }
  title('Progress...')
}
