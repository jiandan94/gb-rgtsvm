ParameterGrid <- function(x){
  xName <- names(x)
  p <- length(xName)
  xLength <- rep(0, p)
  for(i in 1:length(xName)){# obtain the number of values for each variable
    xLength[i] <- length(x[[xName[i]]])
  }
  n <- prod(xLength)
  
  ParameterGridMatrix <- matrix(0, nrow = n, ncol = p)
  ParameterGridMatrix[,1] <- rep(x[[xName[1]]], rep(prod(xLength[2:p]), xLength[1]))
  ParameterGridMatrix[,p] <- rep(x[[xName[p]]], prod(xLength[1:(p - 1)]))
  if(p > 2){
    for(i in 2:(p - 1)){
      T1 <- prod(xLength[(i + 1):p])
      T2 <- prod(xLength[1:(i - 1)])
      ParameterGridMatrix[,i] <- rep(rep(x[[xName[i]]], rep(T1, xLength[i])), T2)
    }
  }
  
  colnames(ParameterGridMatrix) <- xName
  return(ParameterGridMatrix)
}