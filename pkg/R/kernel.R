parzen <- function(k, window.type = c("acf", "spectrum"), M = 10) {
  window.type <- match.arg(window.type)
  if(window.type=="acf"){
    result <- (abs(k)<=M/2)*(1-6*(k/M)^2+6*(abs(k)/M)^3)+((M/2<abs(k))&(abs(k)<=M))*2*(1-abs(k)/M)^3
  }
  if(window.type=="spectrum"){
    Index <- k==0
    result <- 3/8/pi/M^3*(sin(k*M/4)/(1/2*sin(k/2)))^4*(1-2/3*(sin(k/2))^2)
    result[Index] <- 3/8*M/pi
  }
  return(result)
}

daniell <- function(k, window.type = c("acf", "spectrum"), M = 10) {
  window.type <- match.arg(window.type)
  if(window.type=="acf") {
    Index <- k==0
    result <- (abs(k)<=M)*sin(pi*k/M)/(k/M*pi)
    result[Index] <- 1
  }
  if(window.type=="spectrum"){
  result <- (abs(k)<= pi/M)*M/2/pi
  }
  return(result)
}

bartlett <- function(k, window.type = c("acf", "spectrum"), M = 10) {
  window.type <- match.arg(window.type)
  if(window.type=="acf"){
    result <- (abs(k)<=M)*(1-abs(k)/M)
  }
  if(window.type=="spectrum"){
    Index <- k==0
    result <- 1/2/pi/M*(sin(k*M/2)/sin(k/2))^2
    result[Index] <- M/2/pi
  }
  return(result)
}

rectangular <- function(k, window.type = c("acf", "spectrum"), M = 10) {
  window.type <- match.arg(window.type)
  if(window.type=="acf"){
    result <- abs(k)<=M
  }
  if(window.type=="spectrum"){
    Index <- k==0
    result <- 1/2/pi*sin(k*(M+1/2))/sin(k/2)
    result[Index] <- (M+1/2)/pi
  }
  return(result)
}
