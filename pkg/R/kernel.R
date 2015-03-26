parzen <- function(k,window.type="acf",M=10) {
if(window.type=="acf")
return((abs(k)<=M/2)*(1-6*(k/M)^2+6*(abs(k)/M)^3)+((M/2<abs(k))&(abs(k)<=M))*2*(1-abs(k)/M)^3)
if(window.type=="spectrum"){
Index <- k==0
erg <- 3/8/pi/M^3*(sin(k*M/4)/(1/2*sin(k/2)))^4*(1-2/3*(sin(k/2))^2)
erg[Index] <- 3/8*M/pi
return(erg)}
warning("Input type is unkonwn. Using identity as kernel.")
return(k)
}

parzen <- function(k,window.type="acf",M=10) {
if(window.type=="acf")
return((abs(k)<=M/2)*(1-6*(k/M)^2+6*(abs(k)/M)^3)+((M/2<abs(k))&(abs(k)<=M))*2*(1-abs(k)/M)^3)
if(window.type=="spectrum"){
Index <- k==0
erg <- 3/8/pi/M^3*(sin(k*M/4)/(1/2*sin(k/2)))^4*(1-2/3*(sin(k/2))^2)
erg[Index] <- 3/8*M/pi
return(erg)}
warning("Input type is unkonwn. Using identity as kernel.")
return(k)
}

bartlett <- function(k,window.type="acf",M=10) {
if(window.type=="acf")
return((abs(k)<=M)*(1-abs(k)/M))
if(window.type=="spectrum"){
Index <- k==0
erg <- 1/2/pi/M*(sin(k*M/2)/sin(k/2))^2
erg[Index] <- M/2/pi
return(erg)}
warning("Input type is unkonwn. Using identity as kernel.")
return(k)
}

rectangular <- function(k,window.type="acf",M=10) {
if(window.type=="acf")
return((abs(k)<=M))
if(window.type=="spectrum"){
Index <- k==0
erg <- 1/2/pi*sin(k*(M+1/2))/sin(k/2)
erg[Index] <- (M+1/2)/pi
return(erg)}
warning("Input type is unkonwn. Using identity as kernel.")
return(k)
}



