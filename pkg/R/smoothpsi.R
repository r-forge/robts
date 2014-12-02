smoothpsi <- function(x,k=2,l=3) {
a <- (2*k^2*l^2)/(k-l)^3
b <- -(l^3+k*l^2+4*k^2*l)/(k-l)^3
d <- (2*l^2+2*k*l+2*k^2)/(k-l)^3
e <- -(l+k)/(k-l)^3
return(x*(abs(x)<=k)+sign(x)*(a+b*abs(x)+d*x^2+e*abs(x)^3)*((abs(x)>k)&(abs(x)<=l)))  # super-weights
}