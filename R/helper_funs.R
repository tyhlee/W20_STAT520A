cols <- c("naive"="red",
          "orbit_alt" = "blue",
          "orbit_mix" = "brown")

lvls <- c("naive","orbit_alt","orbit_mix")

text_size <- 20

sphereplot3 <- function (x,xx, xxx, col = c("blue","red","brown")) 
{
  if (is.null(col)) 
    col <- rep(1, dim(x)[1])
  rgl::open3d()
  y1 <- x[, 1]
  y2 <- x[, 2]
  y3 <- x[, 3]
  rgl::points3d(y1, y2, y3, col = col[1], radius = 1)
  rgl::points3d(xx[,1], xx[,2], xx[,3], col = col[2], radius = 1)
  rgl::points3d(xxx[,1], xxx[,2], xxx[,3], col = col[3], radius = 1)
  rgl::spheres3d(0, 0, 0, lit = FALSE, color = "white")
  rgl::spheres3d(0, 0, 0, radius = 1, lit = FALSE, color = "black", 
                 front = "lines")
}

# temperature function
# earth temp: 
# http://www.physicalgeography.net/fundamentals/7m.html
# middle: 35C => approx 310K
# poles : -40C => approx 230K
temp <- function(z){
  310 - 80*z^2
}

mse <- function(a,b){
  mean((a-b)^2)
}

# compute the radius of a x-y circle at z
sphere_radius_z <- function(z,r=1){
  tmp_r <- sqrt(r^2-z^2)
  tmp_r
}

# randomly sample a rotation matrix from SO(3)
sample_rot <- function(n=1){
  random_rot <- ruars(n,rangle=rhaar,space="SO3")
  class(random_rot) <- "matrix"
  matrix(random_rot,nrow=3,byrow = F)
}

# transformation of the beta distribution
tranformed_density <- function(a=310,b=-80,alpha=1/2,beta=1){
  ttt <- seq(0.001,0.999,by=0.001)
  transformed_ttt <- a+b*ttt
  den <- (alpha-1)*log( (transformed_ttt-a)/b ) + 
    (beta-1)*(1 - (transformed_ttt-a)/b)+
    log(abs(1/b)) -
    log(beta(alpha,beta))
  return(data.frame(x=transformed_ttt,
                    y=exp(den),
                    type='true'))
}

