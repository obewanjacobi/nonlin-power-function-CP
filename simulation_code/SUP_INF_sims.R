# Jacob Townson
# Simulating the open square around phi to observe examples to get ideas for a 
# proof that the sup(f_i(phi)-f_i(theta_0))^2 and inf(/sum(f_i(phi)-f_i(theta_0))^2)
# can be found at the corners of our open square inscribed around the open ball around phi
# with radius of r_theta

## First we built the sup and inf functions that will find the max and min (since
## in our case they're equivalent) of the squared differences for each case of the 
## A' assumption.

sup_func <- function(x_i, C, b, C_0, b_0){
  sup = max((C*x_i^b - C_0*x_i^b_0)^2)
  return(sup)
}

inf_func <- function(x_i, C, b, C_0, b_0){
  inf = min(sum(C*x_i^b - C_0*x_i^b_0)^2)
  return(inf)
}

## Next up we set up a function to define the open square around phi. It will 
## return the four corners of the open square in the parameter field of (C,b).
## Note that phi here should be entered as a point in the parameter field of
## (C,b) as well, eg (20,.5).

op_square <- function(phi, r_theta){
  t_l = phi + c(r_theta, -r_theta)
  t_r = phi + c(r_theta, r_theta)
  b_l = phi + c(-r_theta, -r_theta)
  b_r = phi + c(-r_theta, r_theta)
  
  square = list("top_left" = t_l,
                "top_right" = t_r,
                "bottom_left" = b_l,
                "bottom_right" = b_r)
  
  return(square)
}

## The next step is to actually try and make some simulations. Users should save
## interesting seeds that give unexpected outcomes.

set.seed(1142020)

x_i = runif(100,0,20)
phi = c(10,.3)
r_theta = .2           # Something to think about here, r_theta needs to be pretty
C_0 = 15                 # small since 0 < b < 1. 
b_0 = .75

if(phi[2] - r_theta < 0 | phi[2] + r_theta > 1){
  print("WARNING: YOUR OPEN SQUARE IS OUT OF BOUNDS OF THE POSSIBLE b VALUES")
}
if(phi[1] - r_theta < 0){
  print("WARNING: YOUR OPEN SQUARE IS OUT OF BOUNDS OF THE POSSIBLE C VALUES")
}

corners = op_square(phi, r_theta)
C = runif(100000, corners[[3]][1], corners[[1]][1])
b = runif(100000, corners[[1]][2], corners[[2]][2])

sup_outcomes = c()
inf_outcomes = c()
for(i in 1:length(C)){
  sup = sup_func(x_i, C[i], b[i], C_0, b_0)
  inf = inf_func(x_i, C[i], b[i], C_0, b_0)
  
  sup_outcomes[i] = sup
  inf_outcomes[i] = inf
}

cornerSup_outcomes = c()
cornerInf_outcomes = c()
for(i in 1:4){
  sup = sup_func(x_i, corners[[i]][1], corners[[i]][2], C_0, b_0)
  inf = inf_func(x_i, corners[[i]][1], corners[[i]][2], C_0, b_0)
  
  cornerSup_outcomes[i] = sup
  cornerInf_outcomes[i] = inf
}


## Compare results:

### sups

max(sup_outcomes)
max(cornerSup_outcomes)

if(max(sup_outcomes) < max(cornerSup_outcomes)){
  print("The corners contain the sup.")
}else{
  print("NOTE: The sup was not in one of the corners.")
}

### infs

min(inf_outcomes)
min(cornerInf_outcomes)

if(min(inf_outcomes) > min(cornerInf_outcomes)){
  print("The corners contain the inf.")
}else{
  print("NOTE: The inf was not in one of the corners.")
}

# Visualizing in 3D

require(rgl)

infF <- function(C, b){
  inf = min(sum((C*x_i^b - C_0*x_i^b_0)^2))
  return(inf)
}
supF <- function(C, b){
  sup = max((C*x_i^b - C_0*x_i^b_0)^2)
  return(sup)
}
infF=Vectorize(infF)
persp3d(infF,xlim=c(phi[1]-(r_theta),phi[1]+(r_theta)),ylim=c(phi[2]-(r_theta),phi[2]+(r_theta)),
        #zlim=c(0,max(sup_outcomes)),
        n=100,xlab="C",ylab="b",zlab = "inf")
supF=Vectorize(supF)
persp3d(supF,xlim=c(phi[1]-(r_theta+.1),phi[1]+(r_theta+.1)),ylim=c(phi[2]-(r_theta+.1),phi[2]+(r_theta+.1)),
        #zlim=c(min(cornerInf_outcomes),max(cornerSup_outcomes)),
        n=100,xlab="C",ylab="b",zlab = "sup")











