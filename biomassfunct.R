library(gridExtra)
library(openair)
library(dplyr)
library(ggplot2)

b <- as.matrix(read.table('coeff.txt',sep=","))
#' Allometric function for Pinus radiata
#'
#' @param x        value measured
#' @param variable variable measured. Use "D" for diameter [cm] at breast height;
#'                 "H" for total height [m] or "D2H" for diameter squared times
#'                 total height.
#' @param type     type of model. When variable is "D" or "H", only "nonlinear" 
#'                 approach is available. When variable is "D2H", both "linear"
#'                 and "nonlinear" approaches are avilable.
#' @param sigma    Variance-covariance matrix. Value equal to 1 assumes no 
#'                 covariance between the different components; 2 assumes
#'                 constant covariance between the different components, and 3
#'                 covariance depends on the different components.
#' @details The functions is part on our article that aims to propose allometric
#'  equations for P. radiata trees in Chile, ensuring compatibility between the 
#'  expected values for each component, but also between equation errors, 
#'  allowing for uncertainty estimation of biomass and carbon capture on pine 
#'  plantations.
#' @return
#' A data frame including the predicted value of each biomass component [kg] 
#' (Stem, Bark, Branches and Leaves) and their variance according to the 
#' selected model. Besides, the biomass expasion factor and its variance are 
#' returned.
#' @export
#' 
#' @examples
#' # Biomass components of a tree with a diameter at breast height of 20 cm.
#' Biomass(20, "D")
#' @author 
#' Sandoval, S; Montes, CR; Mena-Quijada, P; Acuña, E; Olmedo GF.
Biomass <- function(x, variable, type = "nonlinear", sigma = 1) {
  if(variable == "D" & type == "nonlinear" & sigma == 1){
    b <- as.matrix(read.table('coeff.txt',sep=",")[1,])
  } else if(variable == "D" & type == "nonlinear" & sigma == 2){
    b <- as.matrix(read.table('coeff.txt',sep=",")[2,])
  } else if(variable == "D" & type == "nonlinear" & sigma == 3){
    b <- as.matrix(read.table('coeff.txt',sep=",")[3,])
  } else if(variable == "H" & type == "nonlinear" & sigma == 1){
    b <- as.matrix(read.table('coeff.txt',sep=",")[4,])
  } else if(variable == "H" & type == "nonlinear" & sigma == 2){
    b <- as.matrix(read.table('coeff.txt',sep=",")[5,])
  } else if(variable == "H" & type == "nonlinear" & sigma == 3){
    b <- as.matrix(read.table('coeff.txt',sep=",")[6,])
  } else if(variable == "D2H" & type == "linear" & sigma == 1){
    b <- as.matrix(read.table('coeff.txt',sep=",")[7,])
  } else if(variable == "D2H" & type == "linear" & sigma == 2){
    b <- as.matrix(read.table('coeff.txt',sep=",")[8,])
  } else if(variable == "D2H" & type == "linear" & sigma == 3){
    b <- as.matrix(read.table('coeff.txt',sep=",")[9,])
  } else if(variable == "D2H" & type == "nonlinear" & sigma == 1){
    b <- as.matrix(read.table('coeff.txt',sep=",")[10,])
  } else if(variable == "D2H" & type == "nonlinear" & sigma == 2){
    b <- as.matrix(read.table('coeff.txt',sep=",")[11,])
  } else if(variable == "D2H" & type == "nonlinear" & sigma == 3){
    b <- as.matrix(read.table('coeff.txt',sep=",")[12,])
  }
  #Error Propagation
  
  prop_uncertainty_sum <- function(sA,sB,sC,sD) {  
    return (((sA^2)+(sB^2)+(sC^2)+(sD^2))**0.5)
  }
  
  
  prop_uncertainty_fraction <- function(f,A,Sa,B,Sb,COVab){
    return(f*(((A/Sa)^2)+((B/Sb)^2)-(2*COVab/(A*B)))**0.5)
  }
  
  if(type == "nonlinear"){
    stem      <- b[ 1]*x^b[ 2]
    bark      <- b[ 3]*x^b[ 4]
    branches  <- b[ 5]*x^b[ 6]
    leaves    <- b[ 7]*x^b[ 8] 
    tree      <- stem + bark + branches + leaves
    
    stemStd   <- b[ 9]*x^b[10]
    barkStd   <- b[11]*x^b[12]
    branchStd <- b[13]*x^b[14]
    leavesStd <- b[15]*x^b[16]
    treeStd   <- prop_uncertainty_sum(stemStd,barkStd,branchStd,leavesStd)
    
  } else if(type == "linear"){
    stem      <- b[ 1]+x*b[ 2]
    bark      <- b[ 3]+x*b[ 4]
    branches  <- b[ 5]+x*b[ 6]
    leaves    <- b[ 7]+x*b[ 8] 
    tree      <- stem + bark + branches + leaves
    
    stemStd   <- b[ 9]*x^b[10]
    barkStd   <- b[11]*x^b[12]
    branchStd <- b[13]*x^b[14]
    leavesStd <- b[15]*x^b[16]
    
    treeStd   <- prop_uncertainty_sum(stemStd,barkStd,branchStd,leavesStd)
  }
  
  #aca genera la matriz COV
  if (sigma==1){
    covs <- c(0,0,0,0)
  }else if (sigma==2 ){
    covs <- rep(b[17],4)
  }else if (sigma==3){
    covs <- c(b[17],b[18],b[20],b[23])
  }else{
    print("Parameter list error")
  }

  
  output <-data.frame(variable = x, Stem = stem, Bark = bark, Branches = branches, Leaves = leaves, Tree   = tree, 
                      StemStd = stemStd, BarkStd = barkStd, BranchesStd = branchStd, LeavesStd = leavesStd, TreeStd = treeStd)
  
  names(output)[1]<-variable
  
  return(output)
}
