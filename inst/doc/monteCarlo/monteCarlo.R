# load the micEconCES package
library( micEconCES )

# seed for random number generation
set.seed( 1234 )

# number of replications
nRep <- 1000

# number of observations
nObs <- 100

# rho
rho <- 1/3

# variance of the error term
uStdDev <- 1.5

# create data set with explanatory variables
cesData <- data.frame( x1 = rchisq( nObs, 10 ), x2 = rchisq( nObs, 10 ) )

# names of explanatory variables
xxNames <- c( "x1", "x2" )

# variable returns to scale
vrs <- TRUE # FALSE #

# coefficients
cesCoef <- c( gamma = 1, delta = 0.6, rho = rho, nu = 1.1 )[ 1:( 3 + vrs ) ]

# calculate deterministic endogenous variable
cesData$yd <- cesCalc( xNames = xxNames, data = cesData, coef = cesCoef )

# estimation methods
allMethods <- c( "Kmenta", "nls", "LM", "CG", "Newton", "BFGS",
   "Nelder-Mead", "SANN", "DE", "L-BFGS-B", "PORT" )

# objects to store the results
estCoef <- array( NA,
   dim = c( nRep, length( cesCoef ), length( allMethods ) ),
   dimnames = list( 1:nRep, names( cesCoef ), allMethods ) )
convergence <- estCoef[ , 1, ]
rss <- estCoef[ , 1, ]
rSquared <- estCoef[ , 1, ]
iterations <- estCoef[ , 1, ]

## start the monte carlo experiment
for( i in 1:nRep ) {
   cat( i, ": ", sep = "" )
   ptm <- proc.time()

   # adding noise to the endogenous variable
   repeat{
      cesData$ys <- cesData$yd + rnorm( nObs, sd = uStdDev )
      if( min( cesData$ys ) > 0 ) {
         break
      } else {
         cat( "#" )
      }
   }

   # estimate the model using different estimation methods
   for( method in allMethods ) {
      extraArgs <- list()
      if( method == "nls" ) {
         extraArgs <- list(
            control = nls.control( warnOnly = TRUE ) )
      } else if( method == "LM" ) {
         extraArgs <- list(
            control = nls.lm.control( maxiter = 250 ) )
      } else if( method == "Newton" ) {
         extraArgs <- list( iterlim = 250 )
      } else if( method %in% c( "BFGS", "L-BFGS-B" ) ) {
         extraArgs <- list( control = list( maxit = 250 ) )
      } else if( method == "CG" ) {
         extraArgs <- list(
            control = list( maxit = 1000, reltol = 1e-4, type = 2 ) )
      } else if( method == "SANN" ) {
         extraArgs <- list( control = list( maxit = 50000 ) )
      } else if( method == "DE" ) {
         extraArgs <- list(
            control = DEoptim.control( trace = FALSE, itermax = 1000 ) )
      }
      allArgs <- c( list( yName = "ys", xNames = xxNames, data = cesData,
         method = method, vrs = vrs ), extraArgs ) 
      cesResult <- try( do.call( "cesEst", allArgs ) )
      if( class( cesResult )[1] != "try-error" ) {
         # store the estimated coefficients
         estCoef[ i, , method ] <- coef( cesResult )
         # store if the estimation has converged
         if( !is.null( cesResult$convergence ) ) {
            convergence[ i, method ] <- cesResult$convergence
         }
         # sum of squared residuals
         rss[ i, method ] <- cesResult$rss
         # R-squared values
         rSquared[ i, method ] <- summary( cesResult )$r.squared
         # number of iterations
         if( !is.null( cesResult$iter ) ) {
            iterations[ i, method ] <- sum( cesResult$iter )
         }
      }
   }

   ptmNew <- proc.time()
   cat( ptmNew - ptm, "\n" )
   ptm <- ptmNew
}

########### calculate summary results ##############
# differences between the estimated and the true coefficients
diffCoef <- estCoef - aperm(
   array( cesCoef, dim = c( length( cesCoef ), length( allMethods ), nRep ) ),
   c( 3, 1, 2 ) )

# elasticities of substitution and difference between estimates and true value
estSigma <- 1 / ( 1 + estCoef[ ,"rho", ] )
diffSigma <- estSigma - 1 / ( 1 + rho )

# function to calculate summary results of the Monte Carlo simulation
# depending on the selection of replications
calcMcResults <- function( repSelect ) {

   result <- list()

   # biases of the estimated coefficients and elasticity of substitution
   result$bias <- colMeans( diffCoef[ repSelect, , ] )
   # all.equal( bias, colMeans( estCoef ) - matrix( cesCoef, nrow = length( cesCoef ), ncol = length( allMethods ) ) )
   result$bias <- rbind( result$bias,
      sigma = colMeans( diffSigma[ repSelect, ] ) )

   # median deviation of estimated coef. and elast. of subst. from their true values
   result$devMed <- colMedians( diffCoef[ repSelect, , ] )
   result$devMed <- rbind( result$devMed,
      sigma = colMedians( diffSigma[ repSelect, ] ) )

   # root mean squared errors of the estimated coefficients and elasticity of substitution
   result$rmse <- sqrt( colSums( diffCoef[ repSelect, , ]^2 ) / nRep )
   result$rmse <- rbind( result$rmse,
      sigma =  colSums( diffSigma[ repSelect, ]^2 ) / nRep )

   # mean absolute deviations
   result$mad <- colMeans( abs( diffCoef[ repSelect, , ] ) )
   result$mad <- rbind( result$mad,
      sigma = colMeans( abs( diffSigma[ repSelect, ] ) ) )

   # median absolute deviations
   result$adMed <- colMedians( abs( diffCoef[ repSelect, , ] ) )
   result$adMed <- rbind( result$adMed,
      sigma = colMedians( abs( diffSigma[ repSelect, ] ) ) )

   # mean RSS
   result$rssMean <- colMeans( rss[ repSelect, ] )

   # mean R-squared values
   result$rSquaredMean <- colMeans( rSquared[ repSelect, ] )

   return( result )
}

# summary results of *all* replications
resultAll <- calcMcResults( 1:nRep )

# summary results of replications without errors (in any method)
resultNoErr <- calcMcResults( rowSums( is.na( rss ) ) == 0 )

# summary results of replications without errors or non-convergence (in any method)
resultConv <- calcMcResults(
   rowSums( is.na( rss ) | ( !convergence & !is.na( convergence ) ) ) == 0 )


########### create tables for the paper ##############
# general results
tabGeneral <- data.frame( nNoConv =
   colSums( is.na( rss ) | ( !convergence & !is.na( convergence ) ) ) )
tabGeneral$nConv <-
   colSums( !is.na( rss ) & ( convergence | is.na( convergence ) ) )
tabGeneral$rssAll <- resultAll$rssMean
tabGeneral$rssConv <- resultConv$rssMean


############ write tables to disk ##############
library( xtable )
# general results
xTabGeneral <- xtable( tabGeneral, digits = c( rep( 0, 3 ), rep( 7, 2 ) ),
   align = c( "l", rep( "r", 4 ) ) )
print( xTabGeneral, file = "../tables/mcGeneral.tex", floating = FALSE )

# bias
xBias <- xtable( t( resultAll$bias ), digits = rep( 5, 6 ),
   align = c( "l", rep( "r", 5 ) ) )
print( xBias, file = "../tables/mcBias.tex", floating = FALSE )

# root mean square error
xRmse <- xtable( t( resultAll$rmse ), digits = rep( 5, 6 ),
   align = c( "l", rep( "r", 5 ) ) )
print( xRmse, file = "../tables/mcRmse.tex", floating = FALSE )
