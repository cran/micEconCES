# load the micEconCES package
library( "micEconCES" )

# load the data set
data( "GermanIndustry" )

# remove years 1973 - 1975 because of economic disruptions (see Kemfert 1998)
GermanIndustry <- subset( GermanIndustry, year < 1973 | year > 1975, )

# time trend (starting with 0)
GermanIndustry$time <- GermanIndustry$year - 1960

# names of inputs
xNames1 <-  c( "K", "E", "A" )
xNames2 <-  c( "K", "A", "E" )
xNames3 <-  c( "E", "A", "K" )


########## Grid Search for Rho_1 and/or Rho ##############
rhoVec <- c( seq( -1, 1, 0.1 ), seq( 1.2, 4, 0.2 ), seq( 4.4, 14, 0.4 ) )

## BFGS, grid search for rho_1 and rho
cesBfgsGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridRho1 )
plot( cesBfgsGridRho1 )
cesBfgsGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridRho2 )
plot( cesBfgsGridRho2 )
cesBfgsGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridRho3 )
plot( cesBfgsGridRho3 )

# BFGS with grid search estimates as starting values
cesBfgsGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesBfgsGridRho1 ),
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridStartRho1 )
cesBfgsGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesBfgsGridRho2 ), 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridStartRho2 )
cesBfgsGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesBfgsGridRho3 ), 
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgsGridStartRho3 )


## Levenberg-Marquardt, grid search for rho1 and rho
cesLmGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridRho1 )
plot( cesLmGridRho1 )
cesLmGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridRho2 )
plot( cesLmGridRho2 )
cesLmGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridRho3 )
plot( cesLmGridRho3 )

# LM with grid search estimates as starting values
cesLmGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesLmGridRho1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridStartRho1 )
cesLmGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesLmGridRho2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridStartRho2 )
cesLmGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesLmGridRho3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmGridStartRho3 )


## Newton-type, grid search for rho_1 and rho
cesNewtonGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "Newton", iterlim = 500 )
summary( cesNewtonGridRho1 )
plot( cesNewtonGridRho1 )
cesNewtonGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridRho2 )
plot( cesNewtonGridRho2 )
cesNewtonGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, 
   method = "Newton", iterlim = 500 )
summary( cesNewtonGridRho3 )
plot( cesNewtonGridRho3 )

# Newton-type with grid search estimates as starting values
cesNewtonGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesNewtonGridRho1 ),
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridStartRho1 )
cesNewtonGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesNewtonGridRho2 ),
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridStartRho2 )
cesNewtonGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesNewtonGridRho3 ),
   method = "Newton", iterlim = 500, check.analyticals = FALSE )
summary( cesNewtonGridStartRho3 )


## PORT, grid search for rho1 and rho
cesPortGridRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridRho1 )
plot( cesPortGridRho1 )
cesPortGridRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridRho2 )
plot( cesPortGridRho2 )
cesPortGridRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   rho1 = rhoVec, rho = rhoVec, returnGridAll = TRUE, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridRho3 )
plot( cesPortGridRho3 )

# PORT with grid search estimates as starting values
cesPortGridStartRho1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   start = coef( cesPortGridRho1 ),, method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridStartRho1 )
cesPortGridStartRho2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry, 
   start = coef( cesPortGridRho2 ), method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridStartRho2 )
cesPortGridStartRho3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry, 
   start = coef( cesPortGridRho3 ), method = "PORT",
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortGridStartRho3 )


save.image( "kemfert98_grid.RData" )


