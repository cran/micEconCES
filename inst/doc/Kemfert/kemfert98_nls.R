# load the micEconCES package
library( "micEconCES" )

# load the data set
data( "GermanIndustry" )

# remove years 1973 - 1975 because of economic disruptions (see Kemfert 1998)
GermanIndustry <- subset( GermanIndustry, year < 1973 | year > 1975, )

# add a time trend (starting with 0)
GermanIndustry$time <- GermanIndustry$year - 1960

# names of inputs
xNames1 <-  c( "K", "E", "A" )
xNames2 <-  c( "K", "A", "E" )
xNames3 <-  c( "E", "A", "K" )

################# econometric estimation with cesEst ##########################

## Nelder-Mead
cesNm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "NM", control = list( maxit = 5000 ) )
summary( cesNm1 )
cesNm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "NM", control = list( maxit = 5000 ) )
summary( cesNm2 )
cesNm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "NM", control = list( maxit = 5000 ) )
summary( cesNm3 )

## Simulated Annealing
cesSann1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "SANN", control = list( maxit = 2e6 ) )
summary( cesSann1 )
cesSann2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "SANN", control = list( maxit = 2e6 ) )
summary( cesSann2 )
cesSann3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "SANN", control = list( maxit = 2e6 ) )
summary( cesSann3 )

## BFGS
cesBfgs1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgs1 )
cesBfgs2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgs2 )
cesBfgs3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "BFGS", control = list( maxit = 5000 ) )
summary( cesBfgs3 )

## L-BFGS-B
cesBfgsCon1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", control = list( maxit = 5000 )  )
summary( cesBfgsCon1 )
cesBfgsCon2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", control = list( maxit = 5000 )  )
summary( cesBfgsCon2 )
cesBfgsCon3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "L-BFGS-B", control = list( maxit = 5000 )  )
summary( cesBfgsCon3 )

## Levenberg-Marquardt
cesLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm1 )
cesLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm2 )
cesLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm3 )

## Levenberg-Marquardt, multiplicative error term
cesLm1Me <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry, 
   multErr = TRUE, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm1Me )
cesLm2Me <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   multErr = TRUE, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm2Me )
cesLm3Me <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   multErr = TRUE, control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLm3Me )

## Newton-type
cesNewton1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "Newton", iterlim = 500 )
summary( cesNewton1 )
cesNewton2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "Newton", iterlim = 500 )
summary( cesNewton2 )
cesNewton3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "Newton", iterlim = 500 )
summary( cesNewton3 )

## PORT
cesPort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort1 )
cesPort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort2 )
cesPort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort3 )

## PORT, multiplicative error
cesPort1Me <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", multErr = TRUE,
   control = list( eval.max = 2000, iter.max = 2000 ) )
summary( cesPort1Me )
cesPort2Me <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", multErr = TRUE,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort2Me )
cesPort3Me <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", multErr = TRUE,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPort3Me )

## DE
cesDe1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 500,
   itermax = 1e4 ) )
summary( cesDe1 )
cesDe2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 500,
   itermax = 1e4 ) )
summary( cesDe2 )
cesDe3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 500,
   itermax = 1e4 ) )
summary( cesDe3 )

## nls
cesNls1 <- try( cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   vrs = TRUE, method = "nls" ) )
cesNls2 <- try( cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   vrs = TRUE, method = "nls" ) )
cesNls3 <- try( cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   vrs = TRUE, method = "nls" ) )

## NM - Levenberg-Marquardt
cesNmLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   start = coef( cesNm1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesNmLm1 )
cesNmLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   start = coef( cesNm2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesNmLm2 )
cesNmLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   start = coef( cesNm3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesNmLm3 )

## SANN - Levenberg-Marquardt
cesSannLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   start = coef( cesSann1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesSannLm1 )
cesSannLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   start = coef( cesSann2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesSannLm2 )
cesSannLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   start = coef( cesSann3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesSannLm3 )

## DE - Levenberg-Marquardt
cesDeLm1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   start = coef( cesDe1 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesDeLm1 )
cesDeLm2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   start = coef( cesDe2 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesDeLm2 )
cesDeLm3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   start = coef( cesDe3 ),
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesDeLm3 )

## NM - PORT
cesNmPort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesNm1 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesNmPort1 )
cesNmPort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesNm2 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesNmPort2 )
cesNmPort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesNm3 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesNmPort3 )

## SANN - PORT
cesSannPort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesSann1 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesSannPort1 )
cesSannPort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesSann2 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesSannPort2 )
cesSannPort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesSann3 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesSannPort3 )

## DE - PORT
cesDePort1 <- cesEst( "Y", xNames1, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesDe1 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesDePort1 )
cesDePort2 <- cesEst( "Y", xNames2, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesDe2 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesDePort2 )
cesDePort3 <- cesEst( "Y", xNames3, tName = "time", data = GermanIndustry,
   method = "PORT", start = coef( cesDe3 ),
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesDePort3 )



# save the workspace
save.image( "kemfert98_nls.RData" )

