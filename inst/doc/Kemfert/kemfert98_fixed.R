# load the micEconCES package
library( "micEconCES" )

# load the data set
data( "GermanIndustry" )

# remove years 1973 - 1975 because of economic disruptions (see Kemfert 1998)
GermanIndustry <- subset( GermanIndustry, year < 1973 | year > 1975, )

# time trend (starting with 0)
GermanIndustry$time <- GermanIndustry$year - 1960

# removing technological progress using the lambdas of Kemfert (1998)
# (we can do this, because the model has constant returns to scale)
GermanIndustry$K1 <- GermanIndustry$K * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$E1 <- GermanIndustry$E * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$A1 <- GermanIndustry$A * exp( 0.0222 * GermanIndustry$time )
GermanIndustry$K2 <- GermanIndustry$K * exp( 0.0069 * GermanIndustry$time )
GermanIndustry$E2 <- GermanIndustry$E * exp( 0.0069 * GermanIndustry$time )
GermanIndustry$A2 <- GermanIndustry$A * exp( 0.0069 * GermanIndustry$time )
GermanIndustry$K3 <- GermanIndustry$K * exp( 0.00641 * GermanIndustry$time )
GermanIndustry$E3 <- GermanIndustry$E * exp( 0.00641 * GermanIndustry$time )
GermanIndustry$A3 <- GermanIndustry$A * exp( 0.00641 * GermanIndustry$time )

# names of inputs
xNames1 <-  c( "K1", "E1", "A1" )
xNames2 <-  c( "K2", "A2", "E2" )
xNames3 <-  c( "E3", "A3", "K3" )


############# estimation with lambda, rho_1, and rho fixed #####################

## Nelder-Mead, lambda, rho_1, and rho fixed
cesNmFixed1 <- cesEst( "Y", xNames1, data = GermanIndustry,
   method = "NM", rho1 = 0.5300, rho = 0.1813,
   control = list( maxit = 5000 ) )
summary( cesNmFixed1 )
cesNmFixed2 <- cesEst( "Y", xNames2, data = GermanIndustry,
   method = "NM", rho1 = 0.2155, rho = 1.1816,
   control = list( maxit = 5000 ) )
summary( cesNmFixed2 )
cesNmFixed3 <- cesEst( "Y", xNames3, data = GermanIndustry,
   method = "NM", rho1 = 1.3654, rho = 5.8327,
   control = list( maxit = 5000 ) )
summary( cesNmFixed3 )

## BFGS, lambda, rho_1, and rho fixed
cesBfgsFixed1 <- cesEst( "Y", xNames1, data = GermanIndustry,
   method = "BFGS", rho1 = 0.5300, rho = 0.1813,
   control = list( maxit = 5000 ) )
summary( cesBfgsFixed1 )
cesBfgsFixed2 <- cesEst( "Y", xNames2, data = GermanIndustry,
   method = "BFGS", rho1 = 0.2155, rho = 1.1816,
   control = list( maxit = 5000 ) )
summary( cesBfgsFixed2 )
cesBfgsFixed3 <- cesEst( "Y", xNames3, data = GermanIndustry,
   method = "BFGS", rho1 = 1.3654, rho = 5.8327,
   control = list( maxit = 5000 ) )
summary( cesBfgsFixed3 )

## Levenberg-Marquardt, lambda, rho_1, and rho fixed
cesLmFixed1 <- cesEst( "Y", xNames1, data = GermanIndustry, 
   rho1 = 0.5300, rho = 0.1813,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed1 )
cesLmFixed2 <- cesEst( "Y", xNames2, data = GermanIndustry, 
   rho1 = 0.2155, rho = 1.1816,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed2 )
cesLmFixed3 <- cesEst( "Y", xNames3, data = GermanIndustry, 
   rho1 = 1.3654, rho = 5.8327,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed3 )

## Levenberg-Marquardt, lambda, rho_1, and rho fixed, multiplicative error term
cesLmFixed1Me <- cesEst( "Y", xNames1, data = GermanIndustry, 
   rho1 = 0.5300, rho = 0.1813, multErr = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed1Me )
summary( cesLmFixed1Me, rSquaredLog = FALSE )
cesLmFixed2Me <- cesEst( "Y", xNames2, data = GermanIndustry, 
   rho1 = 0.2155, rho = 1.1816, multErr = TRUE,
   control = nls.lm.control( maxiter = 1000, maxfev = 2000 ) )
summary( cesLmFixed2Me )
summary( cesLmFixed2Me, rSquaredLog = FALSE )
cesLmFixed3Me <- cesEst( "Y", xNames3, data = GermanIndustry, 
   rho1 = 1.3654, rho = 5.8327, multErr = TRUE,
   control = nls.lm.control( maxiter = 1024, maxfev = 2000 ) )
summary( cesLmFixed3Me )
summary( cesLmFixed3Me, rSquaredLog = FALSE )

## Newton-type, lambda, rho_1, and rho fixed
cesNewtonFixed1 <- cesEst( "Y", xNames1, data = GermanIndustry,
   method = "Newton", rho1 = 0.5300, rho = 0.1813, iterlim = 500 )
summary( cesNewtonFixed1 )
cesNewtonFixed2 <- cesEst( "Y", xNames2, data = GermanIndustry,
   method = "Newton", rho1 = 0.2155, rho = 1.1816, iterlim = 500 )
summary( cesNewtonFixed2 )
cesNewtonFixed3 <- cesEst( "Y", xNames3, data = GermanIndustry,
   method = "Newton", rho1 = 1.3654, rho = 5.8327, iterlim = 500 )
summary( cesNewtonFixed3 )

## PORT, lambda, rho_1, and rho fixed
cesPortFixed1 <- cesEst( "Y", xNames1, data = GermanIndustry,
   method = "PORT", rho1 = 0.5300, rho = 0.1813,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortFixed1 )
cesPortFixed2 <- cesEst( "Y", xNames2, data = GermanIndustry,
   method = "PORT", rho1 = 0.2155, rho = 1.1816,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortFixed2 )
cesPortFixed3 <- cesEst( "Y", xNames3, data = GermanIndustry,
   method = "PORT", rho1 = 1.3654, rho = 5.8327,
   control = list( eval.max = 1000, iter.max = 1000 ) )
summary( cesPortFixed3 )


# compare RSSs of models with lambda, rho_1, and rho fixed
print( matrix( c( cesNmFixed1$rss, cesBfgsFixed1$rss, cesLmFixed1$rss,
   cesNewtonFixed1$rss, cesPortFixed1$rss ), ncol = 1 ), digits = 16 )
cesFixed1 <- cesLmFixed1
print( matrix( c( cesNmFixed2$rss, cesBfgsFixed2$rss, cesLmFixed2$rss,
   cesNewtonFixed2$rss, cesPortFixed2$rss ), ncol = 1 ), digits = 16 )
cesFixed2 <- cesLmFixed2
print( matrix( c( cesNmFixed3$rss, cesBfgsFixed3$rss, cesLmFixed3$rss,
   cesNewtonFixed3$rss, cesPortFixed3$rss ), ncol = 1 ), digits = 16 )
cesFixed3 <- cesLmFixed3

# save the work space
save.image( "kemfert98_fixed.RData" )

## check if removing the technical progress worked as expected
Y2Calc <- cesCalc( xNames2, data = GermanIndustry, 
   coef = coef( cesFixed2 ), nested = TRUE )
all.equal( Y2Calc, fitted( cesFixed2 ) )
Y2TcCalc <- cesCalc( sub( "[123]$", "", xNames2 ), tName = "time", 
   data = GermanIndustry, 
   coef = c( coef( cesFixed2 )[1], lambda = 0.0069, coef( cesFixed2 )[-1] ), 
   nested = TRUE )
all.equal( Y2Calc, Y2TcCalc )


