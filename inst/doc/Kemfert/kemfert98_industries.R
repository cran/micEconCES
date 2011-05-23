# load the micEconCES package
library( "micEconCES" )

# load the data set
data( "GermanIndustry" )

# remove years with missing or incomplete data
GermanIndustry <- subset( GermanIndustry, year >= 1970 & year <= 1988, )

# remove years 1973 - 1975 because of economic disruptions (see Kemfert 1998)
GermanIndustry <- subset( GermanIndustry, year < 1973 | year > 1975, )

# add a time trend (starting with 0)
GermanIndustry$time <- GermanIndustry$year - 1970

# rhos for grid search
rhoVec <- c( seq( -1, 0.6, 0.2 ), seq( 0.9, 1.5, 0.3 ),
   seq( 2, 10, 1 ), 12, 15, 20, 30, 50, 100 )

# industries (abbreviations)
indAbbr <- c( "C", "S", "N", "I", "V",  "P", "F" )

# names of inputs
xNames <- list()
xNames[[ 1 ]] <-  c( "K", "E", "A" )
xNames[[ 2 ]] <-  c( "K", "A", "E" )
xNames[[ 3 ]] <-  c( "E", "A", "K" )

# list ("folder") for results
indRes <- list()

# names of estimation methods
metNames <- c( "LM", "PORT", "PORT_Grid", "PORT_Start" )

# table for parameter estimates
tabCoef <- array( NA, dim = c( 9, 7, length( metNames ) ),
   dimnames = list( 
      paste( rep( c( "alpha_", "beta_", "m_" ), 3 ), rep( 1:3, each = 3 ),
         " (", rep( c( "rho_1", "rho", "lambda" ), 3 ), ")", sep = "" ),
      indAbbr, metNames ) )

# table for technological change parameters
tabLambda <- array( NA, dim = c( 7, 3, length( metNames ) ),
   dimnames = list( indAbbr, c(1:3), metNames ) )

# table for R-squared values
tabR2 <- tabLambda

# table for RSS values
tabRss <- tabLambda

# table for economic consistency of LM results
tabConsist <- tabLambda[ , , 1, drop = TRUE ]


################# econometric estimation with cesEst ##########################
for( indNo in 1:length( indAbbr ) ) {

   # name of industry-specific output
   yIndName <- paste( indAbbr[ indNo ], "Y", sep = "_" )

   # sub-list ("subfolder") for all models of this industrie
   indRes[[ indNo ]] <- list()
   
   for( modNo in 1:3 ) {

      cat( "\n=======================================================\n" )
      cat( "Industry No. ", indNo, ", model No. ", modNo, "\n", sep = "" )
      cat( "=======================================================\n\n" )

      # names of industry-specific inputs
      xIndNames <- paste( indAbbr[ indNo ], xNames[[ modNo ]], sep = "_" )
      
      # sub-sub-list for all estimation results of this model/industrie
      indRes[[ indNo ]][[ modNo ]] <- list()
   
      ## Levenberg-Marquardt
      indRes[[ indNo ]][[ modNo ]]$lm <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry,
         control = nls.lm.control( maxiter = 1024, maxfev = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$lm ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$lm )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "LM" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "LM" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "LM" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "LM" ] <- tmpSum$rss
      tabConsist[ indNo, modNo ] <- tmpCoef[ "gamma" ] >= 0 &
         tmpCoef[ "delta_1" ] >= 0 & tmpCoef[ "delta_1" ] <= 1 &
         tmpCoef[ "delta" ] >= 0 & tmpCoef[ "delta" ] <= 1 &
         tmpCoef[ "rho_1" ] >= -1 & tmpCoef[ "rho" ] >= -1

      ## PORT
      indRes[[ indNo ]][[ modNo ]]$port <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry, method = "PORT",
         control = list( eval.max = 2000, iter.max = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$port ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$port )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "PORT" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "PORT" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "PORT" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "PORT" ] <- tmpSum$rss

      # PORT, grid search
      indRes[[ indNo ]][[ modNo ]]$portGrid <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry, method = "PORT",
         rho = rhoVec, rho1 = rhoVec,
         control = list( eval.max = 2000, iter.max = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$portGrid ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$portGrid )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "PORT_Grid" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "PORT_Grid" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "PORT_Grid" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "PORT_Grid" ] <- tmpSum$rss

      # PORT, grid search for starting values
      indRes[[ indNo ]][[ modNo ]]$portStart <- cesEst( yIndName, xIndNames, 
         tName = "time", data = GermanIndustry, method = "PORT",
         start = coef( indRes[[ indNo ]][[ modNo ]]$portGrid ),
         control = list( eval.max = 2000, iter.max = 2000 ) )
      print( tmpSum <- summary( indRes[[ indNo ]][[ modNo ]]$portStart ) )
      tmpCoef <- coef( indRes[[ indNo ]][[ modNo ]]$portStart )
      tabCoef[ ( 3 * modNo - 2 ):( 3 * modNo ), indNo, "PORT_Start" ] <- 
         tmpCoef[ c( "rho_1", "rho", "lambda" ) ]
      tabLambda[ indNo, modNo, "PORT_Start" ] <- tmpCoef[ "lambda" ]
      tabR2[ indNo, modNo, "PORT_Start" ] <- tmpSum$r.squared
      tabRss[ indNo, modNo, "PORT_Start" ] <- tmpSum$rss
   }
}

print( round( tabCoef, 3 ) )
print( round( do.call( "cbind", 
   lapply( 1:length( metNames ), function(x) tabR2[,,x] ) ), 3 ) )

print( round( 1e8*( tabR2[ , , "LM" ] - tabR2[ , , "PORT" ] ) ) )
print( round( 1e8*( tabR2[ , , "PORT_Grid" ] - tabR2[ , , "PORT" ] ) ) )
print( round( 1e8*( tabR2[ , , "PORT_Start" ] - tabR2[ , , "PORT" ] ) ) )
print( round( 1e8*( tabR2[ , , "PORT_Start" ] - tabR2[ , , "PORT_Grid" ] ) ) )

print( round( ( tabRss[ , , "PORT" ] - tabRss[ , , "LM" ] ) / 100 ) )
print( round( ( tabRss[ , , "PORT" ] - tabRss[ , , "PORT_Grid" ] ) / 100 ) )
print( round( ( tabRss[ , , "PORT" ] - tabRss[ , , "PORT_Start" ] ) ) )
print( round( ( tabRss[ , , "PORT_Grid" ] - tabRss[ , , "PORT_Start" ] ) ) )

print( tabConsist )

