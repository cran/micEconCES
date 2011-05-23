# load the micEconCES package
library( "micEconCES" )

# load estimation results
load( "kemfert98_nls.RData" )
load( "kemfert98_fixed.RData" )
load( "kemfert98_grid.RData" )

# tables for estimation results
result1 <- matrix( NA, nrow = 24, ncol = 9 )
rownames( result1 ) <- c( "Kemfert (1998)", "fixed",
   "Newton", "BFGS", "L-BFGS-B", "PORT", "LM", 
   "NM", "NM - PORT", "NM - LM", 
   "SANN", "SANN - PORT", "SANN - LM", 
   "DE", "DE - PORT", "DE - LM",
   "Newton grid", "Newton grid start", "BFGS grid", "BFGS grid start", 
   "PORT grid", "PORT grid start", "LM grid", "LM grid start" )
colnames( result1 ) <- c( 
   paste( "$\\", names( coef( cesLm1 ) ), "$", sep = "" ), 
   "c", "RSS", "$R^2$" )
result3 <- result2 <- result1

result1[ "Kemfert (1998)", "$\\lambda$" ] <- 0.0222
result2[ "Kemfert (1998)", "$\\lambda$" ] <- 0.0069
result3[ "Kemfert (1998)", "$\\lambda$" ] <- 0.00641

result1[ "Kemfert (1998)", "$\\rho_1$" ] <- 0.5300
result2[ "Kemfert (1998)", "$\\rho_1$" ] <- 0.2155
result3[ "Kemfert (1998)", "$\\rho_1$" ] <- 1.3654

result1[ "Kemfert (1998)", "$\\rho$" ] <- 0.1813
result2[ "Kemfert (1998)", "$\\rho$" ] <- 1.1816
result3[ "Kemfert (1998)", "$\\rho$" ] <- 5.8327

result1[ "Kemfert (1998)", "$R^2$" ] <- 0.9996
result2[ "Kemfert (1998)", "$R^2$" ] <- 0.786
result3[ "Kemfert (1998)", "$R^2$" ] <- 0.9986

result1[ "fixed", ] <- c( coef( cesFixed1 )[1], 0.0222, 
   coef( cesFixed1 )[-1], cesFixed1$convergence,
   cesFixed1$rss, summary( cesFixed1 )$r.squared )
result2[ "fixed", ] <- c( coef( cesFixed2 )[1], 0.0069,
   coef( cesFixed2 )[-1],cesFixed2$convergence,
   cesFixed2$rss, summary( cesFixed2 )$r.squared )
result3[ "fixed", ] <- c( coef( cesFixed3 )[1], 0.00641,
   coef( cesFixed3 )[-1],cesFixed3$convergence,
   cesFixed3$rss, summary( cesFixed3 )$r.squared )

makeRow <- function( model ) {
   if( is.null( model$multErr ) ) {
      model$multErr <- FALSE
   }
   result <- c( coef( model ), 
      ifelse( is.null( model$convergence ), NA, model$convergence ),
      model$rss, summary( model )$r.squared )
   return( result )
}

result1[ "Newton", ] <- makeRow( cesNewton1 )
result2[ "Newton", ] <- makeRow( cesNewton2 )
result3[ "Newton", ] <- makeRow( cesNewton3 )

result1[ "BFGS", ] <- makeRow( cesBfgs1 )
result2[ "BFGS", ] <- makeRow( cesBfgs2 )
result3[ "BFGS", ] <- makeRow( cesBfgs3 )

result1[ "L-BFGS-B", ] <- makeRow( cesBfgsCon1 )
result2[ "L-BFGS-B", ] <- makeRow( cesBfgsCon2 )
result3[ "L-BFGS-B", ] <- makeRow( cesBfgsCon3 )

result1[ "PORT", ] <- makeRow( cesPort1 )
result2[ "PORT", ] <- makeRow( cesPort2 )
result3[ "PORT", ] <- makeRow( cesPort3 )
   
result1[ "LM", ] <- makeRow( cesLm1 )
result2[ "LM", ] <- makeRow( cesLm2 )
result3[ "LM", ] <- makeRow( cesLm3 )

result1[ "NM", ] <- makeRow( cesNm1 )
result2[ "NM", ] <- makeRow( cesNm2 )
result3[ "NM", ] <- makeRow( cesNm3 )

result1[ "NM - LM", ] <- makeRow( cesNmLm1 )
result2[ "NM - LM", ] <- makeRow( cesNmLm2 )
result3[ "NM - LM", ] <- makeRow( cesNmLm3 )

result1[ "NM - PORT", ] <- makeRow( cesNmPort1 )
result2[ "NM - PORT", ] <- makeRow( cesNmPort2 )
result3[ "NM - PORT", ] <- makeRow( cesNmPort3 )

result1[ "SANN", ] <- makeRow( cesSann1 )
result2[ "SANN", ] <- makeRow( cesSann2 )
result3[ "SANN", ] <- makeRow( cesSann3 )

result1[ "SANN - LM", ] <- makeRow( cesSannLm1 )
result2[ "SANN - LM", ] <- makeRow( cesSannLm2 )
result3[ "SANN - LM", ] <- makeRow( cesSannLm3 )

result1[ "SANN - PORT", ] <- makeRow( cesSannPort1 )
result2[ "SANN - PORT", ] <- makeRow( cesSannPort2 )
result3[ "SANN - PORT", ] <- makeRow( cesSannPort3 )

result1[ "DE", ] <- makeRow( cesDe1 )
result2[ "DE", ] <- makeRow( cesDe2 )
result3[ "DE", ] <- makeRow( cesDe3 )

result1[ "DE - LM", ] <- makeRow( cesDeLm1 )
result2[ "DE - LM", ] <- makeRow( cesDeLm2 )
result3[ "DE - LM", ] <- makeRow( cesDeLm3 )

result1[ "DE - PORT", ] <- makeRow( cesDePort1 )
result2[ "DE - PORT", ] <- makeRow( cesDePort2 )
result3[ "DE - PORT", ] <- makeRow( cesDePort3 )

result1[ "Newton grid", ] <- makeRow( cesNewtonGridRho1 )
result2[ "Newton grid", ] <- makeRow( cesNewtonGridRho2 )
result3[ "Newton grid", ] <- makeRow( cesNewtonGridRho3 )

result1[ "Newton grid start", ] <- makeRow( cesNewtonGridStartRho1 )
result2[ "Newton grid start", ] <- makeRow( cesNewtonGridStartRho2 )
result3[ "Newton grid start", ] <- makeRow( cesNewtonGridStartRho3 )

result1[ "BFGS grid", ] <- makeRow( cesBfgsGridRho1 )
result2[ "BFGS grid", ] <- makeRow( cesBfgsGridRho2 )
result3[ "BFGS grid", ] <- makeRow( cesBfgsGridRho3 )

result1[ "BFGS grid start", ] <- makeRow( cesBfgsGridStartRho1 )
result2[ "BFGS grid start", ] <- makeRow( cesBfgsGridStartRho2 )
result3[ "BFGS grid start", ] <- makeRow( cesBfgsGridStartRho3 )

result1[ "PORT grid", ] <- makeRow( cesPortGridRho1 )
result2[ "PORT grid", ] <- makeRow( cesPortGridRho2 )
result3[ "PORT grid", ] <- makeRow( cesPortGridRho3 )

result1[ "PORT grid start", ] <- makeRow( cesPortGridStartRho1 )
result2[ "PORT grid start", ] <- makeRow( cesPortGridStartRho2 )
result3[ "PORT grid start", ] <- makeRow( cesPortGridStartRho3 )

result1[ "LM grid", ] <- makeRow( cesLmGridRho1 )
result2[ "LM grid", ] <- makeRow( cesLmGridRho2 )
result3[ "LM grid", ] <- makeRow( cesLmGridRho3 )

result1[ "LM grid start", ] <- makeRow( cesLmGridStartRho1 )
result2[ "LM grid start", ] <- makeRow( cesLmGridStartRho2 )
result3[ "LM grid start", ] <- makeRow( cesLmGridStartRho3 )


############## create LaTeX tables ####################
library( xtable )
colorRows <- function( result ) {
   rownames( result ) <- paste( 
      ifelse( !is.na( result[ , "$\\delta_1$" ] ) & (
         result[ , "$\\delta_1$" ] < 0 | result[ , "$\\delta_1$" ] >1 | 
         result[ , "$\\delta$" ] < 0 | result[ , "$\\delta$" ] > 1 ) | 
         result[ , "$\\rho_1$" ] < -1 | result[ , "$\\rho$" ] < -1, 
         "MarkThisRow ", "" ),
      rownames( result ), sep = "" )
   return( result )
}

printTable <- function( xTab, fileName ) {
   tempFile <- file()
   print( xTab, file = tempFile,
      floating = FALSE, sanitize.text.function = function(x){x} )
   latexLines <- readLines( tempFile )
   close( tempFile )
   for( i in grep( "MarkThisRow ", latexLines, value = FALSE ) ) {
      latexLines[ i ] <- sub( "MarkThisRow", "\\\\color{red}" ,latexLines[ i ] )
      latexLines[ i ] <- gsub( "&", "& \\\\color{red}" ,latexLines[ i ] )
   }
   writeLines( latexLines, fileName )
   invisible( latexLines )
}
   
result1 <- colorRows( result1 )
xTab1 <- xtable( result1, align = "lrrrrrrrrr", 
   digits = c( 0, rep( 4, 6 ), 0, 0, 4 ) )
printTable( xTab1, fileName = "../tables/kemfert1Coef.tex" )

result2 <- colorRows( result2 )
xTab2 <- xtable( result2, align = "lrrrrrrrrr", 
   digits = c( 0, rep( 4, 6 ), 0, 0, 4 ) )
printTable( xTab2, fileName = "../tables/kemfert2Coef.tex" )

result3 <- colorRows( result3 )
xTab3 <- xtable( result3, align = "lrrrrrrrrr", 
   digits = c( 0, rep( 4, 6 ), 0, 0, 4 ) )
printTable( xTab3, fileName = "../tables/kemfert3Coef.tex" )

