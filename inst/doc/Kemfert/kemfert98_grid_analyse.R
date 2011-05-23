# load the micEconCES package
library( "micEconCES" )

# load results from grid searches
load( "kemfert98_grid.RData" )

hist( c( cesLmGridRho1$rssArray ),30 )

range( cesLmGridRho1$rssArray )

par( mar = c( 0, 1.5, 3, 0 ) )
plot( cesLmGridRho1 )

plot( cesLmGridRho1, zlim = c( -100000, -3400 ) )

cesLmGridRho1a <- cesLmGridRho1
cesLmGridRho1a$rssArray[ cesLmGridRho1a$rssArray > 5000 ] <- NA
plot( cesLmGridRho1a )

cesLmGridRho1b <- cesLmGridRho1
cesLmGridRho1b$rssArray[ cesLmGridRho1b$rssArray < 6000 ] <- NA
cesLmGridRho1b$rssArray[ cesLmGridRho1b$rssArray > 25000 ] <- NA
plot( cesLmGridRho1b )

cesLmGridRho1c <- cesLmGridRho1
cesLmGridRho1c$rssArray[ cesLmGridRho1c$rssArray < 60000 ] <- NA
plot( cesLmGridRho1c )

a <- sapply( 1:3721, function(x) {
   c( coef( cesLmGridRho1$allRhoFull[[x]] ), 
      rss = cesLmGridRho1$allRhoFull[[x]]$rss ) } )

hist( a[ "rho", a[ "rss", ] < 3700 ] )

plot( a[ "rho", a[ "rho_1", ] == 10 ],  a[ "rss", a[ "rho_1", ] == 10 ] )

a[ , a[ "rho_1", ] == 10 ]

a[ c( "rss", "rho", "delta_1" ), a[ "rho_1", ] == 10 ]

a[ , abs( a[ "delta_1", ] )  > 0.01 & a[ "rss", ] < 6000 ]

range( a[ "delta_1", a[ "rss", ] < 5000 ] )


plot( cesPortGridRho1 )

b <- sapply( 1:3721, function(x) {
   c( coef( cesPortGridRho1$allRhoFull[[x]] ), 
      rss = cesPortGridRho1$allRhoFull[[x]]$rss ) } )

b[ c( "rss", "rho", "delta_1" ), a[ "rho_1", ] == 10 ]

