invisible( '
Code to analyse output of mltest.R and generate plots 
')

parmfns <- list.files( path= simdir, pattern = 'parms.csv', full.name=T)
as.data.frame( t( sapply( parmfns, read.csv ) ) )  -> parmtab
fn2pid <- function(fn){
	x <- regmatches( fn, regexpr( '/[0-9]+', fn ))
	regmatches( x, regexpr( '[0-9]+', x ))
}
rownames( parmtab ) <- sapply( rownames(parmtab ) , fn2pid )

fitdata <- readRDS('mltest0.rds')
fits.tl <- fitdata[[2]]
fits <- fitdata[[1]]
names(fits) <- sapply( names(fits), fn2pid )
names(fits.tl ) <- sapply( names( fits.tl) , fn2pid )

parmsfits <- as.data.frame(  t( exp( sapply( fits, '[[', 'par' ) ) ) )
parmsfits$assortprob <- parmsfits$assortprob / ( 1 + parmsfits$assortprob ) 
parmsfits.tl <- as.data.frame( t( exp( sapply( fits.tl, '[[', 'par' ) ) ) )
parmsfits.tl$assortprob <- parmsfits.tl$assortprob / ( 1 + parmsfits.tl$assortprob ) 

#~ 		plot( parmsfits$wacute, parmtab$wacute[ match( rownames(parmsfits), rownames(parmtab ) ) ] )
#~ 		abline( a = 0,  b = 1 )
#~ 		plot( parmsfits.tl$wacute, parmtab$wacute[ match( rownames(parmsfits.tl), rownames(parmtab ) ) ] )
#~ 		abline( a = 0,  b = 1 )
cat( 'saving wacute.png\n' )
png('wacute.png') 
	matplot(parmtab$wacute[ match( rownames(parmsfits), rownames(parmtab ) ) ] ,  cbind(  parmsfits$wacute, parmsfits.tl$wacute ) , xlab = 'True', ylab = 'Estimated', main = 'wacute' ) 
	abline( a = 0, b = 1)
dev.off() 

cat( 'saving wrisk2.png\n' )
png('wrisk2.png')
	matplot(parmtab$wrisk2[ match( rownames(parmsfits), rownames(parmtab ) ) ] ,  cbind(  parmsfits$wrisk2, parmsfits.tl$wrisk2 ) , ylim = c(0, 20), xlab = 'True', ylab = 'Estimated', main = 'wrisk2') 
	abline( a = 0, b = 1)
dev.off()

cat( 'saving assortprob.png\n')
png( 'assortprob.png')
par(mfrow = c(2,1))
hist( parmsfits$assortprob, main = 'Likelihood 1' );abline ( v= .8, col = 'red')
hist( parmsfits.tl$assortprob, main = 'Likelihood 2' );abline ( v= .8, col = 'red')
dev.off() 

# bias, accuracy, precision 
cat('bias wacute:\n')
cat( 'Likelihood 1: \n')
print ( summary( unlist( parmtab$wacute[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits$wacute ) )
cat( 'Likelihood 2: \n')
print( 
summary( unlist( parmtab$wacute[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits.tl$wacute )
)

cat('\n\n-------')
cat( 'rmse wacute:\n')
cat( 'Likelihood 1: \n')
print( sqrt( mean( ( unlist( parmtab$wacute[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits$wacute )^2 ) ) )
cat( 'Likelihood 2: \n')
print( sqrt( mean( ( unlist( parmtab$wacute[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits.tl$wacute )^2 ) ) )

cat('\n\n-------')
cat('\n\n-------')
cat('bias wrisk2:\n')
cat( 'Likelihood 1: \n')
print( summary( unlist( parmtab$wrisk2[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits$wrisk2 ) )
cat( 'Likelihood 2: \n')
print( summary( unlist( parmtab$wrisk2[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits.tl$wrisk2  ) )
cat('\n\n-------')
cat( 'rmse wacute:\n')
cat( 'Likelihood 1: \n')
print( sqrt( mean( ( unlist( parmtab$wrisk2[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits$wrisk2 )^2 ) ) )
cat( 'Likelihood 2: \n')
print( sqrt( mean( ( unlist( parmtab$wrisk2[ match( rownames(parmsfits), rownames(parmtab ) ) ] ) -    parmsfits.tl$wrisk2 )^2 ) ) )
