require(phydynR.r0.1.2)
require(parallel)


simdir0 <- commandArgs(trailingOnly=TRUE)[1]
simdir <- paste0(simdir0, '/' )
#cargs <- '77525.nwk'
NCPU <- 40
trefns <- list.files( path=simdir, pattern = '[0-9]+.nwk' , full.names=TRUE)

outfnstem <- 'mltest0'

## parms
parms0 <- list( 
 mu = 1/40
 , gamma0 = 1
 , gamma1 = 1/9
 , beta = 1/4
 , wrisk2 = 5
 , wacute = 5
 , init_acute = 1
 , S0 = 5e3
 , prisk2 = .30
 , assortprob = 0.80 
 , agerate = 1/5.
)
demes <- c('acute', 'chron', 'acute2', 'chron2')
m <- 4

## the model 
#~ f = (beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2))
#~ W = (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)
B <- matrix( '0.', nrow=m, ncol=m)
rownames(B) = colnames(B) <- demes
B['acute', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * 
  (wacute * acute / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
  (assortprob + (1-assortprob) * (1-prisk2))'
B['acute', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (wacute * acute / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
  ( (1-assortprob) * (prisk2))'
B['chron', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) *
   (assortprob + (1-assortprob) * (1-prisk2))'
B['chron', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
   ( (1-assortprob) * (prisk2))'

B['acute2', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (wacute *wrisk2 * acute2 / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
  ( (1-assortprob) * (1-prisk2))'
B['acute2', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (wacute*wrisk2 * acute2 / (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) * 
  (assortprob + (1-assortprob) * prisk2)'
B['chron2', 'acute'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron2 *wrisk2/ (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) *
  ( (1-assortprob) * (1-prisk2))'
B['chron2', 'acute2'] <- '(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) * (chron2 *wrisk2/ (wacute * acute + wacute * wrisk2 *acute2 + chron + wrisk2*chron2)	) *
  (assortprob + (1-assortprob) * prisk2)'

M <- matrix( '0.', nrow=m, ncol=m)
rownames(M) = colnames(M) <- demes
M['acute', 'chron'] <- 'gamma0 * acute'
M['acute2', 'chron2'] <- 'gamma0 * acute2'
M['acute2','acute'] <- 'agerate * acute2'
M['chron2', 'chron'] <- 'agerate * chron2'

D <- setNames(rep('0.', m), demes)
D['acute'] <- 'mu * acute'
D['chron'] <- 'mu * chron + gamma1*chron'
D['acute2'] <- 'mu * acute2'
D['chron2'] <- 'mu * chron2 + gamma1 * chron2'

NDD <- c( S = '-(beta * (acute + acute2 + chron + chron2) * S / (S+acute+chron+acute2+chron2)) + mu * S0 - mu *S' )

dm <- build.demographic.process ( B, migrations=M, death=D, nonDeme=NDD, parameter=names(parms0) , rcpp=TRUE)

## LOAD TREE(S)

treefn <- trefns[1] 
#paste0(simdir, cargs[1])
tr <- read.tree( treefn)
n <- length(tr$tip.label)
ssts <- matrix(0, nrow = n, ncol = m)
colnames(ssts) <- demes
rownames(ssts) <- tr$tip.label
annots <- sapply( strsplit( tr$tip.label, '_' ), '[[', 2 )
ssts[ cbind( tr$tip.label, annots) ] <- 1
sts <- setNames( node.depth.edgelength( tr )[1:n], tr$tip.label)
bdt <- DatedTree( tr,sts,  sampleStates = ssts )

bdts <- mclapply( trefns , function(treefn){
	tr <- read.tree( treefn)
	n <- length(tr$tip.label)
	ssts <- matrix(0, nrow = n, ncol = m)
	colnames(ssts) <- demes
	rownames(ssts) <- tr$tip.label
	annots <- sapply( strsplit( tr$tip.label, '_' ), '[[', 2 )
	ssts[ cbind( tr$tip.label, annots) ] <- 1
	sts <- setNames( node.depth.edgelength( tr )[1:n], tr$tip.label)
	DatedTree( tr,sts,  sampleStates = ssts , tol = 1.)
}, mc.cores = NCPU )



## ofuns, estimate beta, wrisk2, wacute, I0
x0 <- c( acute = 4, chron = 0.001, acute2= 0.001, chron2=.001, S=parms0$S0)


#PL; 
of.pl <- function(theta, bdt = bdt ){
	parms1 <- parms0
	parms1$wrisk2 = unname( exp( theta['wrisk2'] ))
	parms1$wacute = unname( exp( theta['wacute'] ))
	parms1$beta = unname( exp( theta['beta'] ))
	ap <- theta['assortprob'] ; parms1$assortprob = unname( exp(ap) / (1 + exp(ap))  )
	x1 <- x0 
	x1['acute'] <- unname( exp( theta['init_acute'] ))
	print( unlist( parms1 ))
	-colik( bdt, theta=parms1, dm, x0 = x1, t0=0, res =200, forgiveAgtY = 1, AgtY_penalty = 1, step_size_res=10)
}
#~ system.time( print( of.pl( theta0, bdts[[1]] ))) 

#new PL; 
of.pl1 <- function(theta, bdt = bdt ){
	parms1 <- parms0
	parms1$wrisk2 = unname( exp( theta['wrisk2'] ))
	parms1$wacute = unname( exp( theta['wacute'] ))
	parms1$beta = unname( exp( theta['beta'] ))
	ap <- theta['assortprob'] ; parms1$assortprob = unname( exp(ap) / (1 + exp(ap))  )
	x1 <- x0 
	x1['acute'] <- unname( exp( theta['init_acute'] ))
	print( unlist( parms1 ))
	-colik( bdt, parms1, dm, x0 = x1, t0=0, res =200, forgiveAgtY = 1, AgtY_penalty = 1, step_size_res=10, test_likelihood=TRUE)
}
#~ system.time( print( of.pl1( theta0, bdts[[1]] ))) 

theta0 <- log( c( wrisk2 = 5, wacute = 5, beta = .25, init_acute=3, assortprob = .75/(1-.75) ))
#~ mll <- of.pl( theta0 )
#~ system.time( {o <- dm( parms0, x0, 0, 250)} )

sim.dm <- function( theta )
{
	parms1 <- parms0
	parms1$wrisk2 = unname( exp( theta['wrisk2'] ))
	parms1$wacute = unname( exp( theta['wacute'] ))
	parms1$beta = unname( exp( theta['beta'] ))
	ap <- theta['assortprob'] ; parms1$assortprob = unname( exp(ap) / (1 + exp(ap))  )
	x1 <- x0 
	x1['acute'] <- unname( exp( theta['init_acute'] ))
	print( unlist( parms1 ))
	dm(  parms1,  x0 = x1, t0=0, t1 = 100, res =50)
}
#~ x <- sim.dm( theta0 )

if (TRUE)
{
	fits1 <- mclapply( bdts, function(.bdt){
		optim( par = theta0, fn = of.pl1,  control=list(abstol=1e-4,trace=6), hessian=FALSE, bdt = .bdt )
	}, mc.cores = NCPU )
	ofn1 <- paste0( outfnstem, '.fitsTL.', simdir0, '.rds')
	saveRDS(fits1, file=ofn1) 
	
	fits0 <- mclapply( bdts, function(.bdt){
		optim( par = theta0, fn = of.pl,  control=list(abstol=1e-4,trace=6), hessian=FALSE, bdt = .bdt )
	}, mc.cores = NCPU )
	ofn0 <- paste0( outfnstem , '.fits.', simdir0, '.rds')
	saveRDS(fits0, file=ofn0) 
	
	names(fits0) <- trefns
	names(fits1) <- trefns 
	saveRDS(list( fits0, fits1 ) , file=paste0(outfnstem, '.rds') )
	
}

