require(ape)

epi2tree <- function( transmission_table, removal_table )
{
# vectors donor, recip, t same length
# vectors pids, tremoval same length; all pids unique, matches donor and recip
# should make this two data frames
# all recip has exactly one donor
# earliest donor should be src
# ~ 	donor, recip, t
print(date())
#~ 	attach( transmission_table)
	donor = transmission_table$donor
	recip = transmission_table$recip
	t = transmission_table$t
	
	donor[is.na(donor)] <- 'src'
	
# ~ 	pids, tremoval, tiplab
#~ 	attach( removal_table )
	pids = removal_table$pids
	tremoval = removal_table$tremoval
	tiplab = removal_table$tiplab

	#TODO should validate input

	N <- length(pids)
	Nnode <- N - 1 + 1 #include source +1 #N + Nnode
	edge <- matrix( NA, nrow = length(t) + length(pids) , ncol = 2)
	edge.length <- rep(NA, length(t) + length(pids))
	recip2donor <- setNames( recip, donor )

	i <- order(t)
	donor <- donor[i]
	recip <- recip[i]
	t <- t[i]
print(c(date(), 'sorted inputs'))
	recip2donorNode <- setNames( rep(NA, length(pids)), pids)
	recip2donorNode[ recip2donor[pids]=='src' ] <- 'src'

	names(tremoval) <- pids

	# add transm edges
	donor2trCounter <- setNames(rep(0, length(pids)+1), c('src', pids) )
	donor2ntr <- sapply( pids, function(pid) sum( donor== pid ))
print(c(date(), 'donor2ntr'))
	node2time <- list()
	node2time[['src']] <- t[1]-1
	k <- 1
print( c(date(), 'frontmatter' ))
	for ( i in 1:length(t)){
		u <- donor[i]
		v <- recip[i]
		donor2trCounter[u] <- donor2trCounter[u] + 1
		trcount <- donor2trCounter[u]
if (is.na(trcount)) browser()
		if (u == 'src'){
			donornode <- u
		} else{
			donornode <- paste(u, trcount, sep='_')
		}
		recip2donorNode[v] <- donornode
		if (u != 'src'){
			edge[k,2] <- donornode
			node2time[[ donornode ]] <- t[i]
			if (trcount==1){
				edge[k, 1] <- 	recip2donorNode[u]
			} else{
				edge[k,1] <- paste(u, trcount-1, sep='_')
			}
			k <- k + 1
		}
		if (i %% 500 == 0) print(c( date(), i))
	}

	# add terminals
	names(tiplab) <- pids
	for ( i in 1:length( pids )){
		pid <- pids[i]
		tl <- tiplab[pid]
		node2time[[tl]] <- tremoval[i]
		if ( donor2trCounter[pid] == 0 ){
			lastnode <- recip2donorNode[pid]
		} else{
			trcount <- donor2trCounter[pid]
			lastnode <- paste(pid, trcount, sep='_')
		}
		edge[k,1] <- lastnode
		edge[k,2] <- tl
		k <- k + 1
	}
	
	internalNodes <- setdiff( unique( as.vector( edge ) ), tiplab)
	
	i_internalNodes <- setNames( N + 1:length(internalNodes), internalNodes )
	i_tips <- setNames( 1:N, tiplab)
	nodemap <- c( i_internalNodes, i_tips)
	
	edge <- edge[!is.na(edge[,1]), ]
	edge2 <- matrix(NA, nrow =nrow(edge), ncol =ncol(edge))
	edge2[,1] <- nodemap[ edge[,1] ]
	edge2[,2] <- nodemap[ edge[,2] ]
	
	edge.length <- rep(NA, nrow(edge2))
	edge.length <- unname( unlist( node2time[edge[,2]] ) - unlist(node2time[ edge[,1] ] ) )
	
	o <- list( edge = edge2, tip.label = tiplab, edge.length = edge.length, Nnode = Nnode )
	class(o) <- 'phylo'
	tre <- multi2di( read.tree(text= write.tree( o )) )
}


if (F)
{ # test on sir0
	transmission_table <- read.csv( 'sir0-transm.csv', header=F)
	colnames(transmission_table ) <-  c('donor', 'recip', 't')
	for (x in colnames(transmission_table)) transmission_table[[x]] <- as.vector( transmission_table[[x]] )
	transmission_table$donor <- as.character( transmission_table$donor )
	transmission_table$recip <- as.character( transmission_table$recip )

	removal_table <- read.csv( 'sir0-rem.csv', header=F)
	colnames(removal_table ) <-  c('pids', 'tremoval', 'tiplab')
	for (x in colnames(removal_table)) removal_table[[x]] <- as.vector( removal_table[[x]] )
	removal_table$pids <- as.character( removal_table$pids ) 

	#~ tre <- epi2tree( transmission_table, removal_table )
}


