library(ape)
library(plotrix)




#PARAMETERS
#tree = phylo object
#scape.size = edge dimension of square landscape
#g.center = phylogenetic signal in species optimal values  phylogenetic attraction (range centers)   See corBlomberg in ape, 1=brownian,<1=rates of evol accelerate, >1=rates decelerate.
#g.range = phylogenetic signal in species niche widths  (range sizes). 1=brownian,<1=rates of evol accelerate, >1=rates decelerate.
#g.repulse = include phylogenetic repulsion
#wd.all = niche width denominator, larger values make larger on average range sizes
#signal.center = T/F simulate with phylosignal in range centers
#signal.range = T/F simulate with phylosignal in range sizes
#same.range = T/F make all range sizes equal
#repulse = T/F include phylogenetic repulsion
#center.scale = adjust the strength of phylogenetic attraction independent of signal (generally leave = 1 this is more for error checking)
#range.scale = adjust the strength of phylogenetic signal in niche width (range) (generally leave = 1 this is more for error checking)
#repulse.scale = adjust the strength of phylogenetic repulsion (generally leave = 1 this is more for error checking)
#site.stoch.scale = adjust the strength of random variation in "carrying capacity" across sites
#sd.center = sd in rnorm() for the range centers, increase to get more variation in center values across species
#sd.range =  sd in rnorm() for the range sizes, increase to get more variation in range sizes across gradients
#rho = vary the overall strength of signal in the phylogeny (Grafen branch adjustment) (probably more for error checking, so leave as NULL)
#th = probability threshold 10^-th above which species are considered present at a site (increase this value to obtain more species at sites in Y (p/a matrix) and thus increase average site SR) 
#env.type = gradient, patchy, or random

scape <- function(tree, 
				scape.size=10, 
				g.center=1, 
				g.range=1, 
				g.repulse=1, 
				wd.all=0.2*(scape.size+1)^2, 
				signal.center=TRUE, 
				signal.range=FALSE, 
				same.range=TRUE, 
				repulse=FALSE, 
				center.scale = 1, 
				range.scale = 1, 
				repulse.scale = 1, 
				site.stoch.scale = .5, 
				sd.center=1, 
				sd.range=1, 
				rho=NULL, 
				K=100,
				th=8,
				extinction = FALSE)
{
	
  #deal with the tree
    if (is(tree)[1] == "phylo")
    {
        if (is.null(tree$edge.length))
        {
            tree <- compute.brlen(tree, 1)    #Note this assigns arbitrary branch-lengths
        }
        V <- vcv.phylo(tree, corr = TRUE)     #Note this will mess with trees that are not ultrametric, such as those with argitrary branch-lengths
    } else {
        V <- tree
    }
    Vinit <- V       
    
    
  #initialize
    nspp <- dim(V)[1]
    bspp2 <- NULL
    Vcomp <- NULL

        Xscale <- 2          #scale the strength of the probability matrix X
        Mscale <- site.stoch.scale        #scale stochasticity in niche distributions 
        Vscale1 <- center.scale         #scale the strength of the optimal values on axis one
        Vscale2 <- center.scale         #scale the strength of the optimal values on axis two
    
    
 # Grafen's rho adjust strength of phylogenetic signal overall.
    if(!is.null(rho)){
     V <- 1-(1-Vinit)^rho
     V <- V/max(V)
    }
 
    #########################################################################################################################
#SIMULATION

    nsites<-scape.size #number of sites for the square landscape
    mx <- t(as.matrix((-(nsites)/2):(nsites/2)))  #env gradient
    m <- length(mx) #new number of sites (equal to nsites + 1)
     
    	      
############
#Establish range centers/niche optima 

    if(signal.center){
      g <- abs(g.center)
      V.a <- vcv(corBlomberg(g, tree),corr=T) #adjust phylogenetic signal for range centre
      iD <- t(chol(V.a))
    } else {
      V.a <- V
      iD <- diag(nspp) #diag = 1
    }
    if(signal.range){
      g <- abs(g.range)
      V.w <- vcv(corBlomberg(g, tree), corr=T) #adjust phylogenetic signal
      iD.w <- t(chol(V.w))
    } else {
      V.w <- V
      iD.w <- diag(nspp)
    }
    
  ##environmental/geographical gradient 1   

    e <- iD %*% rnorm(nspp, sd=sd.center)                                                               #assign optimal values as related to branch lengths, includes variation about mean value. Absolute values meaningless
    e <- Vscale1 * (e - mean(e)) / apply(e, 2, sd)                                           #z-scores and scaling of the optimal values based on environmental signal in phylogeny
    bspp1 <- e
    
#range size autocorrelation    
    if(!same.range){
   spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          wd <- range.scale*iD.w %*% rnorm(nspp,sd=sd.range)
          wd<-wd+(abs(min(wd)))
          wd<-wd/max(wd)

          dif<-sort(wd)[-1]-sort(wd)[-length(wd)]
          rat<-mean(dif/sort(wd)[-1])
          wd[wd==0]<-sort(wd)[2]-sort(wd)[2]*rat                                #Assign the zero with the mean ratio of nearest neighbor distances over the larger item
          wd<-wd.all*wd
          X <- exp(-((spmx - mxsp)^2)/t(matrix(wd,nspp,m))) #Niche distributions
        } else {
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          X <- exp(-((spmx - mxsp)^2)/wd.all) #Niche distributions
    }       
    X <- Xscale * X                                                                       #Scales this initial species x site probability matrix 
    #Xsmooth <- X                                                                         #Distributions without random variation
    X1 <- diag(1 - Mscale * runif(m)) %*% X                                               #Scale and include random variation into the niche distributions

    ##environmental/geographical gradient 2
        e <- iD %*% rnorm(nspp,sd=sd.center)
        e <- Vscale2 * (e - mean(e))/apply(e, 2, sd)
        bspp2 <- e


  if(!same.range){
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          wd <- range.scale*iD.w %*% rnorm(nspp,sd=sd.range)
          wd<-wd+(abs(min(wd)))
          wd<-wd/max(wd)
              #Assign the zero to the nonzero minimum
          dif<-sort(wd)[-1]-sort(wd)[-length(wd)]
          rat<-mean(dif/sort(wd)[-1])
          wd[wd==0]<-sort(wd)[2]-sort(wd)[2]*rat                                #Assign the zero with the mean rato of nearist neighbor distances over the larger item
          wd<-wd.all*wd
          X <- exp(-((spmx - mxsp)^2)/t(matrix(wd,nspp,m))) #Niche distributions     
     } else {
          spmx <- t((array(1, c(nspp, 1))) %*% mx)
          mxsp <- max(mx)*((array(1, c(length(mx), 1))) %*% t(e))
          X <- exp(-((spmx - mxsp)^2)/wd.all) #Niche distributions
     }
     X <- Xscale * X
     X2 <- diag(1 - Mscale * runif(m)) %*% X
     
     

     ##################################        
     #REPULSION
     X.repulse <- NULL
      if (repulse) {
        compscale <- repulse.scale
        b0scale <- 0
        g<-abs(g.repulse)
        V.r<-vcv(corBlomberg(g, tree),corr=T)          #adjust phylogenetic signal 
      #calculate the repulsion matrix
        Vcomp <- solve(V.r, diag(nspp))
        Vcomp <- Vcomp/max(Vcomp)
        Vcomp <- compscale * Vcomp
        iDcomp <- t(chol(Vcomp))
        colnames(Vcomp) <- rownames(Vcomp)
        bcomp <- NULL
        for (i in 1:m) {
          bcomp <- cbind(bcomp, iDcomp %*% rnorm(nspp))
        }
        bcomp0 <- 0
        Xcomp <- exp(bcomp0 + bcomp)/(1 + exp(bcomp0 + bcomp))
        #X <- X * t(Xcomp)
        X1<-X1 * t(Xcomp)
        X2<-X2 * t(Xcomp)
        X.repulse<-t(Xcomp)
      }
      
#JOINT PROBABILITY MATRIX
  X.<-NULL
  spp.Xs<-array(NA,dim=c(m,m,nspp))
  for(i in 1:nspp){
    sppX<-matrix((X1[,i]) %*% t(X2[,i]))
    spp.Xs[,,i]<-sppX
    X.<-cbind(X.,matrix(sppX))
   }
   colnames(X.)<-colnames(X2)
  ######################
  #PA matrix
  m.<- dim(X.)[1]
  Y <- matrix(0, ncol = nspp, nrow = m.)
  Y[10^-th < X.] <- 1     #could also use a hard threshold
  colnames(Y) <- colnames(X.)
  index<-NULL
  
   
  ######################
  #PA matrix
  m. <- dim(X.)[1]
  Y <- matrix(0, ncol = nspp, nrow = m.)
  Y[10^-th < X.] <- 1
  if(extinction==FALSE){
  for(i in  which(colSums(Y)==0)){
  	Y[sample((which(X.[,i]==max(X.[,i]))),1),i] <- 1
  		}
  	}
  colnames(Y) <- tree$tip.label #colnames(X.)
  index <- NULL
  index <- cbind(matrix(sapply(1:m, rep, times=m)), matrix(rep(1:m, times=m)))
  colnames(index) <- c("X1", "X2")
  
  
##WEIGHT by K for abundance matrix- new range of values for each site bounded by 0 (absent), each site sums to K
	Yab <- X.
	Yab[10^-th > X.] <- 0
  
	Yab <- t(apply(Yab, 1, 
  		function(x){
  			if(sum(x)>0){
  				ceiling(x*K/sum(x))
  				}else{
  				x
  				}})) #scale by carrying capacity K			
	Yab[Yab>0] <- floor(sapply(Yab[Yab>0], function(x){runif(1,min=(x-1),max=(x+5))}))
	for(i in which(colSums(Yab)==0)){
  		Yab[(which(Y[,i]==1)),i] <- 1
  		}
	colnames(Yab) <- tree$tip.label
 

##CREATE full environmental matrix (add mx1+mx2)
	env <- matrix(mx, nrow=length(mx), ncol=length(mx), byrow=TRUE)
for(i in 1:nrow(env)){
	env[i,] <- sapply(env[i,],function(x){x+mx[i]})}
 
 
########### OUTPUT  
    	return(list(Y = Y, 
    			    Yab = Yab,
    				index = index, 
    				gradient1 = mx,
    				gradient2 = mx,
                	X.joint = X.,
                	X1 = X1, 
                	X2 = X2, 
                	nichewd = wd.all, 
                	K = K, 
					environ = env,
                	sppXs = spp.Xs, 
                	V.phylo = Vinit, 
                	V.phylo.rho = V, 
                	V.center = V.a, 
                	bsp1 = bspp1, 
                	bspp2 = bspp2
                	))
  
}  #function end



#VALUES
#Y = presence/absence matrix
#Yab = abundance matrix
#index = spatial coordinates for X and Y (stacked columns)
#X.joint = full probabilities of species at sites, used to construct Y 
#X1 = probabilities of species along gradient 1
#X2 = probabilities of species along gradient 2
#gradient1 & gradient2 = env. gradient values (equivalent to lat/long)
#nichewd = average niche width of the assemblage
#K = carrying capacity of each cell
#environ = matrix depicting environmental values over the landscape
#linear, quadratic, stochastic, random, uniform = cost functions based on env

#sppXs = full probabilities of each species as an array arranged in a scape.size X scape.size matrix  
#V.phylo = initial phylogenetic covariance matrix from tree, output of vcv.phylo(tree, corr=T)
#V.phylo.rho = phylogenetic covariance matrix from tree scaled by Grafen if rho is provided, otherwise just an output of vcv.phylo(tree, corr=T)
#V.center = scaled (by g.center) phylo covariance matrix used in the simulations
#bspp1 = species optima for gradient 1
#bspp2 = species optima for gradient 2                               
################END###########################





###EXAMPLES##############################################################
##Load Libraries
#require(ape)


##generate trees with a variety of shapes

#tree<-rtree(128)
#tree <- stree(128, "balanced")

#scape1 <- scape(tree, 
#				scape.size=30, 
#				g.center=1, 
#				g.range=1, 
#				g.repulse=1, 
#				wd.all=0.2*(scape.size+1)^2, 
#				signal.center=TRUE, 
#				signal.range=TRUE, 
#				same.range=TRUE, 
#				repulse=TRUE, 
#				center.scale = 1, 
#				range.scale = 1, 
#				repulse.scale = 1, 
#				site.stoch.scale = .5, 
#				sd.center=1, 
#				sd.range=1, 
#				rho=NULL, 
#				th=8,
#				K=100,
#				env.type="gradient",
#				extinction = FALSE)
#


##Plotting distributions and landscape patterns
plotscape <- function(scape1){
			abundmax<-scape1$K
			PA_mat<-as.matrix(scape1$Y)
			abund_mat<-scape1$Yab
			site.size=nrow(PA_mat)
			species<-ncol(PA_mat)
			mx<-scape1$gradient
			env<-scape1$environ


par(mfrow=c(2,2), oma=c(0,0,2,0))

#plot env gradient
heatcol<-(colorRampPalette(c("lightgoldenrodyellow","yellow","red","darkred")))
image(matrix(env,sqrt(site.size),sqrt(site.size),byrow=TRUE),col=heatcol(length(env)),xaxt="n",yaxt="n",main="Env gradient")
#plot SR
image(matrix(rowSums(PA_mat),sqrt(site.size),sqrt(site.size),byrow=TRUE),col=heatcol(species),xaxt="n",yaxt="n",main="Species Richness")

#species x area
hist(colSums(PA_mat),ylab="Number of species",xlab="Number of sites",main="Species Area Relationship",col="lightgrey")
hist(colSums(abund_mat),ylab="Number of species",xlab="Number of individuals",main="Species Abundance Relationship",col="lightgrey")
#mtext("Env random, phy.signal=0.2, 32 species", outer=TRUE, side=3, cex=1.25)
}

