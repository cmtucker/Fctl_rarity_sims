##Sparseness and distinctiveness calculation functions 
###Measure traits and biogeographical measures

#The output from the scape function is a list including scape$Yab, an abundance site x species matrix; and bsp1, bspp2 - these are 2 sets of trait values associated with the 2 environmental gradients. 


##########Regional Scale Measures of rarity###############

#calculate distance matrix for traits of species in entire regional pool

regional_measures <- function(scapeout){#I've set it up as a single function for all regional measures.
	
regionalDa <- as.matrix(dist(scapeout$bsp1, upper=TRUE, diag=NULL))
regionalDb <- as.matrix(dist(scapeout$bspp2, upper=TRUE, diag=NULL))
regionalD <- (regionalDa + regionalDb)/2
diag(regionalD) <- NA

#rescale trait-distance matrix between 0-1
regionalD <- (regionalD)/max(regionalD, na.rm=TRUE)

#Calculate Regional Uniqueness (U)
U <- apply(regionalD, 1, function(x){min(x, na.rm=TRUE)})
#uniqueness is the closest distance in trait space between species i and all other species in the species pool

#Regional Sparseness SregionalS)
	#species occupancy
	occupancy <- colSums(scapeout$Y)
	max_occ <- max(occupancy)
	
	#scaled between 0-1
	regionalS <- 1-(occupancy/max_occ)

#Calculate the comprehensive rarity measure (regionalR)
regionalR <- (regionalS + U)/2

#prepare the three measures for output
out <- cbind(regionalR, regionalS, U)

#Calculate quantiles for regionalS and U (use 95% and 5%)
qS <- quantile(regionalS, probs=c(0.05, 0.95))
qU <- quantile(U, probs=c(0.05, 0.95)) 

class <- matrix(NA, nrow=nrow(out), ncol=1)

##Classify each species into the 9 possible categories of rarity. I'm just using if and else statements, probably there is a nicer method!

for(i in 1:nrow(class)){
	
if(out[i,"U"]>=qU[2]){ #very distinct
			class[i] <-ifelse(out[i,"regionalS"]<=qS[1], 1, ifelse(out[i,"regionalS"]>=qS[2], 3, 2))
	}else{
		if(out[i,"U"]<=qU[1]){ #not distinct
			class[i] <- ifelse(out[i,"regionalS"]<=qS[1], 7, ifelse(out[i,"regionalS"]>=qS[2], 9, 8))
	}else{
		if(out[i,"U"]>=qU[1]&out[i,"U"]<=qU[2]){ #med distinct
			class[i] <-ifelse(out[i,"regionalS"]<=qS[1], 4, ifelse(out[i,"regionalS"]>=qS[2], 6, 5))
		}
	}}
}
out <- cbind(out, as.vector(class))				
return(out)
}

#some plotting numbers of categories
regional_hist <- function(out){
	class<-as.vector(out[,4])
	class <- class[!is.na(class)]
	class <- as.factor(class)
	return(c(sum(class==1), sum(class==2), sum(class==3), sum(class==4), sum(class==5), sum(class==6), sum(class==7), sum(class==8), sum(class==9)))
	}
	

####################Local scale##########################

#Local scale measures are calculated for each community (x); the simulation produces multiple communities forming a larger landscape. 

#Local Distinctiveness (local_distinct)
local_measures <-function(scapeout){#function to calculate all local measures

#prepare matrix to save all output	
local_distinct <- matrix(NA, nrow=nrow(scapeout$Yab), ncol=ncol(scapeout$Yab))
#make sure the naming is consistent
colnames(local_distinct) <- colnames(scapeout$Yab)
rownames(scapeout$bsp1) <- colnames(scapeout$Yab)
rownames(scapeout$bspp2) <- colnames(scapeout$Yab)

#For each community in the landscape:
for(i in 1:nrow(scapeout$Yab)){
	
		#subset to look only at species present in the comm
		a <- scapeout$Yab[i,which(scapeout$Yab[i,]>0)]
		
		if(length(a)>0){
		#subset trait data to match community data
		b <- scapeout$bsp1[match(names(a), rownames(scapeout$bsp1))]
		b2 <- scapeout$bspp2[match(names(a), rownames(scapeout$bspp2))]
		names(b) <- names(a)
		names(b2) <- names(a)

		#calculate community trait distances
		c1 <- as.matrix(dist(b, upper=TRUE, diag=NULL))
		c2 <- as.matrix(dist(b2, upper=TRUE, diag=NULL))
		c <- as.matrix(dist(b, upper=TRUE, diag=NULL))
		c <- (c1 + c2) / 2
		
		#scale by max to be between 0-1
		c <- c/max(c, na.rm=TRUE)
		
		#calculate distinctivness in reference to the abundances of the other species in the community in d
		#this is the abundance weighted version
		d <- matrix(NA, ncol=length(a), nrow=1)
		for (g in 1:length(d)){
			d[g] <- sum(a[-g]*c[-g, g])/sum(a[-g])
		}
		colnames(d) <- colnames(c)#names matching again
		
		#record in the local_distinct matrix, associate distinctiveness values with the correct species
		local_distinct[i, c(which(match(colnames(local_distinct), colnames(d))!="NA"))] <- d		
		#where species aren't present in a community, their distinctiveness==NA
		
		}}
	#warning message re max, if there are empty patches


#Local Sparseness (local_sparse)
local_sparse <- matrix(NA, nrow=nrow(scapeout$Yab), ncol=ncol(scapeout$Yab)) #matrix for output
colnames(local_sparse) <- colnames(scapeout$Yab)

#For each community (j)
for(j in 1:nrow(scapeout$Yab)){
	
		#subset to species with abundances >0
		a <- scapeout$Yab[j,which(scapeout$Yab[j,]>0)]
		
		if(length(a)>0){
			
		b <- a/sum(a) #relative abundance
		#sparseness calculation
		c <- sapply(b, function(x){exp(-length(a)*log(2)*x)})
		
		#input into the local_sparse matrix
		local_sparse[j, c(which(match(colnames(local_sparse), names(c))!="NA"))] <- t(c)		
		}}

#Calculate comprehensive rarity measurement
localR <- (local_distinct + local_sparse)/2

class <- matrix(NA, nrow=nrow(local_distinct), ncol=ncol(local_distinct))


#classify species per community into the 9 possibl forms of rarity
#quantiles are calculated individually for the range of values in the focal community - is this correct??

for(i in 1:nrow(local_distinct)){#for each community
	
	#calculate quantiles (using 95% and 5%)
qlS <- quantile(local_sparse[i,], probs=c(0.05, 0.95), na.rm=TRUE)
qlD <- quantile(local_distinct[i,], probs=c(0.05, 0.95), na.rm=TRUE) 
	
for(j in 1:ncol(local_distinct)){
		
#If/Else statements to classify based on comparison to the distinctness and sparseness quantile values in that community
if(is.na(local_distinct[i,j])){
	class[i,j] <- NA
	
	}else{
		
		if(local_distinct[i,j]>=qlD[2]){ #very distinct
			class[i,j] <-ifelse(local_sparse[i,j]<=qlS[1], 1, ifelse(local_sparse[i,j]>=qlS[2], 3, 2))
	
	}else{
		
		if(local_distinct[i,j]<=qlD[1]){ #not distinct
			class[i,j] <- ifelse(local_sparse[i,j]<=qlS[1], 7, ifelse(local_sparse[i,j]>=qlS[2], 9, 8))
	
	}else{
	if(local_distinct[i,j]>=qlD[1]&local_distinct[i,j]<=qlD[2])		{ #med distinct
			class[i,j] <-ifelse(local_sparse[i,j]<=qlS[1], 4, ifelse(local_sparse[i,j]>=qlS[2], 6, 5))
		}
	}}}
}}

return(list(localR, local_distinct, local_sparse, class))
}


#some plotting
local_hist<-function(class){
class<-as.vector(class)
class <- class[!is.na(class)]
class <- as.factor(class)
return(c(sum(class==1), sum(class==2), sum(class==3), sum(class==4), sum(class==5), sum(class==6), sum(class==7), sum(class==8), sum(class==9)))
}

