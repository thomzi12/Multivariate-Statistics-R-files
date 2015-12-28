# Necessary functions for the homework ----
library(aplpack)
pairsbg<-
  function (M, bg = TRUE,...) 
  {
    if (is.null(rownames(M))) {
      rownames(M) <- 1:(dim(M)[1])
    }
    pairs(M, panel = function(x, y) {
      if (bg) {
        bagplot(cbind(x, y), add = TRUE,show.whiskers=FALSE)
      }
      text(x, y, abbreviate(row.names(M)), cex = 0.8,...)
      
    })
  }

getshads<-
  function (M,clM,i,j)
  {
    L<-1:dim(M)[1]
    DM<-dist2(M,parameters(clM))
    DMo<-t(apply(DM,1,order))
    DMoi<-DMo[clusters(clM)==i,]
    ni<-dim(DMoi)[1]
    v<-L[(DMo[,1]==i)&(DMo[,2]==j)]
    dii<-DM[(DMo[,1]==i)&(DMo[,2]==j),i]
    dij<-DM[(DMo[,1]==i)&(DMo[,2]==j),j]
    d<-2*dii/(dii+dij)
    ans<-NULL
    ans$v<-v
    ans$ni<-ni
    ans$dii<-dii
    ans$dij<-dij
    ans$shadows<-2*dii/(dii+dij)
    ans$similarity<-sum(d)/ni
    ans
  }



sumcldist<-
  function (M,cl) {
    mv<-NULL; Dt<-0
    t<-unique(cl)
    for (i in t) {
      U<-cbind(M[cl==i,])
      m<-colMeans(U)
      mv<-rbind(mv,m)
      Dm<-dist2(U,m); print(Dm)
      Dt<-Dt+sum(Dm) }
    rownames(mv)<-t
    list(mv[order(t),],Dt)}

#1----
D = cbind(c(0,8,9,13,18,16), c(8,0,1,5,10,8), c(9,4,0, 4,9, 7), c(13,5,4,0,5,3), 
          c(18, 10, 9,5,0,2),c(16,8,7,3,2,0))

ZT <- cbind(c(0,8,9,13,18,16)) #Guessed the one-dimensional coordinates 
rownames(ZT)<-ZT

#1g----
library(flexclust)
#Repeat command until match is found 
M.cl3<-cclust(cbind(0,ZT),k=3,save.data=TRUE);paste(clusters(M.cl3),collapse="") #"122233" 

#Get cluster similarity matrix 
round(clusterSim(M.cl3), 3)
#Get Shadow stars plot 
shadowStars(M.cl3)

#2----

#2a----
MQ2 <- read.csv(file = "HW4Q2.csv", row.names =1, header = TRUE)
head(MQ2)
#2b----
library(scatterplot3d)
scatterplot3d(MQ2, main = "MQ2 3D Scatterplot")

#2c----
#Plot using plot3d and text3d 
library(rgl)
plot3d(MQ2)


#2d----
SFg <- stepFlexclust(MQ2, nrep=10, k = 2:20); plot(SFg, type = "l", ylab = "sum of within cluster distances")
# Skree plot says 6 clusters 

#2e----
# Use while loop to find a minimum intercluster distance 
t <- 445
while(t > 312.03){MQ2.cl3 <- cclust(MQ2,k=6,save.data=TRUE); t<- (info(MQ2.cl3, "distsum")); print(t)}

#Plot the cclust run 
plot(MQ2.cl3, project = princomp(MQ2), main = "Cclust Run") #Clusters projected down to two dimensions 

#2f----
#Plot 3D scatterplot with clustered labeled and in different colors 
plot3d(MQ2, main = "MQ2 3D scatterplot, with labeled clusters", type = "n")
text3d(MQ2, texts = clusters(MQ2.cl3), color = clusters(MQ2.cl3)+1, cex = 0.8)

#2g----
# Obtain barplot for the cclust result 
barplot(MQ2.cl3) 

#2h----
round(centers <- MQ2.cl3@centers,0) 
#Hmm ... these are six points that are vetrices for a cube 

#2i----
# Get similarity matrix for cclust values from 2e
MQ2sm <- round(clusterSim(MQ2.cl3, symmetric = FALSE), 3)

#2j----
#Get shadow values for each cluster
shadow(MQ2.cl3)

## Below no longer needed for hw since not possible ##
# These values tell us how "concentrated" together each cluster is, with a smaller 
# value indicating a cluster that is more apart from other on the plot. For example, 
# cluster 1 has a high cumulative shadow value, and consequently overlaps heavily with 
# cluster 4 and 5. 5 has a low cluster shadow value and only overlaps with 1. 

#2k----
#Find which centroids are closest to 0,0,0 and which one is closest to 1,1,0
centers <- MQ2.cl3@centers
dc <- rbind(centers, c(0,0,0), c(1,1,0))
dist(dc)
# For point 0,0,0, the closest cluster is 5 (0.0608)
# For point 1,1,0, the closest cluster is 3 (0.1197)

#Find the set of points Aa,b whose closest cluster is 6, second closest 3
shads65 <- getshads(MQ2, MQ2.cl3, 5,3) # 9  49  91 211 511
shads65$v
#2l----
# Verify value of MQ2sm(5,3) = (sum of shadow values of points closest to centroid 5, next closest to 3)/
    #(number of points in cluster 5)

shads65[[5]] # Shadow values for each point 
sumshads65 <- sum(shads65[[5]]) # Sum of shadow values for these points 
clust6points <- shads65[[2]]
smval65 <- sumshads65/clust6points
sm65 <-MQ2sm[5,3]
smval65
sm65
all.equal(sm65, smval65) #pretty much the same 

#2m----
# Locate the lines for the shadow values for Aa,b in a shadow stars plot 
getshads(MQ2, MQ2.cl3, 5,3)$shadows #Jeez, these are brutal numbers 
shadowStars(MQ2.cl3, project=prcomp(MQ2), width=2, box=1) # 5-3 shadow values are on bar to the left 

#2n----
#dendrograms 
plot(hclust(dist(MQ2), method = "single"), main = "Cluster Dendrogram, Singly-linked")
plot(hclust(dist(MQ2), method = "complete"), main = "Cluster Dendrogram, Complete")
plot(hclust(dist(MQ2), method = "average"), main = "Cluster Dendrogram, Average")

#2o----
# Find within cluster distances ... ?
MQ2s <- cutree(hclust(dist(MQ2), method = "single"), k = 6)
d1 <- sumcldist(MQ2, MQ2s) #597.0728 ... which is surprising 
d1[[2]]
MQ2a <- cutree(hclust(dist(MQ2), method = "average"), k = 6)
d2<- sumcldist(MQ2, MQ2a) #327.3163
d2[[2]]
MQ2c <- cutree(hclust(dist(MQ2), method = "complete"), k = 6)
d3 <- sumcldist(MQ2, MQ2c) #338.1586
d3[[2]]

#3----
MQ3 <- read.csv(file = "HW4Q3.csv", row.names =1, header = TRUE)
head(MQ3)

#3a----
MQ3 <- MQ3[,c(2,5,6)]

#Scale the crime data 
MQ3 <- scale(MQ3)

#3b----
pairsbg(MQ3)

#3c and 3d ----
#Find the mahalanobis distances for the scaled dataset 
MQ3md <- mahalanobis(MQ3, center = 0, cov = cov(MQ3))
MQ3md <- data.frame(cbind(sort(MQ3md), seq(1,51,1)))
# Remove UT and AK from MQ3
MQ3 <- data.frame(MQ3, seq(1,51,1))
MQ3 <- MQ3[-c(45,50),-4] 
MQ3 <- scale(MQ3) 
head(MQ3)  # matches worksheet 

#3e----
#Obtain principal components for MQ3 
MQ3cov <- cov(MQ3)
MQ3eig <- eigen(MQ3cov)
MQ3eiv <- MQ3eig$vectors # eigenvectors for MQ3 
MQ3pc <- MQ3 %*% MQ3eiv #principal components for MQ3 

#Interpretations for principal components: 
# PC1: states that score highly in PC1 have low levels of crime in these three categories
# PC2: states that score highly in PC2 have high levels of rape but low levels of burglary and theft 
# PC3: states that score highly in PC3 have low levels of burglary but high levels of rape and theft 

#3f----
library(bpca)
plot(bpca(MQ3, method = "gh"), main = "MQ Biplot")
# It should be noted that since this is a (P(lambda)^(-1/2), J) biplot, correlations and projects 
# are preserved in this biplot, while the distances between points are not quite preserved. 
# More generally: States like TX, WA, and OR and FL due to their placement ahead of the principal components 
# for the crimes would appear to be high-crime states, whereas states opposite of the principal components, like 
# KY, PA, WV would appear to be low-crime. WV especially seems to be low theft. 


#3g----
plot(hclust(dist(MQ3), method = "complete"), main = "Cluster Dendrogram, Complete")
#Hmm, tough choice, given this dendrogram. If our cutoff distance is 3, we have four clusters, 
# which seems reasonable to me. If our cutoff distance is 2, then 6 clusters appears. 

#3h----
MQ3c <- cutree(hclust(dist(MQ3), method = "complete"), k = 4)
data.frame(MQ3c)

# Cluster 1: ME, NH, VT, PA, IN, WI, IA, ND, SD, NE, VA, WV, KY, MT, ID, WY
# C2: MA, RI, CT, NY, NJ, OH, IL, MN, MO, KS, NC, TN, AL, MS, AR
# C3: DE, MI, NV, MI, MD, SC, GA, LA, OK, CA, HI
# C4: FL, TX, OR, CO, AZ, DC, WA, DC

#The sum within cluster distances 
del <- sumcldist(MQ3, MQ3c) # 35.36096
del[[2]]

#3i----
#4-clustering of data using cclust with within cluster distances < 34.30
f <- 445

while(f > 34.30){MQ3.cl4 <- cclust(MQ3,k=4,save.data=TRUE); f<- (info(MQ3.cl4, "distsum")); print(f)}

# Similarity matrix with exactly six zeros 
round(clusterSim(MQ3.cl4, symmetric = FALSE), 1)
#States in  each cluster 
sort(clusters(MQ3.cl4))

#3j-----
plot(MQ3.cl4, points = FALSE, main = "Cclust of MQ3")
text(MQ3, rownames(MQ3), cex = 0.7, col = 4)


#Plot with projection to principal components 
plot(MQ3.cl4, points = FALSE, main = "Cclust of MQ3, w/ PC projection", project = prcomp(MQ3))
text(prcomp(MQ3)$x, rownames(MQ3), cex = 0.7, col = 4)  

# The projection to the principal components of MQ3 visually linearizes the data, thus allowing
# us to better interpret it. Roughly, states on the right side of the line, with a higher Pc1 
# value, are high crime as this dataset measures it, where low-scoring states on the left side of the 
# line are low-crime. 

#3k----
#4-clustering of data using cclust with within cluster distances > 39.30
z <- 39
while(z < 39.5){MQ3.cl4n <- cclust(MQ3,k=4,save.data=TRUE); z<- (info(MQ3.cl4n, "distsum")); print(z)}
# z = 39.72294

plot(MQ3.cl4n, points = FALSE, main = "Cclust of MQ3")
text(MQ3, rownames(MQ3), cex = 0.7, col = 4)

#Plot with projection to principal components 
plot(MQ3.cl4n, points = FALSE, main = "Cclust of MQ3, w/ PC projection", project = prcomp(MQ3))
text(prcomp(MQ3)$x, rownames(MQ3), cex = 0.7, col = 4) 

#3l----
barplot(MQ3.cl4) 
#Give centroids 
MQ3.cl4@centers
#The different clusters essentially illustrate different tiers of crime: cluster 1 (c1) has states with low 
# crime in these three categories, whereas C2 to a lesser extent has less crimes in these categories, C3 is 
# almost the opposite of C2, and C4 the opposite of C1: states in this category are high crime in those areas

#3m----
# Create shadowstars plot and locate the lines corresponding to given states 
shadowStars(MQ3.cl4, project=prcomp(MQ3), width=2, box=1)
sort(clusters(MQ3.cl4))
sort(getshads(MQ3, MQ3.cl4, 4,1)$shadows)  
sort(getshads(MQ3, MQ3.cl4, 1,2)$shadows)
sort(getshads(MQ3, MQ3.cl4, 1,4)$shadows)
sort(getshads(MQ3, MQ3.cl4, 2,1)$shadows)
sort(getshads(MQ3, MQ3.cl4, 2,3)$shadows)
sort(getshads(MQ3, MQ3.cl4, 3,2)$shadows)



#4-----

DB4 <- read.csv(file = "HW4Q4.csv", row.names =1, header = TRUE)
#4a -----
data.frame(rownames(DB4)) #Minnesota (MN) is missing! 

#4b----
DB4 <- scale(DB4)

#4c-----
pairsbg(DB4) 

#4d ----
#Find the mahalanobis distances for each state in scaled data set 
DB4md <- mahalanobis(DB4, center = 0, cov = cov(DB4))

#4e----
DB4cov <- cov(DB4)
DB4eig <- eigen(DB4cov)
DB4ev <- DB4eig$vectors

#Principal components 
DB4pc <- DB4 %*% DB4ev

#Give interpretations of the first two principal components
# PC1: States that score highly in PC1 have low amounts of crime in these 5 categories
# PC2: states that score highly in PC2 have high violent crime and murder rates, but 
# are low in Forgery and embezzlment, and to an extent fraud. 

#4f----
# Plot (P,J) biplot (method ="hj")
plot(bpca(DB4, method = "hj"), main = "DB4 (P,J) Biplot")


#4g----
#Produce four stepflexclust plots, to determine what's reasonable for k 
for (i in 1:4){
DB4s <- stepFlexclust(DB4, nrep=30, k = 2:20); plot(DB4s, type = "l", ylab = "sum of within cluster distances")
}

#Honestly, tough to say. However, looking at the plots, 4 or 5 clusters would seem to be a safe interpretation
# of the descriptive data at hand. 

#4h----
#Obtain 4-clustered of set using cclust. Ensure total within cluster distance is below 65.35. 

y <- 70
while(y > 65.35){DB4.cl4 <- cclust(DB4,k=4,save.data=TRUE); y<- (info(DB4.cl4, "distsum")); print(y)}

#List the states in each cluster. 
sort(DB4.cl4@cluster)

#4i----
DB4 <- DB4[,]

plot(DB4.cl4, points = FALSE, main = "Cclust of dDB4, w/ PC projection", project = prcomp(DB4))
text(prcomp(DB4)$x, rownames(DB4), cex = 0.7, col = clusters(DB4.cl4))


#4j----
barplot(DB4.cl4)
#Get cluster centers
DB4.cl4@centers

#4k----
#Get similarity matrices, using both values of `symmetric`
round(clusterSim(DB4.cl4), 3)
round(clusterSim(DB4.cl4, symmetric = TRUE), 3)

#4l----
shadowStars(DB4.cl4, project=prcomp(DB4), width=2, box=1)

#4m----
#Discuss the structure of the data referring back to the cclust results from (c) - (l)

# I don't believe this is necessarily the answer you were looking for, but I think this dataset 
# does a good job of illustrating the messy nature of real-life data. For one thing, the pairs plot 
# would seem to imply that many of the crime variables are positively correlated, and thus colinearity may 
# be an issue in this data. So that's one thing. If you look at the distances between the points, they rise 
# in a steady gradient but as you go up towards GA, LA, AL, NC, DE, MS, AK, and IL, the distances grow large 
# and almost appear as outliers. The first two PCs only differ really when it comes to violent crime and murder 
# rates. Indeed, looking at the biplot, most points are indiscriminately centered around the origin, except for
# a couple high in embezzlement/forgery or violent crime/murder, or both (both, incidently, are many of the 
# states with high distances). Looking at the Skree plot, I know that 4 or 5 plots is a good interpretation, 
# but there's not really a "elbow" to the plot and it's hard to say how many clusters there should truly be. 
# Indeed, the shadow stars plot show most points are between between two clusters, except for those in cluster 
# 4, those relatively low in crime. 

# In summary, this data's structure includes some states low in crime, a majority in the middle, and 
# some crime-ridden outliers. Further review might involve splitting the data set by the murder/violent crime
# and forgery/fraud/embezzlement variables and repeating the analysis. 



