# MTH 362 HW 3 

drawst<-
  function (D, Mst, pl = TRUE, lab = NULL, seg = TRUE, ...) 
  {
    D <- as.matrix(D)
    st <- mst(D)
    if (is.null(colnames(D))) {
      colnames(D) <- 1:dim(Mst)[1]
    }
    if (pl) {
      lab <- colnames(D)
      plot(Mst, type = "n", xlab = "x", ylab = "y", ...)
      text(Mst, labels = lab, ...)
    }
    if (seg) {
      for (i in 1:nrow(D)) {
        w1 <- (1:nrow(D))[st[i, ] == 1]
        segments(Mst[i, 1], Mst[i, 2], Mst[w1, 1], Mst[w1, 
                                                       2])
      }
    }
  }


# Problem 2 ---- 

# Calculate B matrix and 2-variable data matrix with given Euclidean distance matrix 

D <- cbind(c(0,26,26,4,9), c(26,0,4,50,65), c(26,4,0,50,65), c(4,50,50,0,1), c(9,65,65,1,0))

# Form matrix Dc of column averages 
Dc <- t(matrix(apply(D,2,mean), 5,5))
# Form matrix Dr of row averages (transpose of Dc)
Dr <- matrix(apply(D,1,mean), 5,5)
# Form matrix Da of overall average of D 
Da <- matrix(mean(D), 5,5)
# Form B matrix using the below formula 
B <- (-1/2)*(D-Dr-Dc+Da)
B

# Find eigenvalues and vectors of B 
LB <- round(eigen(B)$values, 2)
EB <- round(eigen(B)$vectors, 2)

#Get P-matrix
P <- round(EB%*%diag(LB)^.5,0)

#Double check that P has the same distance matrix as original matrix 
P <- P[,1:2]
D1 <- as.matrix(round(dist(P)^2, 0),5,5)
D <- as.matrix(D)
D==D1 # TRUE 

# #Use cmdscale to calculate the original data set from the Euclidean distances 
# 
# #Find unsquared distance 
# Du <- round(sqrt(D), 0)
# M <- round(cmdscale(Du, k=2),1)
# 
# #Double check that we end up with D when taking squared Euclidean distance of M 
# as.matrix(dist(M)^2)

#3----

SM <- cbind(c(0,36,144,36), c(36,0, 36,16), c(144, 36, 0,36), c(36, 16, 36, 0))

#3a----
#Complete by hand 

# #Find normal Manhattan distance 
# SMu <- sqrt(SM)
# 
# #Find two-dimensional data set from SMu
# M3 <- round(cmdscale(SMu, k=2),1) # For whatever reason, the 6's that should be there are 8's 
# 
# #Double check that this dataset returns SMu
# as.matrix(dist(M3, method = "manhattan")^2) 

#3b ----
#Create matrix from worksheet
a = matrix(mean(SM[,1]), 4, 1)
b = matrix(mean(SM[,2]), 4, 1)
c = matrix(mean(SM[,3]), 4, 1)
d = matrix(mean(SM[,4]), 4, 1)

SMc = cbind(a,b,c,d)
SMr = rbind(t(a),t(b),t(c),t(d))
SMa = matrix(38, 4,4)

SMB = -.5 * (SM - SMc - SMr + SMa)

# Find eigenvalues and eigenvectors for the B matrix 

SMev <- round(eigen(SMB)$values, 2)
SME <- round(eigen(SMB)$vectors, 2)
SME <- SME[,1:2]
SMev <- SMev[1:2]

# Calculate the P matrix 

SMP <- round(SME%*%diag(sqrt(SMev)),2)

# To double-check, find Euclidean distance for SMP
as.matrix(dist(SMP)^2) 

#3c ----
library("ape")
M3C <- cbind(c(-6,0,6,0), c(0,-2,0,2))
drawst(SM, M3C, main = "Minimum Spanning Tree, 3c")

#4----
library("MASS")

#5----
# 
# D5 <- cbind(c(1,1,0,0,0,1), c(1,0,1,0,0,1), c(0,1,1,1,1,1), c(1,1,1,0,1,1), c(0,0,0,1,1,0), 
#             c(1,1,0,1,1,0),c(0,0,1,1,0,0))

D5 <- cbind(c(0,2,4,6,4,2), c(2,0,4,4,2,2), c(4,4,0,4,4,2), c(6,4,4,0,2,6), c(4,2,4,2,0,4), c(2,2,2,6,4,0))
mdD5 <- round(isoMDS(D5, k=1)$points,2)
dmdD5 <- round(dist(mdD5),2)

mdD5 <- round(isoMDS(D5, k=2)$points,2)
dmdD5 <- round(dist(mdD5),2)
#6----
D6 <- rbind(c(2,6,5,1,5), c(10,6,9,8,5), c(9,8,9,9,8), c(6,6,6,4,5), c(5, 10, 6,5,9), 
            c(5,5,4,4,4), c(6,6,6,4,4), c(4,4,4,4,4), c(4,2,5,5,5), c(5,5,5,5,5), c(4,5,5,5,5),
            c(4,4,6,9,10), c(1,2,4,4,6), c(8,5,4,5,1), c(6,6,8,6,8))

#6a----
#Compute squared distance matrix 
D6d <- as.matrix(dist(D6)^2)
D6dd <- dist(D6)

#6b----
#Obtain two-dimensional coordinates from isoMDS
mdD6 <- round(isoMDS(D6d, k=2)$points,2)


plot(mdD6, type = "n") 
text(mdD6, labels=1:15, col = 4)

#6c----

mdD6v <- Shepard(as.dist(D6d), mdD6)
mdD6v <- lapply(mdD6v, round, 2)

#calculate stress
stress <- sqrt(sum((mdD6v$y-mdD6v$yf)^2)/sum(mdD6v$y^2))

#calculate stress using MDSV 
stressMDSV <- isoMDS(D6d, k=2)$stress #THis is a percentage, so off by scale of ten 

#6d----
#Overlay minimum spanning tree 

D6dd <- as.dist(D6d)
drawst(D6d, mdD6, main = "Original Euclidean Distances")
plot(D6)

#6e----
dist_mdD6 <- dist(mdD6)
drawst(dist_mdD6, mdD6, main = "2D Euclidean Distances")

#6f----
# The two spanning trees look VERY different, actually. 

#7a---- 
D6du <- dist(D6)
D7 <- cmdscale(D6du, k=2) #Obtain two-dimensional coordinates 

#7b----
#Draw scatterplot of data
plot(D7, type = "n")
text(D7, labels=1:15,col = 4)

#Draw minimum spanning tree 
drawst(D6du, D7)

#7c ----
#Find distances for coordinates from part (a)
dD7 <- dist(D7)

#Scatterplot and overlay with minimum spanning tree 
plot(D7, type = "n",xlim=c(-7,7), ylim = c(-6,6) )
text(D7, labels=1:15,col = 4)

drawst(dD7, D7)
# The two graphs look similar, although cluster to the right is routed differently. 
# Compared to the spanning tree using the originial Euclidean distances, the "new" one 
# using the two-dimensional coordinates has 11 connecting to 9, not just 10, 6 connects to just 
# 8, not 10 and 14. 1 connects to 13 and instead of just 8 and 14 connects to 7, not 6. 

#8----

QB <- read.csv(file = "QBNFL2013.csv", row.names =1, header = TRUE)
QB

#8c----
QBs <- scale(QB, scale = TRUE)
attr(QBs, c("scaled:center")) <- NULL
attr(QBs, c("scaled:scale")) <- NULL 
QBs

#8d----
#Create distance matrix 
DQBs <- dist(QBs)
DQBsM <- round(as.matrix(DQBs),2)
colnames(DQBsM) <- c("Peyton", "Philip", "Drew", "Russell", "Tony", 
                     "Ben", "Colin", "Matt", "Alex", "Cam", "Andy",
                     "Tom", "Andrew", "Matthew", "Carson")

#Sort the distances from smallest to largest
sort(DQBsM)
# Largest distance: Peyton and Colin (7.23)
# 2nd largest: Peyton and Cam (6.71)
# Smallest distance Alex and Andrew (.82)
# 2nd smallest distance: Tony and Ben (1.04)

#8e----
QBs1 <- isoMDS(DQBs, k = 1)
QBs1$points
QBs1$stress
plot(cbind(0,QBs1$points), type = "n", ylab = "Distance", xlab = "QB name & team")
text(0, QBs1$points, rownames(DQBsM), col = 4 )

#8f----
QBs2 <- isoMDS(DQBs, k = 2)
QBs2$points
QBs2$stress
plot(QBs2$points, type = "n", main = "Two-Dimensional Data", ylab = "Score 1", xlab = "Score 2", xlim = c(-4, 5))
text(QBs2$points, rownames(DQBsM), col = 4, cex = .7 )
drawst(DQBs, QBs2$points)

#data.frame(colnames(DQBsM))
#8g----

#Sort based on Philip River's distances

PR <- DQBsM[order(DQBsM[,2]),]
PR
#Obtain distances from two-dimensional coordinates, order by PR's distances 

PRd <- dist(QBs2$points)
PRdM <- as.matrix(round(PRd, 2), 15, 15)
PRdM <- PRdM[order(PRdM[,2]),]

#For the most part, the order is relatively preserved, although players in the middle are garbled. The top 5 most 
# similar players are the same, with Tony more alike than Ben. The middle 7 are jumbled but in similar places 
# essentially, Tom and Peyton moved up. The three most dissimilar player's distances are preserved 

#8h----

#To answer this problem, I created a for-loop that would print off the dataset ordered by each 
#player. I then visually inspected the datasets to see which ones met the question criteria. 

colnames(PRdM) <- c("Peyton", "Philip", "Drew", "Russell", "Tony", 
                        "Ben", "Colin", "Matt", "Alex", "Cam", "Andy",
                        "Tom", "Andrew", "Matthew", "Carson")
for (i in 1:ncol(DQBsM)){
  Delete1 <- DQBsM[order(DQBsM[,i]),]
  print(data.frame(Delete1[,i]))
  Delete2 <- PRdM[order(PRdM[,i]),]
  print(data.frame(Delete2[,i]))
}

# Carson Palmer: all ranks are incorrect 

#8i----
#Ben R, all correct except for Luck and Rivers

#8j ----
md <- data.frame(mahalanobis(QBs, center = 0, cov =cov(QBs)))
# Two largest distances: Peyton and Russell 

