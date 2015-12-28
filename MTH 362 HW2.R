# MTH 362 homework 
# For matrix multipication, matmult package

# Pairsbg ----
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

#2 ---- Consider this data for 15 students 

C = cbind(c(2, 10, 9, 6, 5, 5, 6, 4, 4, 5, 4, 4, 1, 8, 6), c(6, 6, 8, 6, 10, 5, 6, 4, 2, 5, 5, 4, 2, 5, 6), 
          c(5, 9, 9, 6, 6, 4, 6, 4, 5, 5, 5,6, 4, 4, 8), c(1,8,9,4,5,4,4,4,5,5,5,9,4,5,6),
          c(5, 5,8,5,9, 4, 4, 4, 5, 5, 5, 10, 6, 1, 8))
# (a)----
# Use loop to calculate Ms 
Ms <- C
for (i in 1:ncol(Ms)){
  sum = 0
  for (j in 1:nrow(Ms)){
    sum = sum + Ms[j,i]
  }
  avg = sum/nrow(Ms)
  print(avg)
  for (k in 1:nrow(Ms)){
    Ms[k,i] = Ms[k,i] - avg
  }
}
Mst <- t(Ms)
C4 <- Mst %*% Ms
# Find covariance matrix 
S = cov(C)
C5 <- C4/(nrow(Ms)-1)
#Double check that manual calculations and cov(c) are equal
all.equal(C5, S) # TRUE 

# (b)----
# Find covariance matrix, S, then E matrix of eigenvectors, and eigenvalues e. 
S = cov(C) 
e <- eigen(S)
e_m <- matrix(e$vectors, 5)
ep <- t(e_m)
e1 <- e_m %*% ep
e2 <- ep %*% e_m 

all.equal(e1, e2) # TRUE, EE' = E'E = I 

ev <- e$values
ev <- matrix(ev)
lambda = diag(5)
diag(lambda) <- ev # lambda now has eigenvalues as its diagonal 

s2 <- e_m %*% lambda %*% ep
all.equal(s2, S) #TRUE, S = E(Lambda)E'

round(s2, 6)
#(c)----
P = Ms %*% e_m
covP <- diag(cov(P))
all.equal(covP, ev[,1]) # TRUE, cov(P) = eigenvalues

screeob <- princomp(C)
screeplot(screeob, type = "l", main = "Scree Plot of C") #Saved as Screeplot in MTH362HW2 folder 

#(d)----
# For PC1, S3 had the lowest score (-7.20) and S13 had the highest (4.77)
# For PC2, S12 had the lowest score (-4.76) and S14 had the highest (5.31)
# For PC3, S12 had the lowest score (-2.80) and S5 had the highest (4.19)

# Interestingly, you can look at the eigenvectors to interpret the PC scores. 
# Students scoring high in PC1 have low scores, in particular for test 1 and 4. More generally, 
# lower scores indicate a better students.
# Students scoring high in PC2 have high score on test 1 and low scores on the last test. More generally, 
# higher PC2 scores indicate the student did worse as the semester progressed. 
# Students scoring high in PC3 did well on test 2 but poorly on test 4. 
# Students scoring high in PC4 did poorly on test 3, alright on tests 2 and 3. 
# Students scoring high in PC5 did good on the first and last tests, but not as well on the tests in between. 

#(e)----
# You cannot literally take the negative square root of S, since you have negative values present. 
# Rather, knowing that S = ELE', S^(-.5) = E(1/sqrt(L))E'. 
ev2 <- sqrt(ev)
ev2 <- 1/ev2
lambda2 <- diag(5)
diag(lambda2) <- ev2
s_prime <- e_m %*%lambda2%*%ep 
Ns  = Ms%*% s_prime 
all.equal(cov(Ns), diag(5)) #TRUE  Ns = M*S^(-.5)

#(f)
pairsbg(C, bg=FALSE, main = "f") # Saved in folder 

#(g)
library(stats)
md <- data.frame(mahalanobis(Ms, center = 0, cov =S ))
names(md) <- "Dist"
index <- t(t(seq(1,15,1)))
md <- cbind(index, md)
md %>% select(index,Dist) %>% arrange(desc(Dist))
# Based on the current settings, the three individuals with the highest 
# generalized distances are S5, S12, S14 

# If you look at only two variables at a time: 
for(i in 1:(ncol(Ms)-1)){
  j = i
  loopl <- NULL 
  for (j in j:(ncol(Ms)-1)){
    print(i)
    print(j+1)
    loopl <- data.frame(mahalanobis(Ms[,c(i,j+1)], center = 0, cov =cov(Ms[,c(i, (j+1))] )))
    names(loopl) <- "Dist"
    loopl <- cbind(index, loopl)
    print(t(loopl %>% select(index,Dist) %>% arrange(desc(Dist))))
  }
  
}

#Print off particular combinations 
i = 2
j = i
print(i)
print(j+1)
loopl <- data.frame(mahalanobis(Ms[,c(i,j+1)], center = 0, cov =cov(Ms[,c(i, (j+1))] )))
names(loopl) <- "Dist"
loopl <- cbind(index, loopl)
print((t(round(loopl %>% select(index,Dist) %>% arrange(desc(Dist)), 2))))

Nst = t(Ns)
dists <- Ns%*% Nst
dists <- data.frame(diag(dists))
dists <- cbind(index, dists)
dists %>% select(index, diag.dists.) %>% arrange(desc(diag.dists.))

#(h)----
library(depth)
depths <- data.frame(apply(Ms,1,depth, x=Ms, approx = TRUE)*15)
colnames(depths) <- "Depths"
t(depths)

#(i)----
egv <- sqrt(ev)
srl <- diag(5)
diag(srl) <- egv

j = e_m %*% srl
jt = t(j)
all.equal(S, (j%*%jt)) # JJ' = S

round(j, 4)
Ps = P[,c(1,2)]

library(bpca)
plot(bpca(Ms, method = "hj"), var.scale = FALSE, main = "P,J Biplot for 2(i)")  

#(j)----
jw = j[, c(1,2)]
jwt = t(jw)
sw = jw %*% jwt# Although off by .05 to 1.5+ when compared to S, this is a decent approximation 
sw
jh = j[, c(1,2,3)]
jht = t(jh)
sht = jh %*% jht
sht

#(k)
plot(bpca(Ms, d = 1:3, method = "hj", scale = FALSE), rgl.use = TRUE, main = "P,J Biplot for 2(k)")  


#(l)
round(dist(Ms), 1)
# can't quite comment on bpca plot yet 

#(m)
msm <- P %*% ep
all.equal(msm, Ms) #TRUE, PE' = Ms

Est = t(ep[,c(1,2)])
msms = Ps %*% Est
round(msms, 2) # As far as I can tell, not very close to Ms
round(Ms, 2)
round(msm, 2)

Ph = P[,c(1,2,3)]
Eht = t(ep[,c(1,2,3)]) 
Msh = Ph %*% Eht
round(Msh, 2) # Still not a very good approximation of Ms. Those last two columns most hold a bunch of variation! 

#3----
carData <- read.csv("C:/Users/thomzi12/Documents/MTH362 HW2/CarRoom.csv")
cData <- as.matrix(carData)
carnames <- cData[,1]
cData <- cData[,-1]
rownames(cData) <- carnames
cData <- as.numeric(cData)
cData <- matrix(cData, 20)
cData <- scale(cData)# Scaled Data. Means for columns: 42.075  4.650 29.175  3.550
round(cData, 2)

cD_cor <- cor(cData) # Correlation matrix 
round(cD_cor, 4)
cD_e <- eigen(cD_cor)
cD_eval <- cD_e$values
cD_eval
cD_evec <- cD_e$vectors
cD_evec
#(d)
cD_PCs <- cData %*% cD_evec #Find PC scores 
# The first PC is a general measure of how roomy a car is, although especially regarding headspace throughout the car
## The car that scores the best is the C. Silverado, the car that scores the worst is the C. Spark. This makes sense 
## since you are comparing a large pick-up truck to a smart car. 
# The second PC scores highly cars with plenty of leg space but relatively less head space. 
## The vehicle with the highest PC2 score is VW Passat, lowest PC2 score goes to the C. Spark. 
# The third PC score scores highly cars that are more spacious in the front than rear. 
## The vehicle with the highest PC3 score is the C. Sonic, and the lowest scoring is the Lexus LS. The Ls is pretty spacious!

#(e) 
# I assume that we are using the scaled values here 
plot(bpca(cData, method = "hj"), main = "hj Method Biplot, 3(e)") #stored as 3E biplot 
plot(bpca(cData, method = "jk"), main = "jk Method Biplot, 3(e)") #stored as 3E biplot 
plot(bpca(cData, method = "gh"), main = "gh Method Biplot, 3(e)") #stored as 3E biplot 

#(g)---- 
car_m<- data.frame(mahalanobis(cData, center = 0, cov =cov(cData)))
index2 <- t(t(seq(1,20,1)))
car_m <- cbind(index2, car_m)
colnames(car_m) <- c("Index", "Dists")
car_m %>% select(Index, Dists) %>% arrange(desc(Dists))

carsub <- carData[c(4,7,11),]
carsub
#4----
NBAtm <- read.csv("C:/Users/thomzi12/Documents/MTH362 HW2/NBAtopmin2014.csv")
NBAnames <- NBAtm[,1]
NBAtm <- NBAtm[,-1]
rownames(NBAtm) <- NBAnames

#(a)
cNBAtm <- cor(NBAtm) # Correlation matrix of NBAtm data 

#(b)
NBAtme <- eigen(cNBAtm)
NBA_eval <- NBAtme$values
NBA_evec <- NBAtme$vectors

#Find P 
sNBA <- scale(NBAtm)
NBAp <- sNBA %*% NBA_evec

#(c)
plot(bpca(sNBA, method = "gh"), main = "gh Method Biplot, 4(c))")

#(e)
NBA_d<- data.frame(mahalanobis(sNBA, center = 0, cov =cov(sNBA)))
index2 <- t(t(seq(1,18,1)))
NBA_d <- cbind(index2, NBA_d)
colnames(NBA_d) <- c("Index", "Dists")
NBA_d %>% select(Index, Dists) %>% arrange(desc(Dists))

