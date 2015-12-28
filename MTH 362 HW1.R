# MTH 758 homework ----
gradMat <- cbind(c(.75, .25), c(.25, .75))
delete <- gradMat %^% 200
delete1 <- c((1/6), (5/6))

# MTH 362 Hw 1 ----

#1 ----
a = cbind( c(.5, .1), c(.5, .9))
b = c((1/6), (5/6))
data(USairpollution)
USair <- USairpollution

c = cbind(c(8,5,2,8,2), c(6,4,8,8,4), c(5,9,7,5,9))
M = cbind(c(3,0,-3,3,-3), c(0,-2,2,2,-2), c(-2,2,0,-2,2)) # after "removing mean" of each column
Mm = t(M)
covM <- Mm %*% M 
cov_final <- covM/4 #cov(M)
D = cbind(c((1/3),0,0), c(0,.5,0), c(0,0,.5))
inte = D %*%cov_final
cor_M <- inte %*% D # cor(M)


c = cbind(c(8,5,2,8,2), c(6,4,8,8,4), c(5,9,7,5,9)) 
cor(c) # cor(M)
Mu <- scale(M, scale= TRUE)
cov(Mu) # cov(Mu) = cor(M)


# 2 ----

mat2 = cbind(c(7,1,3,5,9), c(6,8,10,6,10), c(8,4,4,6,8))
colnames(matz) <- c("x", "y", "z")
sub <- mat2[,c(1,3)]
Bp <- bagplot(sub, show.whiskers=FALSE, main = "Bagplot, Variables x and z", 
              xlab = "x", ylab = "z", ylim = c(4, 8.2))
text(sub, labels = Bp$hdepths)
Bp$hdepths # [1] 1 1 1 2 1

#3 ----

set3 = cbind(c(.9, 1, .8, 2, 2.9, 4, 2, 1.7), c(4.9, 3, 1, 2.9, 2.1, 3, 1.9, 4))
colnames(set3) <- c('x', 'y')
rownames(set3) <- c(1,2,3,4,5,6,7,8)
plot(set3, pch = 19)

# saved at 'Scatterplot' in pics or documents 

#4---- 

bp4 <- bagplot(set3, show.whiskers=FALSE, main = "Bagplot (Index/Line Depth)")
text(set3, labels = paste(rownames(set3), "/", bp4$hdepths))
bp4$hdepths #[1] 1 2 1 3 1 1 2 2
# four points with the largest line depth: 4, 8, 2, and 7

#5 ----
chull(set3) #[1] 5 3 1 6

#6 Find correlation between the two variables before and after removing the convex hull ----

cor(set3) # -0.09604193
cor(set3[-chull(set3),]) # -0.3005184

#7 ---- 

z <- c(1.1, 2.0, 1, 4, 5, 10, 4.1, 3)
set7 <- cbind(set3, z)
plot(set7, xlab = "x", ylab = "y", pch = 3, main = "Bubbleplot")
symbols(set7, circles = set7[,3], inches = 0.5, add= TRUE)

#8----

pairs(set7, main = "Pairs Plot")

#9 ----

set9 <- set7[, c(1,3)]
bnd1 <- dpik(set9[,1]) # [1] 0.6633154
bnd2 <- dpik(set9[,2]) # [1] 1.304782

#10 
set10 = c(.9, 1, .8, 2, 2.9, 4, 2, 1.7, 4.9, 3, 1, 2.9, 2.1, 3, 1.9, 4, 1.1, 2, 1, 4, 
          5, 10, 4.1, 3, 9.6, 9.8, 9.7, 11)

# Test different values of h, the bandwidth: 
plot.new()
pm(3, 2)

h <- 0.1 
plot(density(set10, bw = h), main = paste("H =",h), col = 62, lwd=2)
rug(set10, lwd = 1, ticksize = .3)

h <- 0.2 
plot(density(set10, bw = h), main = paste("H =",h), col = 62, lwd=2)
rug(set10, lwd = 1, ticksize = .3)

h <- 0.4 
plot(density(set10, bw = h), main = paste("H =",h), col = 62, lwd=2)
rug(set10, lwd = 1, ticksize = .3)

h <- 0.5 
plot(density(set10, bw = h), main = paste("H =",h), col = 62, lwd=2)
rug(set10, lwd = 1, ticksize = .3)

h <- 0.8
plot(density(set10, bw = h), main = paste("H =",h), col = 62, lwd=2)
rug(set10, lwd = 1, ticksize = .3)

h <- 3 
plot(density(set10, bw = h), main = paste("H =",h), col = 62, lwd=2)
rug(set10, lwd = 1, ticksize = .3)

# Determine a reasonable value for h '
#Based on the graphs produced, I would say either a relatively high h of .5 showing a 
# bimodal distribution would be best for this distribution. Let's see what dpik says... 
h = dpik(set10)
# h = 0.8297808

#11
# Detemine which combination of variables produces Albuquerque, Phoenix, Hartford, and Providence
# as outliers. 

pm(1,1)
for( i in 1:6){
  j = i 
  for (j in j:6){
    # print(paste(c(i, j+1)))
  bagplot(USair[,c(i, j+1)], show.whiskers = FALSE, main = paste("Combination:", i, j+1))
  text(USair[,c(i, j+1)], labels = rownames(USair), cex = .8)
  }
}

#Ah, it's the combination of variables 2 and 6, i.e. temperature and precipitation! 
# (image shows cities on the convex hull)

bg11 <- bagplot(USair[,c(2, 6)])
cities <- chull(USair[,c(2, 6)])
for (i in 1:length(cities)){
  print(rownames(USair)[cities[i]])
}

#12 

for( i in 1:6){
  j = i 
  for (j in j:6){
    # print(paste(c(i, j+1)))
    temp <- compute.bagplot(USair[,c(i, j+1)], show.whiskers = FALSE, main = paste("Combination:", i, j+1))
    
  }
}

outl <- bagplot(USair[,c(1, 5)])
outl$pxy.outlier
