# initialize required library
library(ca)

par(ask=T) # ask before initializing new plot

# read in data table and conduct CA
mydata <- read.table(file=file.choose(),sep=',',header=T,row.names=1)
mydata <- as.matrix(mydata)
data.ca <- ca(mydata)

# define variables
rown <- nrow(mydata) # number of rows in data table
coln <- ncol(mydata) # number of columns in data table
nsim <- 1000 # number of bootstrapped simulations to run
cut <- 0.90 # set convex hull confidence interval for final plot
lab <- colnames(mydata) # define column (type) labels

# create nsim bootstrapped replicates from original data
data.rowsum <- apply(mydata,1,sum)
data.sim <- rmultinom(nsim,data.rowsum[1],prob=mydata[1,]) # create bootstrapped replicate for first row with probability defined by actual proportions of types
for (i in 2:rown) {
data.sim <- cbind(data.sim,rmultinom(nsim,data.rowsum[i],prob=mydata[i,]))} # create bootstrapped replicates for remaining rows

# reorder bootstrapped replicates by number of rows
data.sim <- t(data.sim)
data.sim2 <- matrix(rep(0,nsim*rown*coln),nrow=nsim*rown)
for (k in 1:nsim) {
	for (i in 1:rown) {
		data.sim2[(k-1)*rown+i,] <- data.sim[k+(i-1)*nsim,]}}

# define simulated column principal coordinates using row transition formula
data.rowsc <- data.ca$rowcoord[,1:2]
data.colsim <- t(t(data.rowsc) %*% data.sim2[1:rown,]) / apply(data.sim2[1:rown,],2,sum)
for (k in 2:nsim) {
data.colsim <- rbind(data.colsim, t(t(data.rowsc) %*% data.sim2[((k-1)*rown+1):(k*rown),]) / apply(data.sim2[((k-1)*rown+1):(k*rown),],2,sum))}
data.colsim2 <- matrix(rep(0,nsim*coln*2),nrow=nsim*coln)

# define simulated row principal coordinates using column transition formula
data.colsc <- data.ca$colcoord[,1:2]
data.rowsim <- t(t(data.colsc) %*% t(data.sim2[1:rown,])) / apply(data.sim2[1:rown,],1,sum)
for (k in 2:nsim) {
x <- t(t(data.colsc) %*% t(data.sim2[((k-1)*rown+1):(k*rown),]))
data.rowsim <- rbind(data.rowsim, x / apply(data.sim2[((k-1)*rown+1):(k*rown),],1,sum))}
data.rowsim2 <- matrix(rep(0,nsim*rown*2),nrow=nsim*rown)

# populate row and column simulated point matrices with projected coordinates
for (j in 1:coln) {
	for (k in 1:nsim) {
		data.colsim2[(j-1)*nsim+k,] <- data.colsim[j+(k-1)*coln,]}}
for (j in 1:rown) {
	for (k in 1:nsim) {
		data.rowsim2[(j-1)*nsim+k,] <- data.rowsim[j+(k-1)*rown,]}}

# plot original CA rows and columns
plot(data.ca) 

# set up matrices for recording convex hull sizes
x2 <- matrix(0,nrow(mydata),1)
x3 <- matrix(0,nrow(mydata),1)
x2_a <- matrix(0,nrow(x2),1)
x3_a <- matrix(0,nrow(x3),1)
row.names(x2) <- row.names(mydata)
row.names(x3) <- row.names(mydata)
row.names(x2_a) <- row.names(x2)
row.names(x3_a) <- row.names(x3)

# Define convex hulls around row points and record hull sizes
poly.row <- list()
for (j in 1:nrow(mydata)) {
	points <- data.rowsim2[(nsim*(j-1)+1):(nsim*j),]
	hpts <- chull(points)
	hpts <- c(hpts,hpts[1])
	v <- points[hpts,]
	poly.row[[j]] <- v
	v <- rbind(unique(v),v[1,])
	v2 <- as.matrix(dist(v))
	x <- matrix(0,nrow(v2)-1,1)
	for (i in 1:nrow(x)) {
	x[i,1] <- v2[i+1,i]}
	x2[j,1] <- sum(x)
	x3[j,1] <- max(v2)}

# Define convex hulls around column points
poly.col <- list()
for (j in 1:coln) {
	points <- data.colsim2[(nsim*(j-1)+1):(nsim*j),]
	hpts <- chull(points)
	hpts <- c(hpts,hpts[1])
	poly.col[[j]] <- points[hpts,]}

# plot all bootstrapped data
plot(data.rowsc[,1],data.rowsc[,2],type='n',xlab='CA1',ylab='CA2',main='All Simulated Data')
data.col <- t(t(data.rowsc) %*% mydata) / apply(mydata,2,sum) # set column label positions
# use last two digits of color code to set transparency
for (i in 1:length(poly.row)) { 
polygon(poly.row[[i]],col='#0000ff05',border='#0000ff20')} 
for (i in 1:length(poly.col)) { 
polygon(poly.col[[i]],col='#ff000050',border='#ff000050')
text(data.col[i,1],data.col[i,2],lab[i],font=2,cex=0.5)}

# plot histogram of convex hull perimeter length and set cutoff
hist(main='Convex Hull Perimeter Lengths',x2,breaks=25,col='blue')
abline(v=mean(x2)+(sd(x2[,1])),col='red')
abline(v=mean(x2),col='red',lwd=2,lty=2)
xcut <- mean(x2)+(sd(x2[,1]))

# plot histogram of convex hull maximum dimension and set cutoff
hist(main='Convex Hull Max Dimension',x3,breaks=25,col='blue')
abline(v=mean(x3)+(sd(x3[,1])),col='red')
abline(v=mean(x3),col='red',lwd=2,lty=2)
x3cut <- mean(x3)+(sd(x3[,1]))

# create matrix defining whether sites are retained (1) or removed (0)
for (j in 1:nrow(x2_a)) {
	if (x2[j,] < xcut) {
	x2_a[j,] <- 1}
	if (x3[j,] < x3cut) {
	x3_a[j,] <- 1}}
x_rem <- x2_a*x3_a

# create reduced data set
data.red <- matrix(0,nrow(mydata),coln)
for (i in 1:coln) {
data.red[,i] <- mydata[,i]*x_rem}
row.names(data.red) <- row.names(mydata)
colnames(data.red) <- colnames(mydata)
data.red <- data.red[which(rowSums(data.red) >0),]

# recalculate CA axes for reduced data set
mydata2 <- data.red
data.ca2 <- ca(mydata2)
rown <- nrow(mydata2)
coln <- ncol(mydata2)
lab <- colnames(mydata2)


########################################################################
## Define bootstrapped replicates for reduced data set
data.rowsum <- apply(mydata2,1,sum)
data.sim <- rmultinom(nsim,data.rowsum[1],prob=mydata2[1,])
for (i in 2:rown) {
data.sim <- cbind(data.sim,rmultinom(nsim,data.rowsum[i],prob=mydata2[i,]))}

data.sim <- t(data.sim)
data.sim2 <- matrix(rep(0,nsim*rown*coln),nrow=nsim*rown)
for (k in 1:nsim) {
	for (i in 1:rown) {
		data.sim2[(k-1)*rown+i,] <- data.sim[k+(i-1)*nsim,]}}

data.rowsc2 <- data.ca2$rowcoord[,1:2]
data.colsim <- t(t(data.rowsc2) %*% data.sim2[1:rown,]) / apply(data.sim2[1:rown,],2,sum)
for (k in 2:nsim) {
data.colsim <- rbind(data.colsim, t(t(data.rowsc2) %*% data.sim2[((k-1)*rown+1):(k*rown),]) / apply(data.sim2[((k-1)*rown+1):(k*rown),],2,sum))}
data.colsim2 <- matrix(rep(0,nsim*coln*2),nrow=nsim*coln)

data.colsc2 <- data.ca2$colcoord[,1:2]
data.rowsim <- t(t(data.colsc2) %*% t(data.sim2[1:rown,])) / apply(data.sim2[1:rown,],1,sum)
for (k in 2:nsim) {
x <- t(t(data.colsc2) %*% t(data.sim2[((k-1)*rown+1):(k*rown),]))
data.rowsim <- rbind(data.rowsim, x / apply(data.sim2[((k-1)*rown+1):(k*rown),],1,sum))}
data.rowsim2 <- matrix(rep(0,nsim*rown*2),nrow=nsim*rown)

for (j in 1:coln) {
	for (k in 1:nsim) {
		data.colsim2[(j-1)*nsim+k,] <- data.colsim[j+(k-1)*coln,]}}
for (j in 1:rown) {
	for (k in 1:nsim) {
		data.rowsim2[(j-1)*nsim+k,] <- data.rowsim[j+(k-1)*rown,]}}
########################################################################

# plot CA for reduced data set
plot(data.ca2)

# Define CA row convex hulls with peeled perimeters 
#(peeled proportion defined by 'cut' variable)
poly.row2 <- list()
for (j in 1:nrow(mydata2)) {
	points <- data.rowsim2[(nsim*(j-1)+1):(nsim*j),]
	repeat {
		hpts <- chull(points)
		npts <- nrow(points[-hpts,])
		if(npts/nsim<cut) break
		points <- points[-hpts,]}
	hpts <- c(hpts,hpts[1])
	poly.row2[[j]] <- points[hpts,]}

# Define CA column convex hulls with peeled perimeters 
#(peeled proportion defined by 'cut' variable)
poly.col2 <- list()
for (j in 1:coln) {
	points <- data.colsim2[(nsim*(j-1)+1):(nsim*j),]
	repeat {
		hpts <- chull(points)
		npts <- nrow(points[-hpts,])
		if(npts/nsim<cut) break
		points <- points[-hpts,]}
	hpts <- c(hpts,hpts[1])
	poly.col2[[j]] <- points[hpts,]}

# plot bootstrapped data for reduced data set
plot(data.rowsc2[,1],data.rowsc2[,2],type='n',xlab='CA1',ylab='CA2',main='Reduced Data Set')
data.col2 <- t(t(data.rowsc2) %*% mydata2) / apply(mydata2,2,sum) # set column label positions
# use last two digits of color code to set transparency
for (i in 1:length(poly.row2)) { 
polygon(poly.row2[[i]],col='#0000ff10',border='#0000ff10')} 
for (i in 1:length(poly.col2)) { 
polygon(poly.col2[[i]],col='#ff000050',border='#ff000050')
text(data.col2[i,1],data.col2[i,2],lab[i],font=2,cex=0.5)}


## end script