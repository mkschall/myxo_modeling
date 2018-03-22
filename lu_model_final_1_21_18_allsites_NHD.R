
#Load packages
#starting with a hierarchical generalized linear model
library(lme4)
library(MuMIn)
library(arm)
library(effects)
library(sjPlot)
library(jagsUI)
library(coda)
library(rjags)
library(doBy)

#Load theYOY myxo data in
yoymyx <- read.table(file="YOY_myxo_data.txt", header=TRUE)

#check prevalence distribution
yoymyxs<- read.table(file="yoymyxsites.txt", header=TRUE)
head(yoymyxs)
hist(yoymyxs$prevalence)
hist(logit(yoymyxs$prevalence))

#Check what we have:
names(yoymyx)
str(yoymyx)
head(yoymyx)

summary(yoymyx)

unique(yoymyx$Sample_id)
unique(yoymyx$stream)
unique(yoymyx$loc)

###ALL SITES
##Going to start with all sites and all myxo locations full model
##Need to separate out a few things first
#Only want data from sites that have been reviewed
yoymyx1<-yoymyx[!(yoymyx$review=="N"),]

head(yoymyx1)
dim(yoymyx1)
unique(yoymyx1$stream)
yoymyx1$review


##now need to separate out only the correct size frame (traditional samples-collected July - August)
yoymyx2<-yoymyx1[(yoymyx1$Type=="t"),]

head(yoymyx2)
dim(yoymyx2)
yoymyx2$Type
###REMOVE pOTOMAC SITES
yoymyx2a<-yoymyx2[!(yoymyx2$stream=="Potomac.River"),]
dim(yoymyx2a)
#INDExING FOR SITE AND YEAR VARIABLES
###Some additional cleanup# for index variables
###Site index
yoymyx2a$site <- as.numeric(as.factor(as.numeric(yoymyx2a$stream)))
length(unique(yoymyx2a$site))
##Year index
yoymyx2a$year.num <- as.numeric(as.factor(as.numeric(yoymyx2a$year)))

#Sort alphabetically
dat<-yoymyx2a[order(yoymyx2a$stream),]
yoymyx2a<-yoymyx2a[order(yoymyx2a$stream),]
head(dat)
unique(dat$stream)
summary(dat)

write.csv(dat, file="yoymyxsum.csv")


###LANDUSE AT THE NHD CATCHMENT LEVEL
####################################################################################
#Read in LU data
NHD<-read.table(file="land_use_catchment.txt", header=TRUE)
head(NHD)
summary(NHD)
head(dat)
summary(dat)
unique(NHD$Stream) 
#FirstRun with very differentLU-conococheague

unique(NHD$Stream)




dat3 <- merge(dat, NHD, by.x = "stream", by.y = "Stream")
head(dat)
summary(dat3)
#dat
unique(dat$stream)
unique(dat3$stream)

###Check out NLCD variables
hist(dat3$Cultivated_NHDp)
hist(car::logit(dat3$Cultivated_NHDp))

hist(dat3$Developed_NHDp)
hist(car::logit(dat3$Developed_NHDp))

hist(dat3$Forest_NHDp)
hist(car::logit(dat3$Forest_NHDp))

##TRANSFORM LU VARIABLES
dat3$logitag_NHD<-car::logit(dat3$Cultivated_NHDp)
dat3$logitdev_NHD<-car::logit(dat3$Developed_NHDp)
dat3$logitfor_NHD<-car::logit(dat3$Forest_NHDp)
########################################################################################
#######################################################################################
##TESTING CORRELATION OF LU VARIABLES
####################################################################################

#NHD level values
cor(cbind(dat3$logitag_NHD,dat3$logitdev_NHD,dat3$logitfor_NHD))


###values at the local and NHD catchment scale are correlated - check to see if NHD catchment scale yields the
#same results as at the local catchment. If so maybe don't need both
###ACCUMULATED CATCHMENT##########################################################################
####################################################################################################
###Let's look at accumulated catchment landuse
##Read in accumulated catchment data
accum<-read.table(file="land_use_accum.txt", header=TRUE)
head(accum)
head(dat3)

dat4 <- merge(dat3, accum, by.x = "stream", by.y = "stream")

head(dat4)
summary(dat4)
dat4
unique(dat4$stream)
length(unique(dat$stream))

min(dat4$Cultivated_accump)
max(dat4$Cultivated_accump)
mean(dat4$Cultivated_accump)


min(dat4$Developed_accump)
max(dat4$Developed_accump)
mean(dat4$Developed_accump)

###Check out NLCD variables
hist(dat4$Cultivated_accump)
hist(car::logit(dat4$Cultivated_accump))

hist(dat4$Developed_accump)
hist(car::logit(dat4$Developed_accump))

hist(dat4$Forest_accump)
hist(car::logit(dat4$Forest_accump))

#Transform Accumulated Catchment variables
dat4$logitag_accum<-car::logit(dat4$Cultivated_accump)
dat4$logitdev_accum<-car::logit(dat4$Developed_accump)
dat4$logitfor_accum<-car::logit(dat4$Forest_accump)

#NHD and Catchment values
cor(cbind(dat4$logitag_NHD,dat4$logitdev_NHD,dat4$logitfor_NHD,dat4$logitag_accum, dat4$logitdev_accum, dat4$logitfor_accum))
####################################################################

#GLMER MODEL
#with NHD landuse
M4e<-glmer(m~1+ scale(logitdev_accum) + scale(logitag_accum) + scale(logitdev_NHD)+scale(logitag_NHD)+  (1|stream) + (1|year), family=binomial, data=dat4)
summary(M4e)
##############################################################################################

#############################################################################################################
##JAGS CODE SIMILAR TO ABOVE BUT FOR ACCUMULATED CATCHMENT and NHD - CAUTION AGAIN WITH MULTICOLLINEARITY
#JAGS CODE FOR AG at ACCUMULATED catchment scale#####
##############################################
##############################################################################
agNHD <- as.numeric(by(dat4$logitag_NHD, dat4$site, mean))
agNHD<-as.numeric(scale(agNHD))

devNHD<-as.numeric(by(dat4$logitdev_NHD, dat4$site, mean))
devNHD<-as.numeric(scale(devNHD))

agaccum <- as.numeric(by(dat4$logitag_accum, dat4$site, mean))
agaccum<-as.numeric(scale(agaccum))

devaccum<-as.numeric(by(dat4$logitdev_accum, dat4$site, mean))
devaccum<-as.numeric(scale(devaccum))



###Model statement
sink("model6.txt")
cat("
    model{
    for(i in 1:n){
    
    y[i]~dbin(p.bound[i],1)
    p.bound[i]<-max(0, min(1, p[i]))
    logit(p[i])<-Xbeta[i]
    Xbeta[i]<- b.site[site[i]] + b.year[year[i]]
    
    }
    
    
    for(j in 1:n.site){
    b.site[j]~dnorm(mu.hat[j],tau.site)
    mu.hat[j] <- b.0 + b[1] * agNHD[j] + b[2]*devNHD[j] + b[3]*agaccum[j] +b[4]*devaccum[j]
    }
    
    b.0~dnorm(0,0.001)
    # Bayesian LASSO -  a Laplace (double exponential) prior
    for(i in 1:4){
    b[i] ~ ddexp(0, lambda1)
    }
    
    lambda1 ~ dexp(10)
    
    
    for(j in 1:n.year){
    b.year[j]~dnorm(0, tau.year)
   
    }
    

    tau.site<-pow(sigma.site,-2)
    tau.year<-pow(sigma.year,-2)
    
    sigma.site~dunif(0,100)
    sigma.year~dunif(0,100)
    }
    ",fill = TRUE)
sink()

#########Let's Run in JAGS
##RATIONALE FOR JAGS - want to be able to give probability along with relationship between myxo prevalance
#and different predictors
## Bundle data
data <- list(y=dat4$m, site=dat4$site,agNHD=agNHD, devNHD=devNHD,agaccum=agaccum,devaccum=devaccum,year=dat4$year.num,
             n = nrow(dat4), n.site=length(unique(dat4$stream)), 
             n.year=length(unique(dat4$year)) )

# Create starting values 
inits <- function(){list(b.0=rnorm(1),
                         sigma.site=runif(1), sigma.year=runif(1)
                         
)}

# Parameters to estimate/track
parameters <- c("b.0","b.site","b","b.year","sigma.year", "sigma.site")

# MCMC setting
n.iter <- 90000
n.thin <- 3
n.burn <- 40000
n.chain <- 3
########## JAGS ANALYSIS FOR AG ###########

outall <- jags(data = data, inits = inits, parameters.to.save = parameters, 
                 model.file = "model6.txt", n.chains = n.chain, n.thin = n.thin, n.iter = n.iter, 
                 n.burnin = n.burn, parallel=T)


# Summarize the result
print(outall, digits = 3)
#summary(outall)
#xyplot(outall)
#library(R2jags)


#Probability that the trend with ag posNHD
pagNHD<- mean(outall$sims.list$b[,1] > 0)
pagNHD

pdevNHD <-mean(outall$sims.list$b[,2]<0)
pdevNHD

#Probability that the trend with NHD ag is positive
pagaccum <- mean(outall$sims.list$b[,3] > 0)
pagaccum
#Probability that the trend with NHD development is negative
pdevaccum <- mean(outall$sims.list$b[,4] < 0)
pdevaccum

###Plotting relationship with LU predictors - Just going to do local ag and devleoped
###################################################################################
# Select random intercepts for sites
mean.beta <- outall$mean$b.site

# data to predict
fake1 <- seq(min(agNHD), max(agNHD), length=length(agNHD))

# Obtain parameters of interest-coefficients
#site coefficients
ests<-as.numeric(outall$mean$b.0)
#slope predictor
ests2<-outall$mean$b[1]


# Fitted lines
fit1 <- ests + ests2*fake1

# Obtain 95% CIs for fitted line
est.lineA <- matrix(NA, ncol=length(fake1), nrow=length(outall$sims.list$b.0)) #container for predicted values


for(i in 1:length(outall$sims.list$b.0)){
  for(t in 1:length(fake1)){
    est.lineA[i,t] <- outall$sims.list$b.0[i] + outall$sims.list$b[,1][i] * fake1[t] 
  }
}

# CIs for fitted values
upper.CIA <- apply(est.lineA, 2, quantile, probs=c(0.975) )
lower.CIA <- apply(est.lineA, 2, quantile, probs=c(0.025) )

## Grab 95% CIs for beta's
  u.alpha<- outall$q97.5$b.site 
  l.alpha<- outall$q2.5$b.site


###########################################
####### FIGURE WITH CRI's
res <- 6
name_figure <- "site_NHDAg_up.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...
par(mfrow = c(1,1), mar=c(3,3,1,1)) 


size.labels = 1
size.text = 1

x.label <- 'Logit (proportion NHD agriculture)'
y.label <- 'Logit (myxozoan prevalance)'


plot(mean.beta ~ agNHD,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c(lower.CIA[i.for] , upper.CIA[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(agNHD, mean.beta,pch=16,cex=0.8)

segments(x0=agNHD, x1=agNHD,
         y0=l.alpha, y1=u.alpha, col='black',lwd=1)


lines(fake1,fit1, lwd = 3, col="black", lty = 1)

mtext(x.label, line = 1.5, side = 1, cex = size.text)
mtext(y.label, line = 1.7, side = 2, cex = size.text)

box()
par(def.par)
dev.off()

####### END
##############################################
###SAME FIGURE BUT FOR NHD DEVELOPMENT
#Select random intercepts for sites
mean.beta <- outall$mean$b.site

# data to predict
fake1 <- seq(min(devNHD), max(devNHD), length=length(devNHD))

# Obtain parameters of interest-coefficients
#site coefficients
ests<-as.numeric(outall$mean$b.0)
#slope predictor
ests2<-outall$mean$b[2]


# Fitted lines
fit1 <- ests + ests2*fake1

# Obtain 95% CIs for fitted line
est.lineB <- matrix(NA, ncol=length(fake1), nrow=length(outall$sims.list$b.0)) #container for predicted values


for(i in 1:length(outall$sims.list$b.0)){
  for(t in 1:length(fake1)){
    est.lineB[i,t] <- outall$sims.list$b.0[i] + outall$sims.list$b[,2][i] * fake1[t] 
  }
}

# CIs for fitted values
upper.CIA <- apply(est.lineB, 2, quantile, probs=c(0.975) )
lower.CIA <- apply(est.lineB, 2, quantile, probs=c(0.025) )

## Grab 95% CIs for beta's
u.alpha<- outall$q97.5$b.site 
l.alpha<- outall$q2.5$b.site


###########################################
####### FIGURE WITH CRI's
res <- 6
name_figure <- "site_NHDdev.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...
par(mfrow = c(1,1), mar=c(3,3,1,1)) 


size.labels = 1
size.text = 1

x.label <- 'Logit (proportion NHD development)'
y.label <- 'Logit (myxozoan prevalance)'


plot(mean.beta ~ devNHD,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c(lower.CIA[i.for] , upper.CIA[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(devNHD, mean.beta,pch=16,cex=0.8)

segments(x0=devNHD, x1=devNHD,
         y0=l.alpha, y1=u.alpha, col='black',lwd=1)


lines(fake1,fit1, lwd = 3, col="black", lty = 1)

mtext(x.label, line = 1.5, side = 1, cex = size.text)
mtext(y.label, line = 1.7, side = 2, cex = size.text)

box()
par(def.par)
dev.off()

####### END
##########################################################
######Plotting site and year variability
# site.ran <- matrix(NA, nrow=out$mcmc.info$n.samples, ncol=26)
site.ran <- plogis(outall$sims.list$b.site)
site.means <- apply(site.ran, 2, mean)
site.quants <- apply(site.ran, 2, quantile, c(0.025, 0.975))

#year ran effects
year.ran<-plogis(outall$sims.list$b.0 + outall$sims.list$b.year)
year.means<-apply(year.ran, 2, mean)
year.quants<-apply(year.ran, 2,quantile,c(0.025, 0.975))

##overall intercept
int<-plogis(outall$sims.list$b.0)
mean(int)
quantile(int,0.025)
quantile(int,0.975)


######################################################################
#####################################################################
##################PLOTTING#####################################
##Plotting the effects
##### YEAR EFFECT ####
#Use probabilities
#use means and 95% credible intervals summarized above
year.means<-apply(year.ran, 2, mean)
year.quants<-apply(year.ran, 2,quantile,c(0.025, 0.975))
lower.quant<-year.quants[1,]
upper.quant<-year.quants[2,]
year2<-c(2013,2014,2015,2016)

# Plot effect
#YEAR 
res <- 6
name_figure <- "YEAR_effect.png"
png(filename = name_figure, height = 300*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 1
x.label = 'Year'
y.label = expression(paste('Myxozoan prevalence'))

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(1,1,0.5,1), oma=c(2,2.5,0.1,0.1) )

plot(year.means~year2, ylim=c(0.0,0.9),
     xlim=c(min(dat4$year),max(dat4$year)), axes=F, ylab='', xlab='', type='n')
axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0), at=year2)
axis(side=2,cex.axis=size.text, las=1,mgp=c(0,0.5,0),tck=-0.01)

points(year2, year.means, cex=2, col='black', pch=16)
segments (year2,lower.quant,year2,upper.quant, lwd=1, col="gray10")

mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T)
mtext(y.label, line = 1.2, side = 2, cex = size.text, outer=T)

box()
par(def.par)
dev.off()

###################################################################################
##################################################################################
##SITE VARIABILITY
##Site over the duration
##Plotting the effects
##### Site Variability ####
#Use probabilities
site.means<-apply(site.ran, 2, mean)
summary(site.means)
site.quants<-apply(site.ran, 2,quantile,c(0.025, 0.975))
lower.quant<-site.quants[1,]
upper.quant<-site.quants[2,]
site2<-seq(1,31, by=1)
#Need to reorder to match table
site.means1<-site.means[c(1,10,18,17,6,12,16)]
mean(site.means1[1:4])

site.means2<-site.means[c(2:5,9,11,14,15,27:30,19,20,31,7,13,21,22,23,24,26,8,25)]
mean(site.means2)

site.meanspl<-c(site.means2,site.means1)

lower.quant1<-lower.quant[c(1,10,18,17,6,12,16)]
lower.quant2<-lower.quant[c(2:5,9,11,14,15,27:30,19,20,31,7,13,21,22,23,24,26,8,25)]
mean(lower.quant2)

lower.quantpl<-c(lower.quant2, lower.quant1)

upper.quant1<-upper.quant[c(1,10,18,17,6,12,16)]
upper.quant2<-upper.quant[c(2:5,9,11,14,15,27:30,19,20,31,7,13,21,22,23,24,26,8,25)]
mean(upper.quant2)

upper.quantpl<-c(upper.quant2, upper.quant1)
str_names<-sort(unique(dat$stream))
str_names1<-str_names[c(1,10,18,17,6,12,16)]
str_names2<-str_names[c(2:5,9,11,14,15,27:30,19,20,31,7,13,21,22,23,24,26,8,25)]
str_namespl<-c(as.character(str_names2), as.character(str_names1))
labels<-str_namespl
labels[28]<-'Schuylkill.River.Port.Clinton'
labels[19]<-'Susquehanna.River.Isle.Of.Ques'
labels <- gsub("\\.", " ",labels )

# Plot effect
#Site EFFECT
res <- 6
name_figure <- "Site.png"
png(filename = name_figure, height = 600*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1.5
size.text = 1.5
axissize <- 1
x.label = 'Site'
y.label = 'Myxozoan prevalence'

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(14,1,0.5,1), oma=c(3,3,0.1,0.1) )

plot(site.meanspl~site2, ylim=c(0,1),
     xlim=c(0,32), axes=F, ylab='', xlab='', type='n')
axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0), at=site2, labels = FALSE)
axis(side=2,cex.axis=size.text, las=1,mgp=c(0,0.25,0),tck=-0.01)

text(seq(1,31,by=1), par("usr")[1] + 1.2,labels=labels, srt=75, pos=2, xpd=TRUE, cex=1)
points(site2[1:24], site.meanspl[1:24], cex=2, col='black', pch=16)
points(site2[25:28], site.meanspl[25:28], cex=2,col='gray50', pch=15)
points(site2[29:31], site.meanspl[29:31], cex=2, col='gray70', pch=17)
segments (site2[1:24],lower.quantpl[1:24],site2[1:24],upper.quantpl[1:24], lwd=1, col="black")
segments (site2[25:28],lower.quantpl[25:28],site2[25:28],upper.quantpl[25:28], lwd=1, col="gray50")
segments (site2[29:31],lower.quantpl[29:31],site2[29:31],upper.quantpl[29:31], lwd=1, col="gray70")
#abline(v=24.5, lty=2, lwd=1)
#abline(v=27.5,lty=2,lwd=1)

mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T)
mtext(y.label, line = 1.1, side = 2,adj=0.75, cex = size.text, outer=T)



box()
par(def.par)
dev.off()
#########################################################################################################
