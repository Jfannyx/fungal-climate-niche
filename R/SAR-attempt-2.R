rm(list=ls())
gc()

# loading and filtering
neon_dob <- readRDS("phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
neon <- subset_samples(neon, !is.na(horizon))
#rm(neon_dob)

# neon dataset
d <- sample_data(neon) # sample data data frame

# some exploratory data analysis
# a <- paste(d$siteID,d$sampleTiming,d$horizon,sep=",")
# sort(table(a))
# apply(d, 2, function(x) sum(is.na(x)))
# aggregate(d,list(d$horizon), function(x) sum(is.na(x)))
# a <- aggregate(cbind(d$lon,d$geneticSampleID,d$sampleTiming,d$collectYear, d$horizon),list(d$siteID), function(x) length(unique(x)))
# row.names(a) <- a[,1]
# a <- a[,-1]
# colnames(a) <- c("longitute","geneticSample","samplingTime","collectYear","horizon")
# a <- t(a)
# a
# table(d$collectYear)
# cbind(aggregate(d$collectYear,list(d$siteID),unique),
#       aggregate(d$collectYear,list(d$siteID),function(x) length(unique(x)))$x)
# a <- aggregate(d,list(d$site),function(x) mean(x,na.rm = T))

a <- unique(d$Site) # get all the sites
# create matrix to save result
result.c <- matrix(nrow = 45, ncol = 4)
result.z <- matrix(nrow = 45, ncol = 4)

# computing the relationship between shannon diversity and number of samples
for (i in 1:length(a)){
  # take out one site
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  for (j in 1:(dim1[1])){ 
    # error occurs if compute the diversity of all samples in one site
    # so only take out 1 : (number of samples - 1) to compute 
    # will fix later
    
    # randomly sample j samples in the site 
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(neon_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
    
    # compute the shannon diversity
    # shannon[j] <- vegan::diversity(otu_table(temp), index = "shannon")["TRUE"]
    species[j] <- sum(otu_table(temp)["TRUE"] > 0)
  }
  ex <- as.data.frame(cbind(species, "A"=c(1:dim1[1])))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  result.c[i,] <- temp[1,]
  result.z[i,] <- temp[2,]
}
colnames(result.z) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
colnames(result.c) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
result.z
result.c

neon_site <- merge_samples(neon, "Site")
d_site <- sample_data(neon_site)
summary(lm(result.z ~ d_site$soilInCaClpH + d_site$soilMoisture + d_site$mat_celsius +
     d_site$map_mm + d_site$temp_seasonality))

par(mfrow = c(2,2))
plot(y = result.z[,1], x = d_site$soilInCaClpH)
abline(reg = lm(result.z[,1]~d_site$soilInCaClpH), col = 2)
plot(y = result.z[,1], x = d_site$soilMoisture)
abline(reg = lm(result.z[,1]~d_site$soilMoisture), col = 2)
plot(y = result.z[,1], x = d_site$mat_celsius)
abline(reg = lm(result.z[,1]~d_site$mat_celsius), col = 2)
plot(y = result.z[,1], x = d_site$map_mm)
abline(reg = lm(result.z[,1]~d_site$map_mm), col = 2)

library(ggtrendline)
library(gridExtra)
p1 <- ggtrendline(d_site$soilInCaClpH, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilInCaClpH, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilInCaClpH")+ylab("z")

p2 <- ggtrendline(d_site$soilMoisture, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilMoisture, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$mat_celsius, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$mat_celsius, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z")

p4 <- ggtrendline(d_site$map_mm, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$map_mm, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

p5 <- ggtrendline(d_site$temp_seasonality, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$temp_seasonality, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("temp_seasonality")+ylab("z") 

grid.arrange(p2,p3,p4,p5,ncol = 2)

#dob
dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
rm(neon_dob)

d1 <- sample_data(dob)
# a <- paste(d1$siteID,d$horizon,sep=",")
# sort(table(a))
# table(d1$horizon)
# apply(d, 2, function(x) sum(is.na(x)))
# aggregate(d,list(d$horizon), function(x) sum(is.na(x)))
# a <- aggregate(cbind(d$lon,d$geneticSampleID,d$sampleTiming,d$collectYear, d$horizon),list(d$siteID), function(x) length(unique(x)))
# row.names(a) <- a[,1]
# a <- a[,-1]
# colnames(a) <- c("longitute","geneticSample","samplingTime","collectYear","horizon")
# a <- t(a)
# a
# table(d$collectYear)
# cbind(aggregate(d$collectYear,list(d$siteID),unique),
#       aggregate(d$collectYear,list(d$siteID),function(x) length(unique(x)))$x)

a <- unique(d1$Site)
result.c <- matrix(nrow = 68, ncol = 4)
result.z <- matrix(nrow = 68, ncol = 4)
for (i in 1:length(a)){
  # take out one site
  dob_sub <- subset_samples(dob, Site==a[i])
  dim1 <- dim(otu_table(dob_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  for (j in 1:(dim1[1])){ 
    # error occurs if compute the diversity of all samples in one site
    # so only take out 1 : (number of samples - 1) to compute 
    # will fix later
    
    # randomly sample j samples in the site 
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(dob_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
    
    # compute the shannon diversity
    # shannon[j] <- vegan::diversity(otu_table(temp), index = "shannon")["TRUE"]
    species[j] <- sum(otu_table(temp)["TRUE"] > 0)
  }
  ex <- as.data.frame(cbind(species, "A"=c(1:dim1[1])))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  result.c[i,] <- temp[1,]
  result.z[i,] <- temp[2,]
}
colnames(result.z) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
result.z[,4] > 0.05
hist(result.z[,1])
summary(result.z[,1])

dob_site <- merge_samples(dob, "Site")
d_site <- sample_data(dob_site)
summary(lm(result.z[,1] ~ d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality))

p1 <- ggtrendline(d_site$soilInCaClpH, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilInCaClpH, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilInCaClpH")+ylab("z")

p2 <- ggtrendline(d_site$soilMoisture, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilMoisture, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$mat_celsius, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$mat_celsius, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z")

p4 <- ggtrendline(d_site$map_mm, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$map_mm, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

p5 <- ggtrendline(d_site$temp_seasonality, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$temp_seasonality, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("temp_seasonality")+ylab("z") 

grid.arrange(p2,p3,p4,p5,ncol = 2)
