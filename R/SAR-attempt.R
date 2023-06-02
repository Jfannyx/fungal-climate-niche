rm(list=ls())
gc()
neon_dob <- readRDS("phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))

neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
neon <- subset_samples(neon, !is.na(horizon))
rm(neon_dob)
d <- sample_data(neon)
a <- paste(d$siteID,d$sampleTiming,d$horizon,sep=",")
sort(table(a))
apply(d, 2, function(x) sum(is.na(x)))
aggregate(d,list(d$horizon), function(x) sum(is.na(x)))
a <- aggregate(cbind(d$lon,d$geneticSampleID,d$sampleTiming,d$collectYear, d$horizon),list(d$siteID), function(x) length(unique(x)))
row.names(a) <- a[,1]
a <- a[,-1]
colnames(a) <- c("longitute","geneticSample","samplingTime","collectYear","horizon")
a <- t(a)
a
table(d$collectYear)
cbind(aggregate(d$collectYear,list(d$siteID),unique),
      aggregate(d$collectYear,list(d$siteID),function(x) length(unique(x)))$x)

a <- unique(d$Site)
result.c <- matrix(nrow = 45, ncol = 4)
result.z <- matrix(nrow = 45, ncol = 4)
for (i in 1:length(a)){
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub))
  shannon <- vector(length = dim1[1])
  for (j in 1:(dim1[1]-1)){
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(neon_sub, flag)
    shannon[j] <- vegan::diversity(otu_table(temp), index = "shannon")["TRUE"]
  }
  ex <- as.data.frame(cbind(shannon, "A"=c(1:dim1[1])))
  temp <- summary(nls(shannon~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  result.c[i,] <- temp[1,]
  result.z[i,] <- temp[2,]
}
a <- aggregate(d,list(d$site),function(x) mean(x,na.rm = T))

#dob
dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
rm(neon_dob)

d1 <- sample_data(dob)
a <- paste(d1$siteID,d$horizon,sep=",")
sort(table(a))
table(d1$horizon)
apply(d, 2, function(x) sum(is.na(x)))
aggregate(d,list(d$horizon), function(x) sum(is.na(x)))
a <- aggregate(cbind(d$lon,d$geneticSampleID,d$sampleTiming,d$collectYear, d$horizon),list(d$siteID), function(x) length(unique(x)))
row.names(a) <- a[,1]
a <- a[,-1]
colnames(a) <- c("longitute","geneticSample","samplingTime","collectYear","horizon")
a <- t(a)
a
table(d$collectYear)
cbind(aggregate(d$collectYear,list(d$siteID),unique),
      aggregate(d$collectYear,list(d$siteID),function(x) length(unique(x)))$x)

a <- unique(d1$Site)
result.c <- matrix(nrow = 68, ncol = 4)
result.z <- matrix(nrow = 68, ncol = 4)
for (i in 1:length(a)){
  dob_sub <- subset_samples(dob, Site==a[i])
  dim1 <- dim(otu_table(dob_sub))
  shannon <- vector(length = dim1[1])
  for (j in 1:(dim1[1]-1)){
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(dob_sub, flag)
    shannon[j] <- vegan::diversity(otu_table(temp), index = "shannon")["TRUE"]
  }
  ex <- as.data.frame(cbind(shannon, "A"=c(1:dim1[1])))
  temp <- summary(nls(shannon~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  result.c[i,] <- temp[1,]
  result.z[i,] <- temp[2,]
}
colnames(result.z) <- colnames(summary(nls(shannon~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
