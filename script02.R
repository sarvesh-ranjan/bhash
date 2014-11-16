library(baseN)
library(rGeoHash)
library(iotools)

toLatLon = function(w) {
 lon = w[,seq(1,49,by=2),drop=FALSE] %*% (1 / 2 ^(1:25)) * 360 - 180
 lat = w[,seq(2,50,by=2),drop=FALSE] %*% (1 / 2 ^(1:25)) * 180 - 90
 return(list(lat=lat, lon=lon))
}

getBB = function(hash, col=rgb(0,0,0,0.5)) {
  out = fromGeohash(paste0(hash,c("000000")))
  out2 = fromGeohash(paste0(hash,c("zzzzzz")))
  rect(out$lon, out$lat, out2$lon, out2$lat, border=NA, col=col)
}

fnames = dir("../UScensus2010blk/data/", pattern="*.rda")
fnames = gsub(".blk10.rda", "", fnames, fixed=TRUE)

z = readRDS("output2/bucket_1.Rds")
d = z[,1:5]
z = z[,-(1:5)]

geohash = z %*% (1 / 2 ^(1:50))
index = order(geohash)
z = z[index,]
d = d[index,]
geohash = geohash[index,]
weights = cumsum(d[,1]) / sum(d[,1])

# Calculate the balanced hash
bits = 10
xvals = 1:(2^bits)/(2^bits)
h_inv = approx(weights, geohash, xout=xvals)$y
h_vals = approx(h_inv, xvals, xout=geohash)$y
h_vals[is.na(h_vals)] = 0
balanced_hash = numericToBaseN(h_vals, base=2L, digits=50)

# Entropy calculation
entropy = matrix(NA, ncol=10, nrow=20L)
for (i in 1:nrow(entropy)) {
  ghash = apply(z[,1:i,drop=FALSE],1,paste,collapse="")
  bhash = substr(balanced_hash, 3L, i + 2L)
  for (j in 1:5) {
    gp = ctapply(d[,j], ghash, sum) / sum(d[,j])
    bp = ctapply(d[,j], bhash, sum) / sum(d[,j])
    gp = gp[gp > 0]
    bp = bp[bp > 0]
    entropy[i,j + 0L] = -1 * sum( gp * log(gp,2) )
    entropy[i,j + 5L] = -1 * sum( bp * log(bp,2) )
  }
  print(entropy[i,])
}

# Plot entropy
pdf("fig01.pdf", height=6, width=12)
  par(mfrow=c(1,2))
  plot(entropy[,1], type="l", lwd=2, col="olivedrab", ylim=c(0,40),
        ylab="Shannon entropy",xlab="hash bits")
  lines(entropy[,6], lwd=2, col="salmon")
  lines(entropy[,10], lwd=2, col="darkblue")
  legend(5,35,c("theoretical limit", "balanced hash", "geohash", "balanced hash (renters)"),
      col=c("black","salmon", "olivedrab", "darkblue"), lwd=2, lty=c("dashed", "solid", "solid"))
  abline(0,1,lty="dashed")
  plot(diff(entropy), type="l", lwd=2, col="olivedrab", ylim=c(0,1),
        ylab=expression(paste(Delta, " in Shannon entropy")),xlab="hash bits")
  lines(diff(bentropy), lwd=2, col="salmon")
  abline(h=1,lty="dashed")
dev.off()


# Map balanced_hash break points to geohash
pdf("fig02.pdf", height=12, width=12)
par(mar=c(0,0,0,0))
par(mfrow=c(2,2))
for (bits in c(3,4,5,6)) {
  hash = substr(balanced_hash, 3L, bits + 2L)
  index = setdiff(which(hash[-1] != hash[-length(hash)]),length(hash)-1)
  index = c(1, index, length(hash))
  w = toLatLon(z[index,-1])

  plot(w$lon, w$lat, axes=FALSE)
  box()
  osmap()
  for(i in 1:(length(index)-1)) {
    cl = rainbow(8L, alpha=0.5)[(i %% 8) + 1]
    w2 = toLatLon(z[index[i]:(index[i+1]-1),-1,drop=FALSE])
    cm = lapply(w2, median)
    w2 = toGeohash(w2$lat, w2$lon, 5L)
    getBB(unique(w2), cl)
    text(cm$lon, cm$lat, sprintf("%02d", i))
  }
}
dev.off()

# Cycle over all the states:
fnames = dir("../UScensus2010blk/data/", pattern="*.rda")
fnames = gsub(".blk10.rda", "", fnames, fixed=TRUE)

output = NULL
for (j in 1:length(fnames)) {
  print(paste0(fnames[j], ":"))
  z = readRDS(paste0("/Users/taylor/files/bhash/output2/bucket_", j, ".Rds"))
  d = z[,1:5]
  z = z[,-(1:5)]

  geohash = z %*% (1 / 2 ^(1:50))
  index = order(geohash)
  z = z[index,]
  d = d[index,]
  geohash = geohash[index,]
  weights = cumsum(d[,1]) / sum(d[,1])

  # Calculate the balanced hash
  bits = 10
  xvals = 1:(2^bits)/(2^bits)
  h_inv = approx(weights, geohash, xout=xvals)$y
  h_vals = approx(h_inv, xvals, xout=geohash)$y
  h_vals[is.na(h_vals)] = 0
  balanced_hash = numericToBaseN(h_vals, base=2L, digits=50)

  # Entropy calculation
  bit_vals = c(1L,5L,10L,20L,30L)
  var_vals = c(1L,5L)
  entropy = matrix(NA, ncol=2 * length(var_vals), nrow=length(bit_vals))
  for (i in 1:nrow(entropy)) {
    ghash = apply(z[,1:bit_vals[i],drop=FALSE],1,paste,collapse="")
    bhash = substr(balanced_hash, 3L, bit_vals[i] + 2L)
    for (j in 1:length(var_vals)) {
      gp = ctapply(d[,var_vals[j]], ghash, sum) / sum(d[,var_vals[j]])
      bp = ctapply(d[,var_vals[j]], bhash, sum) / sum(d[,var_vals[j]])
      gp = gp[gp > 0]
      bp = bp[bp > 0]
      entropy[i,j] = -1 * sum( gp * log(gp,2) )
      entropy[i,j + length(var_vals)] = -1 * sum( bp * log(bp,2) )
    }
    print(entropy[i,])
  }

  output = rbind(output, as.numeric(entropy))
}

write.csv(output, "entropy.csv")

#####
library(xtable)
fnames = dir("../UScensus2010blk/data/", pattern="*.rda")
fnames = gsub(".blk10.rda", "", fnames, fixed=TRUE)
output = read.csv("entropy.csv")
row.names(output) = fnames
xtable(output[,-c(1,2,6,7)])