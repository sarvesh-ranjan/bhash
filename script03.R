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

# Calculate various balanced hashes:
bits = c(1,2,3,4,5,6,7,8,10,15,20)
bh = matrix(NA, nrow=length(weights), ncol=length(bits))
for (i in 1:length(bits)) {
  xvals = 1:(2^bits[i])/(2^bits[i])
  h_inv = approx(weights, geohash, xout=xvals)$y
  h_vals = approx(h_inv, xvals, xout=geohash)$y
  h_vals[is.na(h_vals)] = 0
  balanced_hash = numericToBaseN(h_vals, base=2L, digits=50)
  bh[,i] = balanced_hash
}

# Calculate entropy of g_q and b_q^5, b_q^8, b_q^10, b_q^15
for (k in 1:2) {
  var = c(1L,5L)[k]
  bit_vals = c(1L:15L,18L,20L,25L,30L,35L,40L)
  entropy = matrix(NA, ncol= ncol(bh) + 1L, nrow=length(bit_vals))
  for (i in 1:nrow(entropy)) {
    ghash = apply(z[,1L:bit_vals[i],drop=FALSE],1,paste,collapse="")
    gp = ctapply(d[,var], ghash, sum) / sum(d[,var])
    gp = gp[gp > 0]
    entropy[i,1L] = -1 * sum( gp * log(gp,2) )

    for (j in 1L:ncol(bh)) {
      bhash = substr(bh[,j], 3L, bit_vals[i] + 2L)
      bp = ctapply(d[,var], bhash, sum) / sum(d[,var])
      bp = bp[bp > 0]
      entropy[i,j + 1L] = -1 * sum( bp * log(bp,2) )
    }

    print(entropy[i,])
  }
  saveRDS(entropy, paste0("entropy",k,".Rds"))
}
