library(baseN)
library(rGeoHash)
library(iotools)

z = readRDS("output/bucket_36.Rds")

geohash = z[,2:51] %*% (1 / 2 ^(1:50))
index = order(geohash)
z = z[index,]
geohash = geohash[index,]
weights = cumsum(z[,1]) / sum(z[,1])

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


# Entropy
entropy = rep(1,40L)
for (i in 1:length(entropy)) {
  hash = apply(z[,2:(i+1),drop=FALSE],1,paste,collapse="")
  p = ctapply(z[,1], hash, sum) / sum(z[,1])
  entropy[i] = -1 * sum( p * log(p,2) )
  print(entropy[i])
}

# Balanced hash
bits = 10
xvals = 1:(2^bits)/(2^bits)
h_inv = approx(weights, geohash, xout=xvals)$y
h_vals = approx(h_inv, xvals, xout=geohash)$y
h_vals[is.na(h_vals)] = 0
bhash = numericToBaseN(h_vals, base=2L, digits=50)

# Entropy
bentropy = rep(1,40L)
for (i in 1:length(entropy)) {
  hash = substr(bhash, 3L, i + 2L)
  p = ctapply(z[,1], hash, sum) / sum(z[,1])
  bentropy[i] = -1 * sum( p * log(p,2) )
  print(bentropy[i])
}

# Plot entropy
pdf("fig01.pdf", height=6, width=12)
par(mfrow=c(1,2))
plot(entropy, type="l", lwd=2, col="olivedrab", ylim=c(0,40),
      ylab="Shannon entropy",xlab="hash bits")
lines(bentropy, lwd=2, col="salmon")
legend(5,35,c("theoretical", "balanced hash", "geohash"),
    col=c("black","salmon", "olivedrab"), lwd=2, lty=c("dashed", "solid", "solid"))
abline(0,1,lty="dashed")
plot(diff(entropy), type="l", lwd=2, col="olivedrab", ylim=c(0,1),
      ylab=expression(paste(Delta, " in Shannon entropy")),xlab="hash bits")
lines(diff(bentropy), lwd=2, col="salmon")
abline(h=1,lty="dashed")
dev.off()



# Map bhash break points to geohash
pdf("fig02.pdf", height=12, width=12)
par(mar=c(0,0,0,0))
par(mfrow=c(2,2))
for (bits in c(3,4,5,6)) {
  hash = substr(bhash, 3L, bits + 2L)
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


