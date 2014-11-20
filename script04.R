# calculate over all states and dc
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

output = matrix(NA, ncol=6, nrow=51)
for (j in 1L:51L) {
  print(paste0(j, ":"))
  z = readRDS(paste0("/Users/taylor/files/bhash/output2/bucket_", j, ".Rds"))
  d = z[,1:5]
  z = z[,-(1:5)]

  geohash = z %*% (1 / 2 ^(1:50))
  index = order(geohash)
  z = z[index,]
  d = d[index,]
  geohash = geohash[index,]
  weights = cumsum(d[,1]) / sum(d[,1])

  # Calculate various balanced hashes:
  bits = c(5L,10L,15L)
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
    bit_vals = c(15L)
    entropy = rep(NA, 3L)
    for (m in 1L:length(bits)) {
      bhash = substr(bh[,m], 3L, bit_vals + 2L)
      bp = ctapply(d[,var], bhash, sum) / sum(d[,var])
      bp = bp[bp > 0]
      entropy[m] = -1 * sum( bp * log(bp,2) )
    }
    entropy = entropy / (bit_vals)
    print(entropy)
    output[j,1:3 + (k-1)*3] = entropy
  }
}

rownames(output) = c(state.abb[1:8], "DC", state.abb[9:50])
saveRDS(output, "output_state.Rds")

library(snippets)
output = readRDS("output_state.Rds")

# Annoying, but need to initialize index
j = 3
x = output[,j]
y = output[,j+3L]
txt = rownames(output)
plot(x, y, pch=19, cex=0.5, col="white", main=paste0(c("5","10","15")[j],"-bit"),
    xlim=range(output), ylim=range(output))
grid()
l = add.labels(x, y, txt, 0.3*strwidth(txt), 0.3*strheight(txt), mar=0.5)
index = which(x != l$lx | y != l$ly)

pdf("fig04.pdf", height=4.5, width=12)
par(mfrow=c(1,3))
par(mar=c(2,2,2,2))
for (j in 1:3) {
  x = output[,j]
  y = output[,j+3L]
  txt = rownames(output)
  plot(x, y, pch=19, cex=0.5, col="white", main=paste0(c("5","10","15")[j],"-bit"),
      xlim=range(output), ylim=range(output))
  grid()
  l = add.labels(x, y, txt, 0.3*strwidth(txt), 0.3*strheight(txt), mar=0.5)
  if(j == 3) index = which(x != l$lx | y != l$ly)
  #points(x[index],y[index], pch=19, cex=0.5)
  text(x[-index], y[-index], txt[-index], col=4, cex=0.6)
  #text(l$lx[index], l$ly[index], txt[index], col=4, cex=0.6)
  #segments(x[index], y[index], l$lx[index], l$ly[index] - 0.002, col=3)
}
dev.off()
