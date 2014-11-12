library(rGeoHash)
base32 = c(0:9,letters[-c(1,9,12,15)])
bbox = toGeohash(c(40.681500,40.844793), c(-74.064627,-73.880828), 8L)


getBB = function(hash, col=rgb(0,0,0,0.5)) {
  out = fromGeohash(paste0(hash,c("00000000","zzzzzzzz")))
  rect(out$lon[1], out$lat[1], out$lon[2], out$lat[2], border=NA, col=col)
}


par(mar=c(0,0,0,0))
plot(0,0,ylim=c(40.61500,40.844793), xlim=c(-74.064627,-73.880828))
osmap()
for (s in base32[16:32]) {
  for(s2 in base32) {
    getBB(paste0("dr5r",s,s2),grey(runif(1,0.2,0.8),0.5))
  }
}
