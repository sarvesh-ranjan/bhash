load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}

toBitGeoHash = function (lat, lon, bits)
{
  lat_bits = floor(bits/2)
  lon_bits = ceiling(bits/2)
  lat = (lat + 90)/180
  lon = (lon + 180)/360
  lat_val = numericToBaseN(lat, digits = lat_bits)
  lat_val = matrix(unlist(strsplit(lat_val, "")), nrow = lat_bits +
      2)[-(1:2), , drop = FALSE]
  lon_val = numericToBaseN(lon, digits = lon_bits)
  lon_val = matrix(unlist(strsplit(lon_val, "")), nrow = lon_bits +
      2)[-(1:2), , drop = FALSE]
  if ((tot_bits <- lat_bits + lon_bits)%%2 == 1L)
      tot_bits = tot_bits + 1
  ind = as.numeric(t(matrix(1:(tot_bits), ncol = 2)))
  mat = t(rbind(lon_val, lat_val)[ind[1L:(lat_bits + lon_bits)],
      , drop = FALSE])
  return(mat)
}

toLatLon = function(w) {
 lon = w[,seq(1,49,by=2)] %*% (1 / 2 ^(1:25)) * 360 - 180
 lat = w[,seq(2,50,by=2)] %*% (1 / 2 ^(1:25)) * 180 - 90
 return(list(lat=lat, lon=lon))
}

library(baseN)
library(rGeoHash)
library(sp)
fnames = dir("UScensus2010blk/data/", pattern="*.rda")

for(j in 1:length(fnames)) {
  z = load_obj(paste0("UScensus2010blk/data/", fnames[j]))

  mat = matrix(unlist(lapply(z@polygons, function(v) v@labpt)), ncol=2, byrow=TRUE,
                dimnames=list(1:nrow(z)))
  hash = toBitGeoHash(mat[,2], mat[,1], 10L)
  hash = cbind(as.integer(z@data[,c("P0010001")]), apply(hash, 2, as.integer))

  saveRDS(hash[hash[,1] > 0,], paste0("~/Desktop/output/bucket_", j, ".Rds"))
}
################
library(baseN)
library(iotools)
fnames = dir(".", pattern="*[0-9].Rds", full.names=TRUE)
z = NULL
for (j in 1:length(fnames)) {
  z = rbind(z, readRDS(fnames[j]))
}
hash_string = apply(z[,-1], 1, paste, collapse="")
index = order(hash_string)
z = z[index,]
weights = cumsum(z[,1]) / sum(z[,1])

bits = ceiling(-1 * log(min(diff(weights)),2))
bhash = numericToBaseN(weights, digits = bits)

saveRDS(bhash, "/Users/taylor/Desktop/output/bucket_allData_bhash.Rds")

###############
z = readRDS("bucket_38.Rds")
hash_string = apply(z[,-1], 1, paste, collapse="")
index = order(hash_string)
z = z[index,]

weights = cumsum(z[,1]) / sum(z[,1])
y = z[,2:51] %*% (1 / 2 ^(1:50))

plot(y,weights, type="l")
abline(lm(weights ~ y))

################
library(baseN)
library(iotools)
fnames = dir("output/", pattern="*[0-9].Rds", full.names=TRUE)
z = NULL
for (j in 1:length(fnames)) {
  temp = readRDS(fnames[j])
  index = sample(1:nrow(temp), floor(nrow(temp)*0.05))
  z = rbind(z, temp[index,])
}
hash_string = apply(z[,-1], 1, paste, collapse="")
index = order(hash_string)
z = z[index,]

weights = cumsum(z[,1]) / sum(z[,1])
y = z[,2:51] %*% (1 / 2 ^(1:50))

N = 5*5
out = rep(NA, N-1)
for (i in 2:N) {
  p = i
  xvals = 1:(2^p)/(2^p)
  h_inv = approx(weights, y, xout=xvals)$y

  h_vals = approx(h_inv, xvals, xout=y)$y
  out[i-1] = -1 * log(mean(abs(h_vals - weights), na.rm=TRUE),2)
}
plot(2:N, out, pch=19, cex=2)
abline(0,1)

vals = -1 * log(abs(h_vals - weights),2)
vals = vals[!is.na(vals)]

################
library(baseN)
x = readRDS("output/bucket_allData.Rds")
z = readRDS("output/bucket_allData_bhash.Rds")
h = as.numeric(z)

weights = cumsum(x[,1]) / sum(x[,1])

index = sort(sample(1:length(z), 10000L))
plot(h[index], weights[index], type="l")



