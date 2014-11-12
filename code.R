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

library(sp)
library(rGeoHash)
library(sp)
fnames = dir("UScensus2010blk/data/", pattern="*.rda")

for(j in 1:length(fnames)) {
  z = load_obj(paste0("UScensus2010blk/data/", fnames[j]))

  mat = matrix(unlist(lapply(z@polygons, function(v) v@labpt)), ncol=2, byrow=TRUE,
                dimnames=list(1:nrow(z)))
  hash = toBitGeoHash(mat[,2], mat[,1], 50L)
  hash = cbind(as.integer(z@data[,c("P0010001")]), apply(hash, 2, as.integer))

  saveRDS(hash[hash[,1] > 0,], paste0("~/Desktop/output/bucket_", j, ".Rds"))
}
