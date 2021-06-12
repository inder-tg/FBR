
library(raster)
library(mapview)
library(RColorBrewer)

##----------- Preparación de datos ----------##

listFiles <- list.files(path = "C:/Users/inder/OneDrive/Desktop/veranoCiencia2020/data/",
                        pattern = ".tif", full.names = TRUE)

# --- shp file

rootDir <- getwd()

shpDir <- paste0("C:/Users/inder/OneDrive/Desktop/veranoCiencia2020/data/master_MN_LLC")

setwd(shpDir)

dsn <- system.file("vectors", package = "rgdal")

shpMarismas <- shapefile('master_MN_LCC')

shpCalidad <- paste0("C:/Users/inder/OneDrive/Desktop/veranoCiencia2020/data/puntos_MN")

setwd(shpCalidad)

shpCalidad <- shapefile('MN_points_2008082')

calidad_raster <- rasterize(shpCalidad, indice, field="Qual_num")

# plot(calidad_raster)

# ---

zona <- raster(listFiles[1])
fraction <- raster(listFiles[2])
indice <- raster(listFiles[3])

fraction_zona <- fraction
fraction_zona[zona != 1] <- NA
indice_zona <- indice
indice_zona[zona != 1] <- NA

bestQuality <- calidad_raster
bestQuality[calidad_raster!=0] <- NA
bestQuality[zona!=1] <- NA
bestQuality[fraction_zona==0] <- NA

goodQuality <- calidad_raster
goodQuality[calidad_raster!=1] <- NA
goodQuality[zona!=1] <- NA
goodQuality[fraction_zona==0] <- NA

fraction_bestQuality_nonzero <- fraction_zona
fraction_bestQuality_nonzero[calidad_raster!=0] <- NA
fraction_bestQuality_nonzero[fraction_zona==0] <- NA

fraction_goodQuality_nonzero <- fraction_zona
fraction_goodQuality_nonzero[calidad_raster!=1] <- NA
fraction_goodQuality_nonzero[fraction_zona==0] <- NA

indice_bestQuality_nonzero <- indice_zona
indice_bestQuality_nonzero[calidad_raster!=0] <- NA
indice_bestQuality_nonzero[fraction_zona==0] <- NA

indice_goodQuality_nonzero <- indice_zona
indice_goodQuality_nonzero[calidad_raster!=1] <- NA
indice_goodQuality_nonzero[fraction_zona==0] <- NA

# ---

# mpFraction <- mapview(fraction_zona_nozero, map.types="Esri.WorldImagery",
#                       legend=FALSE, layer.name="",homebutton=FALSE, na.color="transparent")

# ---

mpShp <- mapview(shpMarismas,legend=FALSE, layer.name="shp",
                 color="yellow",lwd=4,
                 # col.regions="transparent",
                 alpha.regions=0,homebutton=FALSE,
                 map.types="Esri.WorldImagery")

mpBest <- mapview(bestQuality,legend=FALSE,layer.name="mejorCalidad",
                  col="gray",map.types="Esri.WorldImagery",
                  na.color="transparent", homebutton=FALSE)

mpGood <- mapview(goodQuality,legend=FALSE,layer.name="buenaCalidad",
                  col="gray",map.types="Esri.WorldImagery",
                  na.color="transparent", homebutton=FALSE)

mpBest+mpShp
mpGood+mpShp
mpBest+mpGood+mpShp


mpFractionBest <- mapview(fraction_bestQuality_nonzero*1e-2,
                          col.regions=rev(brewer.pal(9,"Blues"))[1:5],
                          map.types="Esri.WorldImagery",
                          # layer.name="Fracción",
                          na.color="transparent")

mpFractionGood <- mapview(fraction_goodQuality_nonzero*1e-2,
                          col.regions=rev(brewer.pal(9,"Reds"))[1:5],
                          map.types="Esri.WorldImagery",
                          # layer.name="Fracción",
                          na.color="transparent")

mpFractionBest+mpShp
mpFractionGood+mpShp
mpFractionBest+mpFractionGood+mpShp

mpIndiceBest <- mapview(indice_bestQuality_nonzero,
                        col.regions=rev(brewer.pal(9,"Reds"))[1:5],
                        map.types="Esri.WorldImagery",
                        layer.name="Índice",
                        na.color="transparent")

mpIndiceGood <- mapview(indice_goodQuality_nonzero,
                        col.regions=rev(brewer.pal(9,"Reds"))[1:5],
                        map.types="Esri.WorldImagery",
                        layer.name="Índice",
                        na.color="transparent")

mpIndiceBest+mpShp
mpIndiceGood+mpShp
mpIndiceBest+mpIndiceGood+mpShp

# --- plot for paper

mpIndiceGood <- mapview(indice_goodQuality_nonzero,
                        col.regions=rev(brewer.pal(9,"Reds"))[1:5],
                        map.types="Esri.WorldImagery",
                        layer.name="MNDWI6",
                        na.color="transparent")

mpIndiceGood+mpShp

mpFractionGood <- mapview(fraction_goodQuality_nonzero*1e-2,
                          col.regions=rev(brewer.pal(9,"Blues"))[1:5],
                          map.types="Esri.WorldImagery",
                          layer.name="Water fraction",
                          na.color="transparent")

mpFractionGood+mpShp
