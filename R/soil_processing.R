library(soilDB)
library(sf)


ssurgo2sp <- function (id, mapunit = NULL, component = NULL, chorizon = NULL,
                       mapunit.shp = NULL, nmapunit = 1, nsoil = 1, xout = NULL,
                       soil.bottom = 200, method = c("constant", "linear"), nlayers = 10,
                       verbose = FALSE) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    warning("sf package is required for this function")
    return(NULL)
  }
  mapunit2 <- subset(mapunit, select = c("musym", "muname",
                                         "muacres", "farmlndcl", "lkey", "mukey"))
  mapunit2$muacres.percent <- mapunit2$muacres/sum(mapunit2$muacres) *
    100
  mapunit2 <- mapunit2[order(mapunit2$muacres, decreasing = TRUE),
  ]
  mapunit3 <- mapunit2[seq_len(nmapunit), ]
  component2 <- subset(component, select = c("compname", "comppct.r",
                                             "slope.r", "drainagecl", "elev.r", "taxsubgrp", "taxpartsize",
                                             "taxclname", "geomdesc", "mukey", "cokey"))
  if (any(grepl("state", names(component))))
    component2$state <- component$state
  component5 <- NULL
  for (i in mapunit3$mukey) {
    component3 <- component2[component2$mukey == i, ]
    component4 <- component3[order(component3$comppct.r,
                                   decreasing = TRUE), ]
    component4$acres.proportion <- mapunit3[mapunit3$mukey ==
                                              i, ]$muacres.percent/100 * component4$comppct.r/100
    component4$compname.mukey <- as.factor(component4$compname):as.factor(component4$mukey)
    if (inherits(mapunit.shp, "sf")) {
      lonlat <- colMeans(sf::st_coordinates(mapunit.shp["MUKEY" ==
                                                          i]))
      mapunit.shp.d <- as.data.frame(mapunit.shp)
      if (any(grepl("AREASYMBOL", names(mapunit.shp.d)))) {
        component4$state <- rep(unique(strtrim(as.character(mapunit.shp.d[mapunit.shp.d$MUKEY ==
                                                                            i, "AREASYMBOL"]), 2)), nrow(component4))
      }
    }
    else {
      stop("mapunit.shp should be of class sf")
    }
    component4$longitude <- lonlat[1]
    component4$latitude <- lonlat[2]
    component5 <- rbind(component5, component4[seq_len(nsoil),
    ])
  }
  chorizon2 <- chorizon[chorizon$cokey %in% component5$cokey,
  ]
  if (nrow(chorizon2) < 1) {
    stop("component horizon does not match component cokey")
  }
  chorizon3 <- subset(chorizon2, select = c("hzname", "hzdept.r",
                                            "hzdepb.r", "hzthk.r", "sandtotal.r", "silttotal.r",
                                            "claytotal.r", "om.r", "dbthirdbar.r", "partdensity", "ksat.r", "awc.r",
                                            "wthirdbar.r", "wfifteenbar.r", "wsatiated.r", "ph1to1h2o.r",
                                            "pbray1.r", "fragvol.r", "ph1to1h2o.r", "ph01mcacl2.r", "cec7.r",
                                            "cokey", "chkey"))
  chorizon3$hzdepa.r <- (chorizon3$hzdept.r + chorizon3$hzdepb.r)/2
  chorizon3$Depth <- chorizon3$hzdepa.r
  if (any(is.na(chorizon3$hzthk.r))) {
    chorizon3$Thickness <- chorizon3$hzdepb.r - chorizon3$hzdept.r
  }
  else {
    chorizon3$Thickness <- chorizon3$hzthk.r
  }
  if (any(is.na(chorizon3$wthirdbar.r))) {
    chorizon3$DUL <- apsimx:::sr_dul(chorizon3$claytotal.r, chorizon3$sandtotal.r,
                                     chorizon3$om.r)
  }
  else {
    chorizon3$DUL <- chorizon3$wthirdbar.r * 0.01
  }
  if (any(is.na(chorizon3$wfifteenbar.r))) {
    chorizon3$LL15 <- apsimx:::sr_ll(chorizon3$claytotal.r, chorizon3$sandtotal.r,
                                     chorizon3$om.r)
  }
  else {
    chorizon3$LL15 <- chorizon3$wfifteenbar.r * 0.01
  }
  if (any(is.na(chorizon3$wsatiated.r))) {
    DUL_S <- apsimx:::sr_dul_s(chorizon3$claytotal.r, chorizon3$sandtotal.r,
                               chorizon3$om.r)
    chorizon3$SAT <- apsimx:::sr_sat(chorizon3$sandtotal.r, chorizon3$DUL,
                                     DUL_S)
  }
  else {
    chorizon3$SAT <- chorizon3$wsatiated.r * 0.01
  }
  chorizon3$KS <- chorizon3$ksat.r * (60 * 60) * 1e-4
  chorizon3$PH <- chorizon3$ph1to1h2o.r
  # chorizon3$BD <- (1 - chorizon3$SAT) * ifelse(is.na(chorizon3$partdensity),
  #                                              2.65, chorizon3$partdensity)
  chorizon3$BD <- chorizon3$dbthirdbar.r
  chorizon3$AirDry <- chorizon3$LL15 * ifelse(chorizon3$hzdept.r ==
                                                0, 0.5, 1)
  chorizon3$Carbon <- chorizon3$om.r * (1/1.72)
  chorizon3$ParticleSizeClay <- chorizon3$claytotal.r
  chorizon3$ParticleSizeSilt <- chorizon3$silttotal.r
  chorizon3$ParticleSizeSand <- chorizon3$sandtotal.r
  soil.names <- component5$compname.mukey
  soil.list <- vector(mode = "list", length = length(soil.names))
  names(soil.list) <- soil.names
  vars <- c("Depth", "Thickness", "BD", "AirDry", "LL15", "DUL",
            "SAT", "KS", "Carbon", "PH", "ParticleSizeClay", "ParticleSizeSilt",
            "ParticleSizeSand")
  for (sz in 1:length(soil.names)) {
    one.profile <- chorizon3[chorizon3$cokey == component5$cokey[sz],]
    one.component <- component[component$cokey == component5$cokey[sz],]
    if (is.na(soil.bottom))
      soil.bottom <- max(one.profile$hzdepb.r)
    if (nrow(one.profile) < 1) {
      stop("There is no soil horizon for this component")
    }
    
    soil.d <- as.data.frame(one.profile)
    institute <- "SSURGO"
    site <- ""
    year <- 2024
    # id <- as.numeric(format(Sys.Date(), "%j"))
    attr(soil.d, "institute") <- institute
    attr(soil.d, "year") <- year
    attr(soil.d, "id") <- id
    attr(soil.d, "pedon") <- .generatePEDON(.institute = institute,
                                            .site = site,
                                            .year = year,
                                            .id = id)
    attr(soil.d, "source") <- "ssurgo_v1"
    attr(soil.d, "texture") <- .getSoilTexture(mean(soil.d$ParticleSizeSand, na.rm = T),
                                               mean(soil.d$ParticleSizeSilt, na.rm = T),
                                               mean(soil.d$ParticleSizeClay, na.rm = T))
    attr(soil.d, "depth") <- max(soil.d$hzdepb.r)
    attr(soil.d, "description") <- paste("cokey =",
                                         component5$cokey[sz], "- acres percent =", component5$acres.proportion[sz] *
                                           100, "- component percent =", component5$comppct.r[sz],
                                         "- taxonomic classification name =", as.character(component5$taxclname)[sz],
                                         "- drainage class =", as.character(component5$drainagecl)[sz],
                                         "- elevation =", component5$elev.r[sz], "- slope =",
                                         component5$slope.r[sz], "- geomdesc =", as.character(component5$geomdesc)[sz])
    attr(soil.d, "site") <- site
    attr(soil.d, "country") <- "USA"
    attr(soil.d, "latitude") <- component5$latitude[1]
    attr(soil.d, "longitude") <- component5$longitude[1]
    attr(soil.d, "scs_family") <- one.component$taxclname
    attr(soil.d, "scom") <- NA_character_
    attr(soil.d, "salb") <- one.component$albedodry.r
    attr(soil.d, "slu1") <- NA_real_
    attr(soil.d, "sldr") <- 0.47 #Mean value for DSSAT soil default dataset
    attr(soil.d, "slro") <- 72 #Mean value for DSSAT soil default dataset
    attr(soil.d, "slnf") <- 1 #Mean value for DSSAT soil default dataset
    attr(soil.d, "slpf") <- 1 #Mean value for DSSAT soil default dataset
    attr(soil.d, "smhb") <- "SA009"
    attr(soil.d, "smpx") <- "SA002"
    attr(soil.d, "smke") <- NA_real_
    
    soil.d_attributes <- attributes(soil.d)
    
    names(soil.d)[names(soil.d) == "hzdepb.r"] <- "slb"
    soil.d$slmh <- soil.d$slb
    names(soil.d)[names(soil.d) == "LL15"] <- "slll"
    names(soil.d)[names(soil.d) == "DUL"] <- "sdul"
    names(soil.d)[names(soil.d) == "SAT"] <- "ssat"
    soil.d$srgf <- 1
    names(soil.d)[names(soil.d) == "KS"] <- "ssks"
    names(soil.d)[names(soil.d) == "BD"] <- "sbdm"
    names(soil.d)[names(soil.d) == "om.r"] <- "sloc"
    names(soil.d)[names(soil.d) == "claytotal.r"] <- "slcl"
    names(soil.d)[names(soil.d) == "silttotal.r"] <- "slsi"
    names(soil.d)[names(soil.d) == "fragvol.r"] <- "slcf"
    soil.d$slni <- rep(NA_real_, nrow(soil.d))
    names(soil.d)[names(soil.d) == "ph1to1h2o.r"] <- "slhw"
    names(soil.d)[names(soil.d) == "ph01mcacl2.r"] <- "slhb"
    names(soil.d)[names(soil.d) == "cec7.r"] <- "scec"
    soil.d$sadc <- NA_real_
    
    vars <- c("slb", "slmh", "slll", "sdul", "ssat", "srgf", "ssks", "sbdm",
              "sloc", "slcl", "slsi", "slcf", "slni", "slhw", "slhb", "scec")
    soil.d <- soil.d[,vars]
    
    soil.d_attributes$names <- vars
    
    names(soil.d) <- vars
    attributes(soil.d) <- soil.d_attributes
    soil.list[[sz]] <- soil.d
  }
  return(soil.list)
}


get_ssurgo_soil_profile <- function (id, lonlat, shift = -1, nmapunit = 1, nsoil = 1, xout = NULL,
                                     soil.bottom = 200, method = c("constant", "linear"), nlayers = 10,
                                     check = TRUE, fix = FALSE, verbose = FALSE, xargs = NULL) {
  if (!requireNamespace("soilDB", quietly = TRUE)) {
    stop("The soilDB package is required for this function")
    return(NULL)
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("The sp package is required for this function")
    return(NULL)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("The sf package is required for this function")
    return(NULL)
  }
  if (!requireNamespace("spData", quietly = TRUE)) {
    stop("The spData package is required for this function")
    return(NULL)
  }
  if (length(lonlat) != 2 || !is.numeric(lonlat))
    stop("lonlat should be a vector with length equal to 2")
  lon <- lonlat[1]
  lat <- lonlat[2]
  if (requireNamespace("maps", quietly = TRUE)) {
    country <- maps::map.where(x = lon, y = lat)
    if (country != "USA" || is.na(country))
      stop("These coordinates do not correspond to a location in the USA. \n Did you specify the coordinates correctly?")
  }
  if (shift <= 0) {
    spg <- sp::SpatialPoints(cbind(x = lon, y = lat), proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  }
  else {
    shift <- (shift/111) * 0.001
    lonlat.mat <- rbind(lonlat, lonlat + c(shift * 0.75,
                                           0), lonlat + c(shift * 0.75, shift), lonlat + c(0,
                                                                                           shift), lonlat)
    rownames(lonlat.mat) <- NULL
    pg <- sp::Polygon(lonlat.mat)
    spg <- sp::SpatialPolygons(list(sp::Polygons(list(pg),
                                                 "s1")), proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  }
  if (verbose == FALSE) {
    res <- suppressWarnings(soilDB::SDA_spatialQuery(spg,
                                                     what = "mupolygon", geomIntersection = TRUE))
  }
  else {
    res <- soilDB::SDA_spatialQuery(spg, what = "mupolygon",
                                    geomIntersection = TRUE)
  }
  mu.is <- soilDB::format_SQL_in_statement(res$mukey)
  sql <- sprintf("mukey IN %s", mu.is)
  if (verbose == FALSE) {
    fSDA <- suppressWarnings(suppressMessages(fetchSDA(sql, duplicates = TRUE)))
  }
  else {
    fSDA <- fetchSDA(sql, duplicates = TRUE)
  }
  if (nsoil < 0 || is.na(nsoil)) {
    nsoil <- nrow(fSDA@site)
  }
  if (verbose == FALSE) {
    mapunit <- suppressWarnings(suppressMessages(soilDB::get_mapunit_from_SDA(sql)))
  }
  else {
    mapunit <- soilDB::get_mapunit_from_SDA(sql)
  }
  cmpnt <- fSDA@site
  names(cmpnt) <- gsub("_", ".", names(cmpnt), fixed = TRUE)
  cmpnt$geomdesc <- cmpnt$geompos
  if (shift <= 0 || length(unique(mapunit$areasymbol)) == 1) {
    cmpnt$state <- unique(strtrim(mapunit$areasymbol, 2))
  }
  else {
    cmpnt$state <- NA
    warning("This area includes more than one state. \n            I have not though about how to get the state in this case. Please submit an issue \n            with a reproducible example to https://github.com/femiguez/apsimx/issues")
  }
  chrzns <- fSDA@horizons
  names(chrzns) <- gsub("_", ".", names(chrzns), fixed = TRUE)
  if (sum(grepl("partdensity", names(chrzns))) == 0)
    chrzns$partdensity <- NA
  if (sum(grepl("hzthk", names(chrzns))) == 0)
    chrzns$hzthk.r <- NA
  if (sum(grepl("wsatiated", names(chrzns))) == 0)
    chrzns$wsatiated.r <- NA
  if (sum(grepl("wfifteenbar", names(chrzns))) == 0)
    chrzns$wfifteenbar.r <- NA
  if (sum(grepl("wthirdbar", names(chrzns))) == 0)
    chrzns$wthirdbar.r <- NA
  if (shift <= 0) {
    spg.sf <- sf::st_as_sf(spg)
    spg.sf[["MUKEY"]] <- res$mukey
    spg.sf[["AREASYMBOL"]] <- mapunit$areasymbol
    mapunit.shp <- spg.sf
  }
  else {
    mapunit.shp <- sf::st_as_sf(res)
  }
  sp0 <- ssurgo2sp(id, mapunit = mapunit, component = cmpnt, chorizon = chrzns,
                   mapunit.shp = mapunit.shp, nmapunit = nmapunit, nsoil = nsoil,
                   xout = xout, soil.bottom = soil.bottom, method = method,
                   nlayers = nlayers, verbose = verbose)
  ans <- sp0
  
  return(ans)
}