library(soilDB)
library(aqp)
library(spData)
library(sf)


get_component_from_SDA <- function (WHERE = NULL, duplicates = FALSE, childs = TRUE, droplevels = TRUE,
                                    nullFragsAreZero = TRUE, stringsAsFactors = NULL)
{
  if (!missing(stringsAsFactors) && is.logical(stringsAsFactors)) {
    .Deprecated(msg = sprintf("stringsAsFactors argument is deprecated.\nSetting package option with `NASISDomainsAsFactor(%s)`",
                              stringsAsFactors))
    NASISDomainsAsFactor(stringsAsFactors)
  }
  if (!duplicates & grepl(WHERE, pattern = "mukey")[1])
    warning("duplicates is set to FALSE and 'mukey' is in WHERE clause. Note: 'mukey' omitted from result.",
            call. = FALSE)
  es.vars <- "ecoclasstypename, ecoclassref, ecoclassid, ecoclassname"
  co.vars <- "co.cokey, compname, comppct_r, compkind, majcompflag, localphase, runoff, drainagecl, hydricrating, erocl, earthcovkind1, earthcovkind2, elev_r, slope_r, aspectrep, albedodry_r, map_r, airtempa_r, reannualprecip_r, ffd_r, hydgrp,  nirrcapcl, nirrcapscl, irrcapcl, irrcapscl, tfact, wei, weg, corcon, corsteel, frostact, taxclname, taxorder, taxsuborder, taxgrtgroup, taxsubgrp, taxpartsize, taxpartsizemod, taxceactcl, taxreaction, taxtempcl, taxmoistscl, taxtempregime, soiltaxedition"
  vars <- paste0(unlist(strsplit(co.vars, "earthcovkind2,")),
                 collapse = paste0("earthcovkind2, ", es.vars, ","))
  q.component <- paste("SELECT", if (duplicates == FALSE)
    "DISTINCT"
    else "mu.mukey AS mukey,", "mu.nationalmusym,", vars, "FROM\n     legend  l                      INNER JOIN\n     mapunit mu ON mu.lkey = l.lkey INNER JOIN",
    if (duplicates == FALSE) {
      "(SELECT nationalmusym AS nationalmusym2, MIN(mukey) AS mukey2\n      FROM mapunit\n      GROUP BY nationalmusym\n     ) AS mu2 ON mu2.nationalmusym2 = mu.nationalmusym INNER JOIN"
    }
    else "", "(SELECT", vars, ", mukey AS mukey2\n      FROM\n      component  co                        LEFT OUTER JOIN\n      coecoclass ce ON ce.cokey = co.cokey AND\n                       ecoclasstypename IN ('NRCS Rangeland Site', 'NRCS Forestland Site')\n      ) AS co ON co.mukey2 =",
    if (duplicates == FALSE)
      "mu2.mukey2"
    else "mu.mukey", "WHERE", WHERE, "ORDER BY nationalmusym, comppct_r DESC, compname;")
  d.component <- soilDB::SDA_query(q.component)
  if (length(d.component) == 0)
    return(d.component)
  d.component <- soilDB::uncode(d.component, droplevels = droplevels)
  if (duplicates == FALSE) {
    d.component <- cbind(mukey = NA, d.component)
  }
  q.pm <- paste0("SELECT co.cokey, pmg.copmgrpkey, pmgroupname, pmorder, pmkind, pmorigin\n\n    FROM\n    component co                                        LEFT OUTER JOIN\n    copmgrp   pmg ON pmg.cokey         = co.cokey AND\n                       pmg.rvindicator = 'Yes'          LEFT OUTER JOIN\n    copm      pm  ON pm.copmgrpkey     = pmg.copmgrpkey\n\n    WHERE co.cokey IN (",
                 paste0(d.component$cokey, collapse = ", "), ")\n\n    ORDER BY co.cokey, pmg.copmgrpkey, pmorder")
  q.lf <- paste0("SELECT co.cokey, ls.geomfname landscape, lf.geomfname landform, lf.geomfeatid, lf.existsonfeat,\n            geomposmntn mntn, geomposhill hill, geompostrce trce, geomposflats flats,\n            shapeacross, shapedown,\n            hillslopeprof\n\n    FROM\n    component co\n\n    LEFT OUTER JOIN\n        cogeomordesc ls ON ls.cokey       = co.cokey   AND\n                           ls.rvindicator = 'Yes'      AND\n                           ls.geomftname  = 'Landscape'\n    LEFT OUTER JOIN\n        cogeomordesc lf ON lf.cokey       = co.cokey   AND\n                           lf.rvindicator = 'Yes'      AND\n                           lf.geomftname  = 'Landform'\n    LEFT OUTER JOIN\n        cosurfmorphgc  lf_3d ON lf_3d.cogeomdkey = lf.cogeomdkey\n    LEFT OUTER JOIN\n        cosurfmorphss  lf_ss ON lf_ss.cogeomdkey = lf.cogeomdkey\n    LEFT OUTER JOIN\n        cosurfmorphhpp lf_2d ON lf_2d.cogeomdkey = lf.cogeomdkey\n\n    WHERE co.cokey IN (",
                 paste0(d.component$cokey, collapse = ", "), ")\n\n    ORDER BY cokey, ls.geomftname, ls.geomfeatid, ls.existsonfeat, lf.geomftname, lf.geomfeatid, lf.existsonfeat")
  if (childs == TRUE) {
    d.pm <- soilDB::SDA_query(q.pm)
    d.cogmd <- soilDB::SDA_query(q.lf)
    d.cosrf <- soilDB:::.get_cosurffrags_from_SDA(unique(d.component$cokey),
                                                  nullFragsAreZero = nullFragsAreZero)
    d.cm <- soilDB::uncode(soilDB:::.get_comonth_from_SDA(d.component$cokey))
    d.pm <- soilDB:::.copm_prep(d.pm, db = "SDA")
    d.cogmd <- soilDB:::.cogmd_prep(d.cogmd, db = "SDA")
    d.component <- merge(d.component, d.pm, by = "cokey",
                         all.x = TRUE)
    d.component <- merge(d.component, d.cogmd, by = "cokey",
                         all.x = TRUE)
    d.component <- merge(d.component, d.cosrf, by = "cokey",
                         all.x = TRUE)
    d.component <- merge(d.component, d.cm, by = "cokey",
                         all.x = TRUE)
    idx <- grepl("surface_", names(d.component))
    if (nullFragsAreZero == nullFragsAreZero) {
      d.component[idx] <- lapply(d.component[idx], function(x) ifelse(is.na(x),
                                                                      0, x))
    }
  }
  idx <- table(d.component$cokey) > 1
  if (any(idx)) {
    cokeys <- as.integer(names(idx[idx == TRUE]))
    idx <- d.component$cokey %in% cokeys
    assign("component.ecosite.problems", value = cokeys,
           envir = soilDB::get_soilDB_env())
    message("-> QC: multiple ecosites linked to 1 component\n\tUse `get('component.ecosite.problems', envir = get_soilDB_env())` for component keys (cokey)")
    nodups <- {
      . <- d.component[idx, ]
      . <- split(., .$cokey)
      . <- lapply(., function(x) {
        temp = x[1, ]
        temp$ecoclassname = paste0(x$ecoclassname, collapse = ", ")
        temp$ecoclasstypename = paste0(x$ecoclasstypename,
                                       collapse = ", ")
        temp$ecoclassref = paste0(x$ecoclassref, collapse = ", ")
        temp$ecoclassid = paste0(x$ecoclassid, collapse = ", ")
        return(temp)
      })
      . <- do.call("rbind", .)
    }
    d.component <- rbind(d.component[!idx, ], nodups)
    d.component <- with(d.component, d.component[order(nationalmusym,
                                                       -comppct_r, compname), ])
  }
  return(d.component)
}


get_chorizon_from_SDA <- function (WHERE = NULL, duplicates = FALSE, childs = TRUE, nullFragsAreZero = TRUE,
                                   droplevels = TRUE, stringsAsFactors = NULL)
{
  if (!missing(stringsAsFactors) && is.logical(stringsAsFactors)) {
    .Deprecated(msg = sprintf("stringsAsFactors argument is deprecated.\nSetting package option with `NASISDomainsAsFactor(%s)`",
                              stringsAsFactors))
    NASISDomainsAsFactor(stringsAsFactors)
  }
  q.chorizon <- paste("\n  SELECT", if (duplicates == FALSE) {
    "DISTINCT"
  }, "hzname, hzdept_r, hzdepb_r, texture, texcl, lieutex,\n     fragvol_l, fragvol_r, fragvol_h, \n     sandtotal_l, sandtotal_r, sandtotal_h, \n     silttotal_l, silttotal_r, silttotal_h, \n     claytotal_l, claytotal_r, claytotal_h,\n     om_l, om_r, om_h, \n     dbthirdbar_l, dbthirdbar_r, dbthirdbar_h,\n     ksat_l, ksat_r, ksat_h,  \n     awc_l, awc_r, awc_h, \n     lep_r, sar_r, ec_r, cec7_r, sumbases_r, \n     ph1to1h2o_l, ph1to1h2o_r, ph1to1h2o_h,\n     ph01mcacl2_l, ph01mcacl2_r, ph01mcacl2_h,\n     pbray1_l, pbray1_r, pbray1_h,\n     caco3_l, caco3_r, caco3_h, \n     kwfact, kffact, c.cokey, ch.chkey\n  FROM legend l INNER JOIN\n       mapunit mu ON mu.lkey = l.lkey",
  if (duplicates == FALSE) {
    paste(" INNER JOIN\n  (SELECT MIN(nationalmusym) nationalmusym2, MIN(mukey) AS mukey2\n   FROM mapunit\n   GROUP BY nationalmusym) AS mu2 ON mu2.mukey2 = mu.mukey\n  ")
  }
  else {
    paste(" INNER JOIN\n   (SELECT nationalmusym, mukey\n    FROM mapunit) AS mu2 ON mu2.mukey = mu.mukey\n   ")
  }, "INNER JOIN\n   component    c    ON c.mukey      = mu.mukey   LEFT JOIN\n   chorizon     ch   ON ch.cokey     = c.cokey    LEFT OUTER JOIN\n   chtexturegrp chtg ON chtg.chkey   = ch.chkey AND rvindicator = 'Yes' RIGHT JOIN\n   chtexture    cht  ON cht.chtgkey  = chtg.chtgkey\n\n   LEFT OUTER JOIN\n       (SELECT SUM(fragvol_l) fragvol_l, SUM(fragvol_r) fragvol_r, SUM(fragvol_h) fragvol_h, ch2.chkey\n        FROM chorizon ch2\n        INNER JOIN chfrags chf ON chf.chkey = ch2.chkey\n\n        GROUP BY ch2.chkey) chfrags2  ON chfrags2.chkey = ch.chkey\n\n  WHERE",
  WHERE, "ORDER BY c.cokey, hzdept_r ASC;")
  d.chorizon <- soilDB::SDA_query(q.chorizon)
  metadata <- NULL
  load(system.file("data/metadata.rda", package = "soilDB")[1])
  if (!is.null(d.chorizon) && nrow(d.chorizon) > 0) {
    d.chorizon <- within(d.chorizon, {
      nationalmusym = NULL
      texture = tolower(texture)
      if (getOption("stringsAsFactors", default = FALSE) ||
          getOption("soilDB.NASIS.DomainsAsFactor", default = FALSE)) {
        texcl = factor(texcl, levels = metadata[metadata$ColumnPhysicalName ==
                                                  "texcl", "ChoiceLabel"])
      }
      if (droplevels == droplevels && is.factor(texcl)) {
        texcl = droplevels(texcl)
      }
    })
    if (childs == TRUE) {
      WHERE = paste0("WHERE co.cokey IN (", paste0(unique(d.chorizon$cokey),
                                                   collapse = ","), ")")
      q.chfrags <- paste("\n  \n                            -- find fragsize_r\n                            CREATE TABLE #RF1 (cokey INT, chkey INT, chfragskey INT, fragvol_r REAL,\n                            shape CHAR(7), para INT, nonpara INT, fragsize_r2 INT);\n  \n                            INSERT INTO  #RF1 (cokey, chkey, chfragskey, fragvol_r, shape, para, nonpara, fragsize_r2)\n                            SELECT             co.cokey, ch.chkey, chfragskey, fragvol_r,\n                            -- shape\n                            CASE WHEN fragshp = 'Nonflat' OR fragshp IS NULL THEN 'nonflat' ELSE 'flat' END shape,\n                            -- hardness\n                            CASE WHEN fraghard IN ('Extremely weakly cemented', 'Very weakly cemented', 'Weakly cemented', 'Weakly cemented*', 'Moderately cemented', 'Moderately cemented*', 'soft')                     THEN 1 ELSE NULL END para,\n                            CASE WHEN fraghard IN ('Strongly cemented', 'Strongly cemented*', 'Very strongly cemented', 'Extremely strong', 'Indurated', 'hard')   OR fraghard IS NULL THEN 1 ELSE NULL END nonpara,\n                            -- fragsize_r\n                            CASE WHEN fragsize_r IS NOT NULL THEN fragsize_r\n                            WHEN fragsize_r IS NULL     AND fragsize_h IS NOT NULL AND fragsize_l IS NOT NULL\n                            THEN (fragsize_h + fragsize_l) / 2\n                            WHEN fragsize_h IS NOT NULL THEN fragsize_h\n                            WHEN fragsize_l IS NOT NULL THEN fragsize_l\n                            ELSE NULL END\n                            fragsize_r2\n  \n                            FROM\n                            component co                        LEFT OUTER JOIN\n                            chorizon  ch ON ch.cokey = co.cokey LEFT OUTER JOIN\n                            chfrags   cf ON cf.chkey = ch.chkey",
                         WHERE = WHERE, "\n                            ORDER BY co.cokey, ch.chkey, cf.chfragskey;\n  \n  \n                            -- compute logicals\n                            CREATE TABLE #RF2 (\n                            cokey INT, chkey INT, chfragskey INT, fragvol_r REAL, para INT, nonpara INT,\n                            fine_gravel INT, gravel INT, cobbles INT, stones INT, boulders INT, channers INT, flagstones INT,\n                            unspecified INT\n                            );\n                            INSERT INTO  #RF2 (\n                            cokey, chkey, chfragskey, fragvol_r, para, nonpara,\n                            fine_gravel, gravel, cobbles, stones, boulders, channers, flagstones,\n                            unspecified\n                            )\n                            SELECT\n                            cokey, chkey, chfragskey, fragvol_r, para, nonpara,\n                            -- fragments\n                            CASE WHEN   fragsize_r2 >= 2   AND fragsize_r2 < 5   AND shape = 'nonflat' THEN 1 ELSE NULL END fine_gravel,\n                            CASE WHEN   fragsize_r2 >= 2   AND fragsize_r2 < 75  AND shape = 'nonflat' THEN 1 ELSE NULL END gravel,\n                            CASE WHEN   fragsize_r2 >= 75  AND fragsize_r2 < 250 AND shape = 'nonflat' THEN 1 ELSE NULL END cobbles,\n                            CASE WHEN ((fragsize_r2 >= 250 AND fragsize_r2 < 600 AND shape = 'nonflat') OR\n                            (fragsize_r2 >= 380 AND fragsize_r2 <= 600 AND shape = 'flat'))\n                            THEN 1 ELSE NULL END stones,\n                            CASE WHEN   fragsize_r2 >= 600 THEN 1 ELSE NULL END boulders,\n                            CASE WHEN   fragsize_r2 >= 2   AND fragsize_r2 < 150 AND shape = 'flat' THEN 1 ELSE NULL END channers,\n                            CASE WHEN   fragsize_r2 >= 150 AND fragsize_r2 < 380 AND shape = 'flat' THEN 1 ELSE NULL END flagstones,\n                            CASE WHEN   fragsize_r2 IS NULL                                         THEN 1 ELSE NULL END unspecified\n  \n                            FROM\n                            #RF1\n  \n                            ORDER BY cokey, chkey, chfragskey;\n  \n  \n                            -- summarize rock fragments\n                            SELECT\n                            chkey,\n                            -- nonpara rock fragments\n                            SUM(fragvol_r * fine_gravel * nonpara)  fine_gravel,\n                            SUM(fragvol_r * gravel      * nonpara)  gravel,\n                            SUM(fragvol_r * cobbles     * nonpara)  cobbles,\n                            SUM(fragvol_r * stones      * nonpara)  stones,\n                            SUM(fragvol_r * boulders    * nonpara)  boulders,\n                            SUM(fragvol_r * channers    * nonpara)  channers,\n                            SUM(fragvol_r * flagstones  * nonpara)  flagstones,\n                            -- para rock fragments\n                            SUM(fragvol_r * fine_gravel * para)     parafine_gravel,\n                            SUM(fragvol_r * gravel      * para)     paragravel,\n                            SUM(fragvol_r * cobbles     * para)     paracobbles,\n                            SUM(fragvol_r * stones      * para)     parastones,\n                            SUM(fragvol_r * boulders    * para)     paraboulders,\n                            SUM(fragvol_r * channers    * para)     parachanners,\n                            SUM(fragvol_r * flagstones  * para)     paraflagstones,\n                            -- unspecified\n                            SUM(fragvol_r * unspecified)            unspecified,\n                            -- total_frags_pct_para\n                            SUM(fragvol_r               * nonpara)  total_frags_pct_nopf,\n                            -- total_frags_pct\n                            SUM(fragvol_r)                          total_frags_pct\n  \n                            FROM #RF2\n  \n                            GROUP BY cokey, chkey\n  \n                            ORDER BY cokey, chkey;\n  \n  \n                            -- cleanup\n                            DROP TABLE #RF1;\n                            DROP TABLE #RF2;\n                            ")
      d.chfrags <- soilDB::SDA_query(q.chfrags)
      idx <- !names(d.chfrags) %in% "chkey"
      if (nullFragsAreZero == TRUE) {
        d.chfrags[idx] <- lapply(d.chfrags[idx], function(x) ifelse(is.na(x),
                                                                    0, x))
      }
      d.chorizon <- merge(d.chorizon, d.chfrags, all.x = TRUE,
                          by = "chkey")
    }
  }
  return(d.chorizon)
}



fetchSDA <- function (WHERE = NULL, duplicates = FALSE, childs = TRUE, nullFragsAreZero = TRUE,
                      rmHzErrors = FALSE, droplevels = TRUE, stringsAsFactors = NULL)
{
  if (!missing(stringsAsFactors) && is.logical(stringsAsFactors)) {
    .Deprecated(msg = sprintf("stringsAsFactors argument is deprecated.\nSetting package option with `NASISDomainsAsFactor(%s)`",
                              stringsAsFactors))
    NASISDomainsAsFactor(stringsAsFactors)
  }
  f.component <- get_component_from_SDA(WHERE, duplicates = duplicates,
                                        childs = childs, droplevels = droplevels, nullFragsAreZero = TRUE)
  if (is.null(f.component)) {
    stop("WHERE clause returned no components.", call. = FALSE)
  }
  f.chorizon <- get_chorizon_from_SDA(paste0("c.cokey IN",
                                             soilDB::format_SQL_in_statement(unique(f.component$cokey))),
                                      duplicates = duplicates, droplevels = droplevels)
  f.mapunit <- soilDB::get_mapunit_from_SDA(paste("mu.nationalmusym IN",
                                                  soilDB::format_SQL_in_statement(unique(f.component$nationalmusym))))
  f.diag <- soilDB:::.get_diagnostics_from_SDA(f.component$cokey)
  f.restr <- soilDB:::.get_restrictions_from_SDA(f.component$cokey)
  if (rmHzErrors) {
    f.chorizon.test <- aqp::checkHzDepthLogic(f.chorizon,
                                              hzdepths = c("hzdept_r", "hzdepb_r"), idname = "cokey",
                                              fast = TRUE)
    good.ids <- as.character(f.chorizon.test$cokey[which(f.chorizon.test$valid)])
    bad.ids <- as.character(f.chorizon.test$cokey[which(!f.chorizon.test$valid)])
    f.chorizon <- f.chorizon[which(f.chorizon$cokey %in%
                                     good.ids), ]
    if (length(bad.ids) > 0)
      assign("component.hz.problems", value = bad.ids,
             envir = soilDB::get_soilDB_env())
  }
  aqp::depths(f.chorizon) <- cokey ~ hzdept_r + hzdepb_r
  aqp::site(f.chorizon) <- f.component
  aqp::site(f.chorizon) <- f.mapunit
  if ("chkey" %in% aqp::horizonNames(f.chorizon) && all(!is.na("chkey"))) {
    aqp::hzidname(f.chorizon) <- "chkey"
  }
  aqp::hzdesgnname(f.chorizon) <- "hzname"
  aqp::hztexclname(f.chorizon) <- "texture"
  if (is.data.frame(f.diag)) {
    aqp::diagnostic_hz(f.chorizon) <- f.diag
  }
  if (is.data.frame(f.restr)) {
    aqp::restrictions(f.chorizon) <- f.restr
  }
  if (exists("component.hz.problems", envir = soilDB::get_soilDB_env()))
    message("-> QC: horizon errors detected, use `get('component.hz.problems', envir=get_soilDB_env())` for component keys (cokey)")
  return(f.chorizon)
}

get_ssurgo_tables <- function (lonlat, shift = -1, aoi, verbose = FALSE)
{
  if (!requireNamespace("soilDB", quietly = TRUE)) {
    warning("The soilDB package is required for this function")
    return(NULL)
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    warning("The sp package is required for this function")
    return(NULL)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    warning("The sf package is required for this function")
    return(NULL)
  }
  if (!requireNamespace("spData", quietly = TRUE)) {
    warning("The spData package is required for this function")
    return(NULL)
  }
  if (!missing(lonlat) && !missing(aoi))
    stop("Either use 'lonlat' or 'aoi', but not both", call. = FALSE)
  if (missing(aoi)) {
    lon <- lonlat[1]
    lat <- lonlat[2]
    if (lon < -180 || lon > 180)
      stop("longitude should be between -180 and 180")
    if (lat < -90 || lat > 90)
      stop("latitude should be between -90 and 90")
    if (shift <= 0) {
      spg <- sp::SpatialPoints(cbind(x = lon, y = lat),
                               proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
    }
    else {
      shift <- (shift/111) * 0.001
      lonlat.mat <- rbind(lonlat, lonlat + c(shift * 0.75,
                                             0), lonlat + c(shift * 0.75, shift), lonlat +
                            c(0, shift), lonlat)
      rownames(lonlat.mat) <- NULL
      pg <- sp::Polygon(lonlat.mat)
      spg <- sp::SpatialPolygons(list(sp::Polygons(list(pg),
                                                   "s1")), proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
    }
  }
  else {
    if (inherits(aoi, "sf"))
      aoi <- sf::as_Spatial(aoi, "Spatial")
    if (!inherits(aoi, "SpatialPolygons"))
      stop("'aoi' should be of class 'SpatialPolygons'.",
           call. = TRUE)
    spg <- aoi
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
    fSDA <- suppressWarnings(suppressMessages(fetchSDA(sql,
                                                       duplicates = TRUE)))
  }
  else {
    fSDA <- fetchSDA(sql, duplicates = TRUE)
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
  states <- spData::us_states
  states <- sf::st_transform(states, crs = 3857)
  pspg <- sf::st_transform(sf::st_as_sf(spg), crs = 3857)
  ii <- as.integer(sf::st_intersects(pspg, states))
  cmpnt$state <- states[["NAME"]][[ii]]
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
  if (shift <= 0 && missing(aoi)) {
    spg.sf <- sf::st_as_sf(spg)
    spg.sf[["MUKEY"]] <- res$mukey
    spg.sf[["AREASYMBOL"]] <- mapunit$areasymbol
    mapunit.shp <- spg.sf
  }
  else {
    mapunit.shp <- sf::st_as_sf(res)
  }
  ans <- list(mapunit = mapunit, component = cmpnt, chorizon = chrzns,
              mapunit.shp = mapunit.shp)
  ans
}
