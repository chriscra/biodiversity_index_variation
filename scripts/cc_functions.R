### Functions used for Zambia biodiversity project

# ------------------------- #
# Install only missing packages
# ------------------------- #
install_missing_packages <- function(list_of_packages) {
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[ , "Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, repo = 'https://cloud.r-project.org/')
  }
  sapply(list_of_packages, library, character.only = TRUE)
}

# sizes of things in the environment
# run this as env_size(ls())
env_size <- function(workspace) {
  size = 0
  for (x in workspace){
    thisSize = object_size(get(x))
    size = size + thisSize
    message(x, " = ", appendLF = F); print(thisSize, units='auto')
  }
  message("total workspace is ", appendLF = F); print(size, units='auto')
}



invert = function(x){
  return(1-x)}

normalize <- function(x) {
  (x - cellStats(x,"min")) / (cellStats(x,"max") - cellStats(x,"min"))
}


gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                             pypath=NULL, readpoly=TRUE, quiet=TRUE) {
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- "/Users/christophercrawford/Google Drive/_Projects/Zambia/agroEcoTradeoff/python/gdal_polygonize.py" #Sys.which('gdal_polygonize.py')
  }
  if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists))
      stop(sprintf('File already exists: %s',
                   toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                  pypath, rastpath, gdalformat, outshape)))
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
    return(shp)
  }
  return(NULL)
}

check_crs <- function (x, y) {
  x_crs <- if (class(x)[1] == "sf") {st_crs(x)$proj4string} else {proj4string(x)}
  y_crs <- if (class(y)[1] == "sf") {st_crs(y)$proj4string} else {proj4string(y)}
  if (x_crs == y_crs) {
    "Match: the two projections are the same. Woot!"
  } else {
    "Don't Match: the two projections are not the same. Try reprojecting."
    }
}

cc_write_reload_raster <- function(raster, name, directory) {
  # This function saves and reloads a raster,
  # 1. sets the layer name of that raster to the input character vector "name"
  # 2. creates the file path from that name
  # 3. writes the raster to that filepath
  # 4. Finally, reloads the raster

  names(raster) <- name
  #names(named_raster) <- deparse(quote(named_raster)) # failed usage
  path <- paste0(directory, .Platform$file.sep, names(raster),".tif")
  writeRaster(raster, path, overwrite=TRUE)
  raster(path)
}

# glwd_africa_r <- cc_write_reload_raster(glwd_africa_r, "glwd_africa_r", p_basemaps)



cc_make_valid <- function (sf, add_reasons = TRUE) {
  if (class(sf)[1] != "sf") {stop("The function requires the input to be an sf object.")}# 0. Check to make sure input is an sf object.
  if (add_reasons == FALSE) {
    sf <- st_make_valid(sf) # 2a. run st_make_valid to make the geometries valid. This usually takes a long time.
    return(sf)
  } else {
    sf <- mutate(sf, pre_fix_reasons = st_is_valid(sf, reason = TRUE))
    # 1. Add new column listing whether each geometry is valid or not (along with reasons why not; NA means corrupt geometry)

    sf <- st_make_valid(sf)
    # 2a. run st_make_valid to make the geometries valid. This usually takes a long time.

    sf <- mutate(sf, post_fix_reasons = st_is_valid(sf, reason = TRUE))
    # 3. Check geometries again, and add another column to reasons_df

    return(sf)
  }
}



cc_rescale <- function(raster, factor, run_mask = TRUE, mask, fun = mean) {
  raster_rescaled <- raster::aggregate(raster, fact = factor, fun = fun)
  raster_output <- disaggregate(raster_rescaled, fact = factor)

  if (run_mask) {
    raster_output_masked <- raster_output %>%
      crop(extent(mask)) %>%
      raster::mask(mask)
    raster_output_masked
  } else {
    raster_output
  }
}

cc_prep_zambia <- function(raster) {
  raster <- crop(raster, extent(msk))
  raster[is.na(raster)] <- 0 # replace NAs with 0s
  raster <- raster::mask(raster, msk)
}

# aggregate, then disaggregate:
# cc_rescale <- function(raster, factor, run_mask = TRUE, mask) {
#   output <- list(length(factor))
#   for (i in factor) {
#     raster_rescaled <- aggregate(raster, fact = factor[i])
#     raster_output <- disaggregate(raster_rescaled, fact = factor[i])
#     if (run_mask) {
#       raster_output_masked <- raster_output %>%
#         crop(extent(mask)) %>%
#         raster::mask(mask)
#       output[[i]] <- raster_output_masked
#     } else {
#       output[[i]] <- raster_output
#     }
#   }
# }

cc_percentilize <- function(input_raster, projection = "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") {
  dt <- as.data.table(input_raster, na.rm = TRUE)
  names(dt)[3] <- "layer"
  dt[, key := 1:.N] # set a key with the original row order (.N is the number of rows in the data.table)
  setorder(dt, layer) # change the order of the data.table, permanently, ordering from smallest to greatest.
  dt[, percentiles := 1:.N]  # create a new column listing the order from smallest to greatest.
  setorder(dt, key) # put back in original order
  dt[, c("layer", "key") := NULL]  # drop the key column
  dt_to_raster(dt, CRSobj = sp::CRS(projection))
}




# function to rasterize, etc.
cc_make_raster8 <- function (input_sf, run_filter = TRUE,
                             filter_presence = TRUE, presence_code = 1,
                             filter_origin = TRUE, origin_code = c(1, 2),
                             filter_marine = TRUE, marine_code = "False",
                             filter_seasonal = TRUE, seasonal_code = c(1, 2, 3),
                             filter_EX_EW = TRUE,
                             filter_category = TRUE, category_code = c("CR", "EN", "VU"),
                             filter_range_size = TRUE, range_threshold = 0.5,
                             filter_range_size_zam = TRUE,
                             odd_n_global, odd_n_zam,
                             filter_both = TRUE,
                             run_reproject = TRUE,
                             run_clip = FALSE, clip_area,
                             run_extract = FALSE, run_cast = FALSE,
                             projection = st_crs("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0
                                                 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                             fasterize_field = NULL, fasterize_fun = 'sum', #template_raster = msk,
                             run_mask = TRUE, mask = msk, prep_zambia = TRUE) {
  # if run_filter = FALSE, set all the filter switches to false. This is used to turn off this entire filter step, so that make_raster can be used for just regular polygon-to-raster conversion.
  if (run_filter == FALSE) {
    filter_presence = FALSE
    filter_origin = FALSE
    filter_marine = FALSE
    filter_seasonal = FALSE
    filter_EX_EW = FALSE
    filter_category = FALSE
    filter_range_size = FALSE
    filter_range_size_zam = FALSE
    filter_both = FALSE
  }

  # if run_clip is TRUE, then clip the polygons to the specified area, in clip_area, which are unioned using st_union. This should be the same projection as the input sf object.
  clipped_sf <- input_sf %>%
    {if (run_clip)        st_intersection(., st_union(clip_area)) else .} %>%
    {if (run_extract)     st_collection_extract(., "POLYGON") else .} %>%
    {if (run_cast)        st_cast(., "MULTIPOLYGON",) else .}

  # Species range maps need to be filtered before richness maps are calculated. In particular, we need to
  # ------ filter to only extant species, native or reintroduced species, filter out passage and
  # ------ seasonal occurence uncertain ranges, AND EW or EX species. Then, we filter based on threatened
  # ------ status, small-ranged, etc.

  # If filter switches are set to TRUE, the line will be run.
  # 1a. filter to presence, defaults to code 1 (extant species only). Other unused codes are:
  # ------ 2. Probably Extant (discontinued, ambiguous), 3. Possibly Extant,
  # ------ 4. Possibly Extinct, 5. Extinct, 6. Presence Uncertain.
  # 1b. filter to origin, defaults to codes 1 & 2 (native and reintroduced species respectively).
  # ------ Other unused codes are: 3. Introduced, 4. Vagrant, 5. Origin Uncertain, 6. Assisted Colonisation.
  # 1c. filter by marine, defaults to "FALSE" (all species but marine ones)
  # 1d. filter species with range seasonality. Defaults to c(1,2,3).
  # ------ Codes are: "Resident" (1), "Breeding" (2), "Non-breeding Season" (3),
  # ------ Passage (4), and Seasonal Occurrence Uncertain (5).
  # 1e. filter out Extinct or Extinct in the Wild species ("EX", "EW")

  # 2. filter by redlist category. Defaults to c("CR", "EN", "VU")),
  # ------ i.e. "threatened" species. Other categories include LC [Least Concern],
  # ------ EX [Extinct], EW [Extinct in the Wild], NT [Near Threatened], & DD [Data Deficient].
  # 3. filter by total global range size. Defaults to selecting small ranged species,
  # ------ i.e. those with ranges smaller than the global median. Note that the range_size_quantile
  # ------ was calculated by ordering the total range sizes from smallest to largest and
  # ------ normalizing to 1. Total range size is actually just the sum of the parts
  # ------ of the range we're using (i.e. extant, native, reintroduced, breeding,
  # ------ non-breeding, resident, non-marine, etc.)
  # 4. filter to only small-ranged species, based on total range *within* Zambia. So, this prioritizes
  # ------ species that are below the median in their Zambia clipped range.
  # 5. filter to any species that has a threatened redlist category or
  # ------ is small-ranged (total range size below median).

  filtered_all <- clipped_sf %>%
    {if (filter_presence) filter(., presence == presence_code) else .} %>%
    {if (filter_origin)   filter(., origin %in% origin_code) else .} %>%
    {if (filter_marine)   filter(., marine == marine_code) else .} %>%
    {if (filter_seasonal) filter(., seasonal %in% seasonal_code) else .} %>%
    {if (filter_EX_EW)    filter(., !category %in% c("EW", "EX")) else .}

  filtered_threat <- filtered_all %>%
    {if (filter_category) filter(., category %in% category_code) else .}

  # 5527 mammal species globally (odd), 252 mammal species in Zambia (even)
  # 10936 bird species globally (even), 738 bird species in Zambia (even)
  # 6593 amphibians species globally (odd), 94 amphibians species in Zambia (even)
  # 6358 reptile species globally (even), 37 reptile species in Zambia (odd)
  filtered_small <- filtered_all %>% # use < if odd globally (mam, amp), <= if even (bird, rep)
    {if (filter_range_size) {
      if (odd_n_global) {
        filter(., range_size_quantile < range_threshold)
        } else {
          filter(., range_size_quantile <= range_threshold)
          }
      } else {
        .
      }
      }

  filtered_small_zam <- filtered_all %>% # use < if odd in Zambia (rep), <= if even (mam, bird, amp)
    {if (filter_range_size_zam) {
      if (odd_n_zam) {
        filter(., range_size_quantile_zambia < range_threshold)
      } else {
        filter(., range_size_quantile_zambia <= range_threshold)
      }
    } else {
      .
    }
      }


  filtered_threat_small <- filtered_all %>%
    {if (filter_both) {
      if (odd_n_global) {
        filter(., category %in% category_code | range_size_quantile < range_threshold)
        } else {
          filter(., category %in% category_code | range_size_quantile <= range_threshold)
          }
      } else {
        .
    }
    }

  # if run_reproject is TRUE, then reproject sf object to crs in projection
  if (run_reproject) {
  reprojected_all <- filtered_all %>% st_transform(., projection)
  reprojected_threat <- filtered_threat %>% st_transform(., projection)
  reprojected_small <- filtered_small %>% st_transform(., projection)
  reprojected_small_zam <- filtered_small_zam %>% st_transform(., projection)
  reprojected_threat_small <- filtered_threat_small %>% st_transform(., projection)
}
  # Create template raster by extending the input raster (mask) to the extent of the input polygons.
  mask_filled <- mask
  mask_filled[is.na(mask_filled)] <- 0 # replace NAs with 0s

  template_raster_extent <-
    {if (run_clip) clip_area else input_sf} %>%
    st_transform(., projection) %>%
    {if (class(.)[1] == "sf") as_Spatial(.) else .} %>%
    extent(.)

  template_raster <- extend(mask_filled, template_raster_extent, value = 0)

  # use fasterize to create a raster from prepped input polygons.
  # The raster is created using the function fasterize_fun, which by default 'sum'. (To just produce a 0-1 raster of polygons, use the function 'last'.).
  # If run_mask = TRUE, then the new raster will be cropped and then masked to the input "mask".
  all_richness <-
    {if (run_reproject) reprojected_all else filtered_all} %>%
    fasterize(., template_raster, field = fasterize_field, fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  endemism_richness <-
    {if (run_reproject) reprojected_all else filtered_all} %>%
    fasterize(., template_raster, field = 'inverse_range_glob', fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  endemism_zam_richness <-
    {if (run_reproject) reprojected_all else filtered_all} %>%
    fasterize(., template_raster, field = 'inverse_range_zam', fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  threat_richness <-
    {if (run_reproject) reprojected_threat else filtered_threat} %>%
    fasterize(., template_raster, field = fasterize_field, fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  threat_weighted_richness <-
    {if (run_reproject) reprojected_all else filtered_all} %>%
    fasterize(., template_raster, field = 'threat_weight', fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  small_richness <-
    {if (run_reproject) reprojected_small else filtered_small} %>%
    fasterize(., template_raster, field = fasterize_field, fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  small_zam_richness <-
    {if (run_reproject) reprojected_small_zam else filtered_small_zam} %>%
    fasterize(., template_raster, field = fasterize_field, fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  small_threat_richness <-
    {if (run_reproject) reprojected_threat_small else filtered_threat_small} %>%
    fasterize(., template_raster, field = fasterize_field, fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}


  # if the prep_zambia element is TRUE, then the output raster will be cropped to mask, NAs filled with 0s, and masked to mask. This preps the raster for use in the tradeoff model.
  if (prep_zambia == TRUE) {
    all_richness <- cc_prep_zambia(all_richness)
    endemism_richness <- cc_prep_zambia(endemism_richness)
    endemism_zam_richness <- cc_prep_zambia(endemism_zam_richness)
    threat_richness <- cc_prep_zambia(threat_richness)
    threat_weighted_richness <- cc_prep_zambia(threat_weighted_richness)
    small_richness <- cc_prep_zambia(small_richness)
    small_zam_richness <- cc_prep_zambia(small_zam_richness)
    small_threat_richness <- cc_prep_zambia(small_threat_richness)
  }

  names(all_richness) <- deparse(quote(all_richness))
  names(endemism_richness) <- deparse(quote(endemism_richness))
  names(endemism_zam_richness) <- deparse(quote(endemism_zam_richness))
  names(threat_richness) <- deparse(quote(threat_richness))
  names(threat_weighted_richness) <- deparse(quote(threat_weighted_richness))
  names(small_richness) <- deparse(quote(small_richness))
  names(small_zam_richness) <- deparse(quote(small_zam_richness))
  names(small_threat_richness) <- deparse(quote(small_threat_richness))

  brick <- brick(
    all_richness,
    endemism_richness,
    endemism_zam_richness,
    threat_richness,
    threat_weighted_richness,
    small_richness,
    small_zam_richness,
    small_threat_richness
  )

  # create a list with four elements, and the naming follows this convension: list(name_of_list_item = object, ...). The last element in a function is returned by default.
  list(clipped_sf = clipped_sf,

       filtered_all = filtered_all,
       filtered_threat = filtered_threat,
       filtered_small = filtered_small,
       filtered_small_zam = filtered_small_zam,
       filtered_threat_small = filtered_threat_small,

       reprojected_all = reprojected_all,
       reprojected_threat = reprojected_threat,
       reprojected_small = reprojected_small,
       reprojected_small_zam = reprojected_small_zam,
       reprojected_threat_small = reprojected_threat_small,

       brick = brick
       # all_richness = all_richness,
       # endemism_richness = endemism_richness,
       # endemism_zam_richness = endemism_zam_richness,
       # threat_richness = threat_richness,
       # threat_weighted_richness = threat_weighted_richness,
       # small_richness = small_richness,
       # small_zam_richness = small_zam_richness,
       # small_threat_richness = small_threat_richness
       )
}





# single versions
cc_make_raster <- function (input_sf, run_filter = TRUE,
                             filter_presence = TRUE, presence_code = 1,
                             filter_origin = TRUE, origin_code = c(1, 2),
                             filter_marine = TRUE, marine_code = "False",
                             filter_seasonal = TRUE, seasonal_code = c(1, 2, 3),
                             filter_category = FALSE, category_code = c("CR", "EN", "VU"),
                             filter_range_size = FALSE, range_threshold = 0.5,
                             filter_both = FALSE,
                             run_clip = TRUE, clip_area,
                             run_reproject = TRUE, run_extract = TRUE, run_cast = FALSE,
                             projection = st_crs("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0
                                             +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                             fasterize_field = NULL, fasterize_fun = 'sum',
                             run_mask = TRUE, mask = msk, prep_zambia = FALSE) {
  # if run_filter = FALSE, set all the filter switches to false. This is used to turn off this entire filter step, so that make_raster can be used for just regular polygon-to-raster conversion.
  if (run_filter == FALSE) {
    filter_presence = FALSE
    filter_origin = FALSE
    filter_marine = FALSE
    filter_seasonal = FALSE
    filter_category = FALSE
    filter_range_size = FALSE
    filter_both = FALSE
  }

  # if run_clip is TRUE, then clip the polygons to the specified area, in clip_area, which are unioned using st_union. This should be the same projection as the input sf object.
  clipped_sf <- input_sf %>%
    {if (run_clip)        st_intersection(., st_union(clip_area)) else .} %>%
    {if (run_extract)     st_collection_extract(., "POLYGON") else .} %>%
    {if (run_cast)        st_cast(., "MULTIPOLYGON",) else .}

  # If filter switches are set to TRUE, the line will be run.
  # 1. filter to presence, defaults to code 1 (extant species only). Other unused codes are: 2. Probably Extant (discontinued, ambiguous), 3. Possibly Extant, 4. Possibly Extinct, 5. Extinct, 6. Presence Uncertain.
  # 2. filter to origin, defaults to codes 1 & 2 (native and reintroduced species respectively). Other unused codes are: 3. Introduced, 4. Vagrant, 5. Origin Uncertain, 6. Assisted Colonisation.
  # 3. filter by marine, defaults to "FALSE" (all species but marine ones)
  # 4. filter species with range seasonality. Defaults to c(1,2,3). Codes are: "Resident" (1), "Breeding" (2), "Non-breeding Season" (3), Passage (4), and Seasonal Occurrence Uncertain (5).
  # 5. filter by redlist category. Defaults to c("CR", "EN", "VU")), i.e. "threatened" species. Other categories include LC [Least Concern], EX [Extinct], EW [Extinct in the Wild], NT [Near Threatened], & DD [Data Deficient].
  # 6. filter by total range size. Defaults to selecting small ranged species, i.e. those with ranges smaller than the median. Note that the range_size_quantile was calculated by ordering the total range sizes from smallest to largest and normalizing to 1. Total range size is actually just the sum of the parts of the range we're using (i.e. extant, native, reintroduced, breeding, non-breeding, resident, non-marine, etc.)
  # 7. filter to any species that has a threatened redlist category or is small-ranged (total range size below median).
  filtered_sf <- clipped_sf %>%
    {if (filter_presence) filter(., presence == presence_code) else .} %>%
    {if (filter_origin)   filter(., origin %in% origin_code) else .} %>%
    {if (filter_marine)   filter(., marine == marine_code) else .} %>%
    {if (filter_seasonal) filter(., seasonal %in% seasonal_code) else .} %>%
    {if (filter_category) filter(., category %in% category_code) else .} %>%
    {if (filter_range_size) filter(., range_size_quantile < range_threshold) else .} %>%
    {if (filter_both) filter(., category %in% category_code |
                               range_size_quantile < range_threshold) else .}

  # if run_reproject is TRUE, then reproject sf object to crs in projection
  reprojected_sf <- filtered_sf %>%
    {if (run_reproject)   st_transform(., projection) else .}

  # use fasterize to create a raster from prepped input polygons. Note that the template raster is created by extending the input raster (mask) to the extent of the input polygons. The raster is created using the function fasterize_fun, which by default 'sum' (To just produce a 0-1 raster of polygons, use the function 'last'.). If run_mask = TRUE, then the new raster will be cropped and then masked to the input "mask".
  mask_filled <- mask
  mask_filled[is.na(mask_filled)] <- 0 # replace NAs with 0s
  template_raster_extent <-
    {if (run_clip) clip_area else input_sf} %>%
    st_transform(., projection) %>%
    {if (class(.)[1] == "sf") as_Spatial(.) else .} %>%
    extent(.)

  template_raster <- extend(mask_filled, template_raster_extent, value = 0)

  output_raster <- reprojected_sf %>%
    fasterize(., template_raster, field = fasterize_field, fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  # if the prep_zambia element is TRUE, then the output raster will be cropped to mask, NAs filled with 0s, and masked to mask. This preps the raster for use in the tradeoff model.
  if (prep_zambia == TRUE) {
    output_raster <- crop(output_raster, extent(mask))
    output_raster[is.na(output_raster)] <- 0 # replace NAs with 0s
    output_raster <- raster::mask(output_raster, mask)
  }

  names(output_raster) <- deparse(quote(output_raster))

  # create a list with four elements, and the naming follows this convension: list(name_of_list_item = object, ...). The last element in a function is returned by default.
  list(filtered_sf = filtered_sf,
       clipped_sf = clipped_sf,
       reprojected_sf = reprojected_sf,
       output_raster = output_raster)
}



cc_tradeoff_mod_outputs <- function(mod_toff, input_key, cbetas) {

  # x y coordinates to reconstruct rasters.
  bd_input_dt <- fread(paste0(p_dat, .Platform$file.sep,
                     input_key, .Platform$file.sep,
                     input_key, "-cons-priorities.csv"))

  msk_csv <- fread(paste0(p_dat, .Platform$file.sep,
                          input_key, .Platform$file.sep,
                          input_key, "-mask.csv"))

  bd_input <- bd_input_dt %>%
    cbind(msk_csv, .) %>%
    dt_to_raster(., CRSobj) %>%
    dropLayer(., c(1))

  zambia_no_pas <- mod_toff$inputs$mask[, .(x, y)]
  CRS_obj <- mod_toff$inputs$sp$crs

  # standard model outputs
  bd_input_no_pas <- mod_toff$inputs$cons$cons.priorities %>%
    cbind(zambia_no_pas, .) %>%
    dt_to_raster(., CRSobj)

  conv_r <- dt_to_raster(mod_toff$conv, CRSobj) # raster brick of areas that get converted, to maize and to soy
  # mod_conv_r has two layers:
  # mod_conv_r$maize
  # mod_conv_r$soy
  conv_r_ms <- sum(dt_to_raster(mod_toff$conv, CRSobj) * c(1, 2)) # produces conversion raster with value of 1 for maize, 2 for soy
  conv_r_all <- sum(conv_r) # one raster, with 1 for conversion

  # for visualizations:
  # mod_conv_r_maize[mod_conv_r_maize == 0] <- NA
  # mod_conv_r_soy[mod_conv_r_soy == 0] <- NA # to show only areas converted to soy

  ## create conversion probabilities raster
  convprob_r <- mod_toff %>%
    list(Y = .$inputs$y_std,
         C = .$inputs$carbon_p,
         BD = .$inputs$cons_p,
         COST = .$inputs$cost_p) %>%
    constraints(inlist = ., cbetas = cbetas, silent = FALSE) %>%
    cbind(zambia_no_pas, .) %>%
    dt_to_raster(., CRSobj)

  ## just converted maize areas gradient - NOTE: THESE ARENT WORKING, NEED TO FIX.
  #  but, in the meantime, just use the convprob_r$maize or $soy
  conv_grad_maize <- conv_r$maize * convprob_r$maize
  conv_grad_soy <- conv_r$soy * convprob_r$soy

  list(
    #name = object,
    bd_input = bd_input,
    bd_input_no_pas = bd_input_no_pas,
    conv_r = conv_r,
    conv_r_ms = conv_r_ms,
    conv_r_all = conv_r_all,
    convprob_r = convprob_r,
    conv_grad_maize = conv_grad_maize,
    conv_grad_soy = conv_grad_soy
  )
}

cc_write_bd_to_dt <- function(bd_input, input_key) {
  bd_input_dt <- as.data.table(bd_input, xy = FALSE, keep.rownames=FALSE)
  names(bd_input_dt) <- "cons.priorities"
  fwrite(na.omit(bd_input_dt), # must omit all NA values
         file = paste0(p_dat, .Platform$file.sep,
                       input_key, .Platform$file.sep,
                       input_key, "-cons-priorities.csv"))
}

cc_raster_stats <- function(raster1, raster2) {
# based in part on jaccard by Joona Lehtomäki (https://gist.github.com/jlehtoma/3369793)

  # Calculate the intersection of the two rasters, this is given by adding
  # the binary rasters together -> 2 indicates intersection
  combination <- raster1 + raster2*2
  intersection <- combination == 3

  # Union is all the area covered by the both rasters
  union <- combination >= 1

  overlap_r1 <-
    (cellStats(intersection, stat = "sum") /
       cellStats(raster1, stat = "sum"))
  overlap_r2 <-
    (cellStats(intersection, stat = "sum") /
       cellStats(raster2, stat = "sum"))

  jaccard <-
    (cellStats(intersection, stat = "sum") /
       cellStats(union, stat = "sum"))

  list(
    combination = combination,
    intersection = intersection,
    union = union,
    overlap_r1 = overlap_r1,
    overlap_r2 = overlap_r2,
    jaccard = jaccard
  )
  # list(
  #   overlap_r1 = overlap_r1,
  #   jaccard = jaccard
  # )
}

cc_raster_stats_lite <- function(raster1, raster2) {
  # based in part on jaccard by Joona Lehtomäki (https://gist.github.com/jlehtoma/3369793)

  # Calculate the intersection of the two rasters, this is given by adding
  # the binary rasters together -> 2 indicates intersection
  combination <- raster1 + raster2
  intersection <- combination == 2

  # Union is all the area covered by the both rasters
  union <- combination >= 1

  overlap_r1 <-
    (cellStats(intersection, stat = "sum") /
       cellStats(raster1, stat = "sum"))
  overlap_r2 <-
    (cellStats(intersection, stat = "sum") /
       cellStats(raster2, stat = "sum"))

  jaccard <-
    (cellStats(intersection, stat = "sum") /
       cellStats(union, stat = "sum"))

  list(
    overlap_r1 = overlap_r1,
    overlap_r2 = overlap_r2,
    jaccard = jaccard
  )
}

cc_dt_stats <- function(dt1, dt2) {
  combination <- dt1 + dt2
  intersection <- combination == 2
  union <- combination >= 1

  overlap_r1 <-
    (sum(intersection) /
       sum(dt1))
  overlap_r2 <-
    (sum(intersection) /
       sum(dt2))

  jaccard <-
    (sum(intersection) /
       sum(union))
  list(
    overlap_r1 = overlap_r1,
    overlap_r2 = overlap_r2,
    jaccard = jaccard
  )
}


cc_var <- function(x, y) {
  var_from_control <- sum((x - y)^2)/length(x)
  return(var_from_control)
}

# x <- test
# y <- results_eq$bd_conv_norm[1]
# x - y
# sum((x - y)^2)/33
# results_eq_norm_var[1,]
# length(x)
#
# test_v <- vector(length = 33)
# for (i in 1:33) {
#   test_v[i] <- cc_var(as.numeric(results_eq_norm[i, -c(1, 35:39)]),
#                       results_eq$bd_conv_norm[i])
# }
# test_v
# cc_var(as.numeric(results_eq_norm[1, -c(1, 35:39)]),
#                       results_eq$bd_conv_norm[1])
#
# results_eq_norm[, c(2:34)]
#
#
# transpose(results_eq_norm[1:33, c(2:34)])
# results_eq$bd_conv_norm
# head(results_eq_norm)[, 1:5]
# test_var_ap = apply(results_eq_norm[, c(2:34)], MARGIN = 1, FUN = function(x) sqrt(var(x)/length(x)))
#
# test_mapply <- mapply(FUN = cc_var,
#   transpose(results_eq_norm[1:33, c(2:34)]),  # transposed df, so that each conv result is a column, rather than a row. This is iterated through by column
#   results_eq$bd_conv_norm[1:33] # second vector, with the bd loss in controlling layer
#   )

# from Joona Lehtomäki (https://gist.github.com/jlehtoma/3369793)
jaccard <- function(raster1, raster2, threshhold, warn.uneven=TRUE) {

  # Get the values above the threshhold
  raster1.bin <- raster1 >= threshhold
  raster2.bin <- raster2 >= threshhold

  if (warn.uneven) {
    raster1.size <- count(raster1.bin, 1)
    raster2.size <- count(raster2.bin, 1)
    # Sort from smaller to larger
    sizes <- sort(c(raster1.size, raster2.size))
    if (sizes[2] / sizes[1] > 20) {
      warning("The extents of raster values above the threshhold differ more than 20-fold: Jaccard coefficient may not be informative.")
    }
  }

  # Calculate the intersection of the two rasters, this is given by adding
  # the binary rasters together -> 2 indicates intersection
  combination <- raster1.bin + raster2.bin
  intersection <- combination == 2

  # Union is all the area covered by the both rasters
  union <- combination >= 1

  jaccard <- count(intersection, 1) / count(union, 1)

  return(count(intersection, 1) / count(union, 1))
}




# old ------------------------
cc_make_raster_old <- function (input_sf, run_filter = TRUE,
                             filter_presence = TRUE, presence_code = 1,
                             filter_origin = TRUE, origin_code = c(1, 2),
                             filter_marine = TRUE, marine_code = "False",
                             filter_seasonal = TRUE, seasonal_code = c(1, 2, 3),
                             filter_category = TRUE, category_code = c("CR", "EN", "VU"),
                             filter_range_size = FALSE, range_threshold = 0.5,
                             filter_both = FALSE,
                             run_clip = TRUE, clip_area,
                             run_reproject = TRUE, run_extract = TRUE, run_cast = FALSE,
                             projection = st_crs("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0
                                             +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                             fasterize_field = NULL, fasterize_fun = 'sum', #template_raster = msk,
                             run_mask = FALSE, mask = msk, prep_zambia = TRUE) {
  # if run_filter = FALSE, set all the filter switches to false. This is used to turn off this entire filter step, so that make_raster can be used for just regular polygon-to-raster conversion.
  if (run_filter == FALSE) {
    filter_presence = FALSE
    filter_origin = FALSE
    filter_marine = FALSE
    filter_seasonal = FALSE
    filter_category = FALSE
    filter_range_size = FALSE
    filter_both = FALSE
  }

  # if run_clip is TRUE, then clip the polygons to the specified area, in clip_area, which are unioned using st_union. This should be the same projection as the input sf object.
  clipped_sf <- input_sf %>%
    {if (run_clip)        st_intersection(., st_union(clip_area)) else .} %>%
    {if (run_extract)     st_collection_extract(., "POLYGON") else .} %>%
    {if (run_cast)        st_cast(., "MULTIPOLYGON",) else .}

  # If filter switches are set to TRUE, the line will be run.
  # 1. filter to presence, defaults to code 1 (extant species only). Other unused codes are: 2. Probably Extant (discontinued, ambiguous), 3. Possibly Extant, 4. Possibly Extinct, 5. Extinct, 6. Presence Uncertain.
  # 2. filter to origin, defaults to codes 1 & 2 (native and reintroduced species respectively). Other unused codes are: 3. Introduced, 4. Vagrant, 5. Origin Uncertain, 6. Assisted Colonisation.
  # 3. filter by marine, defaults to "FALSE" (all species but marine ones)
  # 4. filter species with range seasonality. Defaults to c(1,2,3). Codes are: "Resident" (1), "Breeding" (2), "Non-breeding Season" (3), Passage (4), and Seasonal Occurrence Uncertain (5).
  # 5. filter by redlist category. Defaults to c("CR", "EN", "VU")), i.e. "threatened" species. Other categories include LC [Least Concern], EX [Extinct], EW [Extinct in the Wild], NT [Near Threatened], & DD [Data Deficient].
  # 6. filter by total range size. Defaults to selecting small ranged species, i.e. those with ranges smaller than the median. Note that the range_size_quantile was calculated by ordering the total range sizes from smallest to largest and normalizing to 1. Total range size is actually just the sum of the parts of the range we're using (i.e. extant, native, reintroduced, breeding, non-breeding, resident, non-marine, etc.)
  # 7. filter to any species that has a threatened redlist category or is small-ranged (total range size below median).
  filtered_sf <- clipped_sf %>%
    {if (filter_presence) filter(., presence == presence_code) else .} %>%
    {if (filter_origin)   filter(., origin %in% origin_code) else .} %>%
    {if (filter_marine)   filter(., marine == marine_code) else .} %>%
    {if (filter_seasonal) filter(., seasonal %in% seasonal_code) else .} %>%
    {if (filter_category) filter(., category %in% category_code) else .} %>%
    {if (filter_range_size) filter(., range_size_quantile < range_threshold) else .} %>%
    {if (filter_both) filter(., category %in% category_code |
                               range_size_quantile < range_threshold) else .}

  # if run_reproject is TRUE, then reproject sf object to crs in projection
  reprojected_sf <- filtered_sf %>%
    {if (run_reproject)   st_transform(., projection) else .}

  # use fasterize to create a raster from prepped input polygons. Note that the template raster is created by extending the input raster (mask) to the extent of the input polygons. The raster is created using the function fasterize_fun, which by default 'sum' (To just produce a 0-1 raster of polygons, use the function 'last'.). If run_mask = TRUE, then the new raster will be cropped and then masked to the input "mask".
  mask_filled <- mask
  mask_filled[is.na(mask_filled)] <- 0 # replace NAs with 0s
  template_raster <- extend(mask_filled, extent(st_transform(clip_area, projection)), value = 0)

  output_raster <- reprojected_sf %>%
    fasterize(., template_raster, field = fasterize_field, fun = fasterize_fun) %>%
    {if (run_mask) crop(., extent(mask)) else .} %>%
    {if (run_mask) raster::mask(., mask) else .}

  # if the prep_zambia element is TRUE, then the output raster will be cropped to mask, NAs filled with 0s, and masked to mask. This preps the raster for use in the tradeoff model.
  if (prep_zambia == TRUE) {
    output_raster <- crop(output_raster, extent(mask))
    output_raster[is.na(output_raster)] <- 0 # replace NAs with 0s
    output_raster <- raster::mask(output_raster, mask)
  }

  names(output_raster) <- deparse(quote(output_raster))

  # create a list with four elements, and the naming follows this convension: list(name_of_list_item = object, ...). The last element in a function is returned by default.
  list(filtered_sf = filtered_sf,
       clipped_sf = clipped_sf,
       reprojected_sf = reprojected_sf,
       output_raster = output_raster)
}
