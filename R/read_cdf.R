# TODOs:
# - leap year handling
# - handling/replacing additional_parameters in meta_data
# - check handling of non-default nc's

# utility function to determine time resolution of a netcdf file
get_timestep <- function(file_nc) {
  tunit <- file_nc$dim$time$units
  if (is.null(tunit)) {
    return(NULL)
  }
  time_values <- file_nc$dim$time$vals
  # default
  timestep <- 1

  if (grepl("[[:digit:]]", tunit)) {
    time_cf <- CFtime::CFtime(definition = tunit, 
                              calendar = file_nc$dim$time$calendar,
                              offsets = ncdf4::ncvar_get(nc = file_nc,varid = "time")
    )
    rows <- length(time_cf)
    dates <- matrix(as.integer(unlist(strsplit(
      x = CFtime::as_timestamp(time_cf),split = "-"))),ncol = 3, byrow = T,
      dimnames = list(1:rows, c("year","month","day")))
    
    if (rows > 1) {
      ddiff_all <- dates[2:rows,1:3] - dates[1:(rows-1),1:3]
      ddiff <- apply(X = ddiff_all, MARGIN = c(2),FUN = median)
    } else {
      ddiff <- dates - dates + c(1,0,0) # will always be 1/1/1 and thus evaluated as daily
    }
    if (CFtime::unit(time_cf) == "years") {
      time_res <- "annual"
      first_year <- dates[1,1]
      if (ddiff[,"year"] > 1) timestep <- ddiff[,"year"]
    } else if (CFtime::unit(time_cf) == "months") {
      time_res <- "monthly"
      first_year <- dates[1,1]
      if (ddiff[,"year"] > 0) {
        time_res <- "annual"
      } 
      if (ddiff[,"month"] > 0) {
        time_res <- "monthly"
      }
      if (ddiff[,"month"] > 1) timestep <- ddiff[,"year"]
    } else if (CFtime::unit(time_cf) == "days") {
      first_year <- dates[1,1]
      time_res <- ""
      if (ddiff[,"year"] > 0) {
        time_res <- "annual"
      }
      if (ddiff[,"month"] > 0) {
        time_res <- "monthly"
      }
      if (ddiff[,"day"] > 0) {
        time_res <- "daily"
      }
      if (time_res == "") {
        stop("Automatic detection of firstyear and time resolution failed.")
      }
    } else {
        stop("Automatic detection of firstyear and time resolution failed.")
    }
  } else {
    if (grepl("year", tunit, ignore.case = TRUE)) {
      time_res <- "annual"
      first_year <- time_values[1]
      timestep <- time_values[2] - time_values[1]
    } else {
      stop("Automatic detection of firstyear and time resolution failed.")
    }
  }

  return(
    list(
      time_res = time_res,
      first_year = as.numeric(first_year),
      timestep = timestep
    )
  )
}

# utility function to guess the main variable of a netcdf file
get_main_variable <- function(file_nc) {
  for (variable_name in names(file_nc$var)) {
    
    dim_names <- sapply(file_nc$var[[variable_name]]$dim, function(x) x$name)

    if (grepl("lon", paste(dim_names, collapse = " "), ignore.case = TRUE) &&
        grepl("lat", paste(dim_names, collapse = " "), ignore.case = TRUE)
    ) {
      return(variable_name)
    }
  }
  stop("None of the variables could certainly be identified as main variable.")
}

#' Reads netcdf and returns it header info
#'
#' Reads an arbitrary netcdf + returns the header values to read as lpjml object
#'
#' @param filename netcdf file name
#' @param variable_name optional variable to be read, in case automatic
#' detection does not work as intended or several variables are stored within
#' the file
#'
#' @return header data
#'
read_cdf_meta <- function(
  filename,
  variable_name = NULL
) {
  file_type <- lpjmlkit::detect_io_type(filename = filename)

  file_nc <- ncdf4::nc_open(filename = filename)
  varnames <- names(file_nc$var)
  if (is.null(variable_name)) variable_name <- get_main_variable(
    file_nc = file_nc
  )
  
  ## spatial dimensions
  # get lon/lat information
  total_dim_names <- names(file_nc$dim)
  latdim <- grep("lat", total_dim_names, ignore.case = TRUE)
  londim <- grep("lon", total_dim_names, ignore.case = TRUE)
  lon <- file_nc$dim[[londim]]$vals
  lat <- file_nc$dim[[latdim]]$vals
  nlon <- length(lon)
  nlat <- length(lat)
  
  # now only the main variable dimensions
  dim_names <- sapply(file_nc$var[[variable_name]]$dim, function(x) x$name)
  
  # take the median difference between the cells as the resolution
  if (nlon > 1) {
  spatial_difference_lon <- lon[2:length(lon)] - lon[1:(length(lon) - 1)]
  resolution_lon <- abs(stats::median(spatial_difference_lon))
  }
  if (nlat > 1) {
  spatial_difference_lat <- lat[2:length(lat)] - lat[1:(length(lat) - 1)]
  resolution_lat <- abs(stats::median(spatial_difference_lat))
  }
  # if we have only a subset, we need to make some assumptions for the resolution
  if (nlat == 1 && nlon == 1) {
    resolution_lat <- 0.5
    resolution_lon <- 0.5
  } else if (nlat == 1 && nlon > 1) {
    resolution_lat <- resolution_lon
  } else if (nlat > 1 && nlon == 1) {
    resolution_lon <- resolution_lat
  }
  global_attributes <- ncdf4::ncatt_get(file_nc, 0)
  var_attributes <- ncdf4::ncatt_get(file_nc, variable_name)

  meta_list <- list()
  
  ## time
  time_values <- file_nc$dim$time$vals
  time_info <- get_timestep(file_nc = file_nc)

  if (!is.null(time_info)) {
    if (grepl("[[:digit:]]", file_nc$dim$time$units)) {
      time_cf <- CFtime::CFtime(
                        definition = file_nc$dim$time$units, 
                        calendar = file_nc$dim$time$calendar,
                        offsets = ncdf4::ncvar_get(nc = file_nc,varid = "time") 
                                )
    }
    time_res <- time_info$time_res
    meta_list$firstyear <- time_info$first_year
    meta_list$timestep <- time_info$timestep
    
    if (time_res == "annual") {
      meta_list$nyear <- length(time_values)
      meta_list$nstep <- 1
    } else if (time_res == "monthly") {
      meta_list$nyear <- length(time_values) / 12
      meta_list$nstep <- 12
    } else if (time_res == "daily") {
      # here we don't have to catch tunit == "years/months since X" 
      dates <- CFtime::parse_timestamps(t = time_cf,
                                        x = CFtime::as_timestamp(time_cf))[1:3]
      rows <- dim(dates)[1]
      dates[rows,1] - dates[1,1]
      meta_list$nyear <- dates[rows,1] - dates[1,1] + 1
      meta_list$nstep <- 365 # per definition we can currently only handle this
      
      if (!all.equal(round(meta_list$nyear, 0), meta_list$nyear)) {
        stop(
          "Number of time records not a multiple of 365 - pls check calendar"
        )
      }
    }
    
    meta_list$lastyear <- meta_list$firstyear + 
                                     (meta_list$nyear - 1) * meta_list$timestep

  } else {
    
    if (!variable_name %in% c("cellid", "grid", "LPJGRID")) {
      print("Time information could not be extracted from the netcdf file. Using defaults.")
    }
    meta_list$firstyear <- 0
    meta_list$lastyear <- 0
    meta_list$nyear <- 1
    meta_list$nstep <- 1
    meta_list$timestep <- 1
  }
  
  # get the first of the variable names that has not been identified as the
  # main variable or one of time, lat or lon band names
  lat_name <- dim_names[grep("lat", dim_names, ignore.case = TRUE)]
  lon_name <- dim_names[grep("lon", dim_names, ignore.case = TRUE)]
  time_name <- dim_names[grep("time", dim_names, ignore.case = TRUE)]
  
  dimnames_reduced <- setdiff(dim_names, c(variable_name, lat_name, lon_name, time_name))
  
  if (length(dimnames_reduced) > 0) {
    bands_var_name <- dimnames_reduced[1]
    if (length(dimnames_reduced) > 1) {
      print(paste("Several potential band dimensions detected. Picking:",bands_var_name))
    }
    bands <- ncdf4::ncvar_get(nc = file_nc, varid = bands_var_name)
  } else {
    bands_var_name <- ""
    bands <- 1
  }

  meta_list$sim_name <- global_attributes$title
  meta_list$source <- global_attributes$source
  meta_list$history <- global_attributes$history
  meta_list$name <- tolower(variable_name)
  meta_list$variable <- variable_name
  meta_list$descr <- var_attributes$long_name
  meta_list$unit <- var_attributes$units
  meta_list$nbands <- length(bands)
  meta_list$band_names <- bands

  # Changed for testing
  meta_list$ncell <- nlon * nlat
  # Changed for testing
  meta_list$firstcell <- 0

  meta_list$cellsize_lon <- resolution_lon
  # Can be negative if flipped in cdf -> will be treated afterwards and reset
  meta_list$cellsize_lat <- resolution_lat

  meta_list$format <- file_type
  meta_list$filename <- basename(filename)
  meta_list$subset <- FALSE
  meta_list$datatype <- file_nc$var[[variable_name]]$prec

  meta_list$scalar <- 1
  meta_list$order <- "cellseq"
  meta_list$bigendian <- FALSE

  ncdf4::nc_close(file_nc)

  meta_data <- lpjmlkit::LPJmLMetaData$new(
    x = meta_list,
    data_dir = dirname(filename)
  )
  meta_data
}

#' Reads netcdf and returns it as array
#'
#' Reads an arbitrary netcdf and returns the values at the lon/lat locations..
#'
#' @param filename netcdf file name
#' @param nc_header header data, read in from either meta file or data in
#'        netcdf file
#' @param subset list object defining which subset of the data to be read
#'
#' @return array with netcdf's data, dim=c(nlon,nlat,bands,steps (months/days),years)
#'
read_cdf <- function(
  filename,
  nc_header,
  subset = list()
) {
  variable_name <- nc_header$variable
  file_nc <- ncdf4::nc_open(filename = filename)

  # Determine all years in the file
  years_only <- seq(
    from       = default(nc_header$firstyear, 1901),
    by = default(nc_header$timestep, 1),
    length.out = default(nc_header$nyear, 1)
  )
  # Determine all years in the file
  years_nstep <- rep(
    years_only,
    each = default(nc_header$nstep, 1)
  )

  if (!is.null(names(subset))) {
    if (any(c("lon", "lat", "coords", "coordinates") %in% names(subset))) {
      stop("Subsetting by lon/lat not supported for reading NetCDF files.")
    }

    if (!all(names(subset) %in% c("year", "band"))) {
      stop(
        "Subset must be a list with elements of the array dimensions ",
        "'lon', 'lat', 'year', 'band'."
      )
    }
  }

  # Years to read
  if (!is.null(names(subset)) && "year" %in% names(subset)) {
    if (is.numeric(subset[["year"]])) {
      subset[["year"]] <- years_only[subset[["year"]]]
    }
    timesteps <- which(years_nstep %in% subset[["year"]])
    years <- years_nstep[timesteps]

  } else {
    timesteps <- seq_along(years_nstep)
    years <- years_nstep
  }

  # bands to read
  if (!is.null(names(subset)) && "band" %in% names(subset)) {
    if (nc_header$nbands == 1) stop("Can't extract bands from single band input.")
    if (is.numeric(subset[["band"]])) {
      band_subset_ids <- subset[["band"]]
    } else {
      band_subset_ids <- match(subset[["band"]], nc_header$band_names)
    }
  } else {
    band_subset_ids <- seq_len(nc_header$nbands)
  }
  nbands <- length(band_subset_ids)

  # get lon/lat information
  total_dim_names <- names(file_nc$dim)
  latdim <- grep("lat", total_dim_names, ignore.case = TRUE)
  londim <- grep("lon", total_dim_names, ignore.case = TRUE)
  lon <- file_nc$dim[[londim]]$vals
  lat <- file_nc$dim[[latdim]]$vals
  nlon <- length(lon)
  nlat <- length(lat)
  
  # now only the main variable dimensions
  dim_names <- sapply(file_nc$var[[variable_name]]$dim, function(x) x$name)
  
  # create empty data array to fill in next step
  outdata <- array(
    NA,
    dim = c(
      lon = nlon,
      lat = nlat,
      time = length(timesteps),
      band = nbands
    ),
    dimnames = list(
      lon = lon,
      lat = lat,
      time = create_time_names(
        nstep = default(nc_header$nstep, 1),
        years = unique(years)
      ),
      band = nc_header$band_names[band_subset_ids]
    )
  )
  
  time_idx <- 1
  for (i_time in timesteps) {
    band_idx <- 1
    for (i_band in band_subset_ids) {
      
      if (length(dim_names) == 4) { # lon, lat, time, band
        
        outdata[, , time_idx, band_idx] <- 
          ncdf4::ncvar_get(
            nc = file_nc,
            varid = variable_name,
            count = c(-1, -1, 1, 1),
            start = c(1, 1, i_band, i_time),
            collapse_degen = FALSE # avoid dropping of the lat dim, if len=1
          )[,,1,1] # drop band and time
        
      } else if (length(dim_names) == 3) { # time, lon, lat or bands, lon, lat
        
        if ("time" %in% dim_names) { # lon, lat, time 
          outdata[, , time_idx, band_idx] <- 
            ncdf4::ncvar_get(
              nc = file_nc, varid = variable_name, 
              count = c(-1, -1, 1),
              start = c(1, 1, i_time),
              collapse_degen = FALSE
            )[,,1] # drop time
        } else { # lat, lon, band
          outdata[, , time_idx, band_idx] <- 
            ncdf4::ncvar_get(
              nc = file_nc,
              varid = variable_name,
              count = c(-1, -1, 1),
              start = c(1, 1, i_band),
              collapse_degen = FALSE
            )[,,1] # drop band
        }
        
      }else if (length(dim_names) == 2) {
        
        outdata[, , time_idx, band_idx] <- 
          ncdf4::ncvar_get(
            nc = file_nc, varid = variable_name, 
            count = c(-1, -1),
            start = c(1, 1),
            collapse_degen = FALSE # avoid dropping of the lat dim, if len=1
          ) # don't drop anything
        
      } else {
        
        stop("Less than 2 spatial dimensions in data array information found in netcdf file.")
      }
      band_idx <- band_idx + 1
    }
    time_idx <- time_idx + 1
  }
  # check if the data is in correct (LPJmL) lon and lat order: if not, transpose
  outdata <- transpose_lon_lat(outdata)

  ncdf4::nc_close(file_nc)
  
  # todo: delete leap days here
  return(outdata)
}
