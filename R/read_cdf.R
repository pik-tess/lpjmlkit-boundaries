# utility function to determine time resolution of a netcdf file
get_timestep <- function(file_nc) {
  tunit <- file_nc$dim$time$units

  if (is.null(tunit)) {
    return(NULL)
  }

  time_values <- file_nc$dim$time$vals
  # default
  timestep <- 1

  if (grepl("years since", tunit, ignore.case = TRUE)) {
    time_res <- "annual"
    base_year <- as.integer(
      unlist(
        strsplit(
          unlist(
            strsplit(
              tunit,
              split = " ",
              fixed = TRUE
            )
          )[3],
          split = "-",
          fixed = TRUE
        )
      )[1]
    )
    offset_year <- time_values[1]
    first_year <- base_year + offset_year

  } else if (grepl("year", tunit, ignore.case = TRUE)) {
    time_res <- "annual"
    first_year <- time_values[1]
    timestep <- time_values[2] - time_values[1]

  } else if (grepl("days since", tunit, fixed = TRUE)) {
    ddiff <- time_values[2] - time_values[1]
    base_year <- as.integer(
      unlist(
        strsplit(
          unlist(
            strsplit(
              tunit,
              split = " ",
              fixed = TRUE
            )
          )[3],
          split = "-",
          fixed = TRUE
        )
      )[1]
    )

    offset_year <- floor(time_values[1] / 365)
    first_year <- base_year + offset_year

    if (ddiff > 27 && ddiff < 32) {
      time_res <- "monthly"
    } else if (ddiff > 364 && ddiff < 367) {
      time_res <- "annual"
    } else if (ddiff == 1) {
      time_res <- "daily"
    } else {
      stop("Automatic detection of firstyear and time resolution failed.")
    }
  } else {
    stop("Automatic detection of firstyear and time resolution failed.")
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
    ndims <- file_nc$var[[variable_name]]$ndims
    dim_names <- c()
    for (d in 1:ndims) {
      dim_names <- append(dim_names, file_nc$var[[variable_name]]$dim[[d]]$name)
    }
    if (grepl("lon", paste(dim_names, collapse = " "), ignore.case = TRUE) &&
        grepl("lat", paste(dim_names, collapse = " "), ignore.case = TRUE)
    ) {
      return(variable_name)
    }
  }
  message(
    "None of the variables could certainly be identified as main variable,",
    "guessing the last one: ",
    variable_name
  )
  return(variable_name)
}

#' Reads netcdf and returns it header info
#'
#' Reads an arbitrary netcdf + returns the header values to read as lpjml object
#' Todos:
#' - check calendar for irregularities (leapdays)?
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

  # get the first of the variable names that has not been identified as the
  # main variable
  if (length(varnames) > 1) {
    bands_var_name <- varnames[-match(variable_name, varnames)][1]
    bands <- ncdf4::ncvar_get(nc = file_nc, varid = bands_var_name)
  } else {
    bands_var_name <- ""
    bands <- 1
  }

  # get lon/lat information
  dim_names <- names(file_nc$dim)
  latdim <- which(grepl("lat", dim_names, ignore.case = TRUE))
  londim <- which(grepl("lon", dim_names, ignore.case = TRUE))
  lon <- file_nc$dim[[londim]]$vals
  lat <- file_nc$dim[[latdim]]$vals
  nlon <- length(lon)
  nlat <- length(lat)

  # take the median difference between the cells as the resolution
  spatial_difference_lon <- lon[2:length(lon)] - lon[1:(length(lon) - 1)]
  resolution_lon <- stats::median(spatial_difference_lon)
  spatial_difference_lat <- lat[2:length(lat)] - lat[1:(length(lat) - 1)]
  resolution_lat <- stats::median(spatial_difference_lat)

  global_attributes <- ncdf4::ncatt_get(file_nc, 0)
  var_attributes <- ncdf4::ncatt_get(file_nc, variable_name)

  meta_list <- list()

  time_values <- file_nc$dim$time$vals
  time_info <- get_timestep(file_nc = file_nc)

  if (!is.null(time_info)) {
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
      meta_list$nyear <- length(time_values) / 365
      meta_list$nstep <- 365
      if (!all.equal(round(meta_list$nyear, 0), meta_list$nyear)) {
        stop(
          "Number of time records not a multiple of 365 - pls check calendar"
        )
      }
    }
    if (meta_list$timestep == 1) {
      meta_list$lastyear <- meta_list$first_year + meta_list$nyear - 1
    } else {
      meta_list$lastyear <- time_values[length(time_values)]
    }
  } else if (!variable_name %in% c("cellid", "grid", "LPJGRID")) {
    stop("Time information could not be extracted from the netcdf file.")
  } else {
    meta_list$firstyear <- 0
    meta_list$lastyear <- 0
    meta_list$nyear <- 1
    meta_list$nstep <- 1
    meta_list$timestep <- 1
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
  # Can be negative if flipped in cdf
  meta_list$cellsize_lat <- resolution_lat

  meta_list$format <- file_type
  meta_list$filename <- basename(filename)
  meta_list$subset <- FALSE
  meta_list$datatype <- file_nc$var[[variable_name]]$prec

  # Can netcdf be scaled?
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
#' Todos:
#' - subsetting cells/days/months does not work yet
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

  # get lon/lat information - not part of header yet
  dim_names <- names(file_nc$dim)
  latdim <- which(grepl("lat", dim_names, ignore.case = TRUE))
  londim <- which(grepl("lon", dim_names, ignore.case = TRUE))
  lon <- file_nc$dim[[londim]]$vals
  lat <- file_nc$dim[[latdim]]$vals
  nlon <- length(lon)
  nlat <- length(lat)
  time_idx <- 1

  if (length(dim_names) > 3) {
    outdata <- array(
      NA,
      dim = c(
        lon = nlon,
        lat = nlat,
        time = length(timesteps) * default(nc_header$nstep, 1),
        band = nbands
      ),
      dimnames = list(
        lon = lon,
        lat = lat,
        time = create_time_names(
          nstep = default(nc_header$nstep, 1),
          years = years
        ),
        band = nc_header$band_names[band_subset_ids]
      )
    )

    for (i_time in timesteps) {
      if (nc_header$nbands == 1) {
        data <- ncdf4::ncvar_get(
          nc = file_nc, varid = variable_name, count = c(-1, -1, 1),
          start = c(1, 1, i_time)
        )
        outdata[, , time_idx, 1] <- data
      } else {
        data <- ncdf4::ncvar_get(
          nc = file_nc,
          varid = variable_name,
          count = c(-1, -1, -1, 1),
          start = c(1, 1, 1, i_time)
        )
        outdata[, , time_idx, ] <- data[, , band_subset_ids]
      } # end if nbands == 1
      time_idx <- time_idx + 1
    }

  } else if (length(dim_names) == 3) {
    outdata <- array(
      NA,
      dim = c(
        lon = nlon,
        lat = nlat,
        time = length(timesteps) * default(nc_header$nstep, 1),
        band = nbands
      ),
      dimnames = list(
        lon = lon,
        lat = lat,
        time = create_time_names(
          nstep = default(nc_header$nstep, 1),
          years = years
        ),
        band = nc_header$band_names[band_subset_ids]
      )
    )

    for (i_time in timesteps) {
      if (nc_header$nbands == 1) {
        data <- ncdf4::ncvar_get(
          nc = file_nc, varid = variable_name, count = c(-1, -1, 1),
          start = c(1, 1, i_time)
        )
        outdata[, , time_idx, 1] <- data
      } else {
        data <- ncdf4::ncvar_get(
          nc = file_nc,
          varid = variable_name,
          count = c(-1, -1, -1, 1),
          start = c(1, 1, 1, i_time)
        )
        outdata[, , time_idx, ] <- data[, , band_subset_ids]
      } # end if nbands == 1
      time_idx <- time_idx + 1
    }

  }else if (length(dim_names) == 2) {

    outdata <- array(
      NA,
      dim = c(
        lon = nlon,
        lat = nlat,
        band = nbands
      ),
      dimnames = list(
        lon = lon,
        lat = lat,
        band = nc_header$band_names[band_subset_ids]
      )
    )

    if (nc_header$nbands == 1) {
      data <- ncdf4::ncvar_get(
        nc = file_nc, varid = variable_name, count = c(-1, -1),
        start = c(1, 1)
      )
      outdata[, , 1] <- data
    } else {
      data <- ncdf4::ncvar_get(
        nc = file_nc,
        varid = variable_name,
        count = c(-1, -1, 1),
        start = c(1, 1, 1)
      )
      outdata[, , ] <- data[, , band_subset_ids]
    } # end if nbands == 1

  } else {
    stop("No spatial information found in netcdf file.")
  }

  # check if the data is in correct (LPJmL) lon and lat order: if not, transpose
  outdata <- transpose_lon_lat(outdata)

  ncdf4::nc_close(file_nc)

  return(outdata)
}
