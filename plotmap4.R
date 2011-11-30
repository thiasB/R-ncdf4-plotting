`plotmap4` <-
  function(file, file.small=NULL, sponge=8, sponge.small=15,
           varname = NULL, lsm.file=NULL, lsm.file.small=NULL,
           col=NULL, levels=NULL, sea.col=NULL, timestep=1,
           alt.contour=FALSE, alt.lev=NULL, signif.transp=50,
           signif.contour=FALSE, signif.lev=NULL,file.signif=NULL,
           grid=TRUE, grid.txt=TRUE, grid.lty=4,cities=FALSE,
           i.time=1, i.lev=1, map.lwd=1.0, box.outer=TRUE,
           cex.axis=1, cex.lab=1, cex.main=2, xlim=NULL, ylim=NULL,
           main='', xlab='', ylab='', add=FALSE,
           colourplot=TRUE, hires=TRUE, interior=FALSE,
           alt.poli=TRUE, signif.poli=TRUE,
           nlongrid=10, nlatgrid=5, lon.ind, lat.ind, myunits, myaxis='all',
           projection=NULL)
{
  # load hi-resolution map-data only when hires is set
  # in order to avoid complications on a system without
  # the mapdata package
  if (hires) library('mapdata')
  library('mapproj')

  # read in the data from file
  nc <- nc_open(file)

  # read in the coordinates of the rotated pole

  if( !is.null(nc$dim$rlon)){
    lon    <- nc$dim$rlon$vals
    lat    <- nc$dim$rlat$vals
  }
  else {
    lon    <- nc$dim$lon$vals
    lat    <- nc$dim$lat$vals
  }
  nlon   <- length(lon)
  nlat   <- length(lat)

  if (lon[1] > lon[nlon] & lon[2] > lon[1]){
    lon[lon < 0] <- lon[lon < 0] + 360
  }

  # find the variable name (if not set)
  if (is.null(varname)){
    noread  <- c('lon', 'lat')
    var.n   <- setdiff(names(nc$var), noread)
    if (any(var.n == 'HSURF')){
      varname <- 'HSURF'
    } else {
      dims    <- lapply(nc$var[names(nc$var) %in% var.n], function(x) x$size)
      which.v <- sapply(dims, function(x) all(x[1:2] == c(nlon, nlat)))
      varname <- names(which.v)[which.max(which.v)]
    }
  }

  data    <- ncvar_get(nc, varname)

  if (length(dim(data)) == 3){
    data    <- data[,,i.time]
  } else if (length(dim(data)) == 4){
    data    <- data[,,i.lev,i.time]
  }

  # read in longname and units
  longname    <- nc$var[[varname]]$longname
  if ( myunits == '' ) {
    units       <- nc$var[[varname]]$units
  } else { units <- myunits }
  flag_values <- ncatt_get(nc, nc$var[[varname]],
                           'flag_values')
  if (flag_values$hasatt){
    flag_values <- flag_values$value
  } else {
    flag_values <- NULL
  }
  flag_meanings <- ncatt_get(nc, nc$var[[varname]],
                             'flag_meanings')
  if (flag_meanings$hasatt){
    flag_meanings <- unlist(strsplit(flag_meanings$value, ' '))
  } else {
    flag_meanings <- NULL
  }

  nc_close(nc)

  # if HSURF is plotted, the lsm.file is the same as file
  if (varname == 'HSURF'){
    lsm.file    <- file
    if (!is.null(file.small)){
      lsm.file.small  <- file.small
    }
  }

  # mask out the sponge zone
  if (!is.na(sponge) & sponge > 0) {
    data[c(1:sponge,(nlon-sponge+1):nlon),] <- NA
    data[,c(1:sponge,(nlat-sponge+1):nlat)] <- NA
  }

  # read in the land sea mask
  if (is.null(lsm.file)){
    lsm         <- array(TRUE, dim(data))
    alt.contour <- FALSE
  } else {
    nc.lsm  <- nc_open(lsm.file)
    lsm     <- ncvar_get(nc.lsm, 'FR_LAND')
    lsm     <- lsm > 0.5
    if (alt.contour & any(names(nc.lsm$var) == 'FR_LAND')){
      alt <- ncvar_get(nc.lsm, 'FR_LAND')
    }
    nc_close(nc.lsm)
  }
  # read in the significance mask
  if (is.null(file.signif)){
    signif         <- array(TRUE, dim(data))
    signif.contour <- FALSE
  } else {
    nc.signif  <- nc_open(file.signif)
    signif     <- ncvar_get(nc.signif, 'var1')
    #    signif     <- signif > 0.5
    if (signif.contour & any(names(nc.signif$var) == 'var1')){
      signif <- ncvar_get(nc.signif, 'var1')
    }
    nc_close(nc.signif)
  }

  # read in the data from the nested region
  if (!is.null(file.small)){
    nc.small    <- nc_open(file.small)
    data.small  <- ncvar_get(nc.small, varname)
    lon.small  <- nc.small$dim$rlon$vals
    lat.small  <- nc.small$dim$rlat$vals
    nlon.small  <- length(lon.small)
    nlat.small  <- length(lat.small)
    #    if (lon.small[1] > lon.small(nlon.small) & lon.small[2] > lon.small[1]){
    #      lon.small[lon.small < 0] <- lon.small[lon.small < 0] + 360
    #    }
    nc_close(nc.small)

    if (length(dim(data.small)) == 3){
      data.small    <- data.small[,,i.time]
    } else if (length(dim(data.small)) == 4){
      data.small    <- data.small[,,i.lev,i.time]
    }

    # mask out the sponge zone
    if (!is.na(sponge.small) & sponge.small > 0) {
      data.small[c(1:sponge.small,(nlon.small-sponge.small+1):nlon.small),] <- NA
      data.small[,c(1:sponge.small,(nlat.small-sponge.small+1):nlat.small)] <- NA
    }

    # read in the land sea mask
    if (is.null(lsm.file.small)){
      lsm.small     <- array(TRUE, dim(data.small))
    } else {
      nc.lsm.small  <- nc_open(lsm.file.small)
      lsm.small     <- ncvar_get(nc.lsm.small, 'FR_LAND')
      lsm.small     <- lsm.small > 0.5
      if (alt.contour & any(names(nc.lsm.small$var) == 'FR_LAND')){
        alt.small <- ncvar_get(nc.lsm, 'FR_LAND')
      }
      nc_close(nc.lsm.small)
    }
  }

  # set levels
  if (is.null(levels)){
    if (varname == 'HSURF'){
      levs <- c(-200,0,100,200,500,1000,1500,2000,3000,10000)
    } else {
      if (is.null(flag_values)){
        if (exists('data.small')){
          levs    <- pretty(c(data, data.small), 20)
        } else {
          levs    <- pretty(data, 20)
        }
      } else {
        levs    <- approx(seq(along=flag_values), flag_values,
                          0.5 + 0:length(flag_values), yleft=min(flag_values)-diff(flag_values[1:2]),
                          yright=max(flag_values)+diff(flag_values[length(flag_values) - 0:1]))$y
      }
    }
  } else {
    levs        <- levels
  }

  #  print(paste('plotting variable:',varname))
  # set the colours and levels
  ncols <- length(levs)-1
  if (is.null(col)){
    if (colourplot){
      if (varname == 'HSURF'){
        colours <- .colseq(length(levs)-1, .hsurf2, smooth=0)
        sea.col <- .water
      } else if (varname == 'SOILTYP'){
        colours <- .soil[flag_values+1]
      } else if (varname %in% c('TOT_PREC', 'precip', 'pr')){
        colours <- .colseq(length(levs)-1, .gpcc, smooth=0)
      } else if (varname %in% c('TOT_PREC_DIFF')){
        colours <- rbfuninv(ncols)
      } else if (varname %in% c('NGAUGES')){
        colours <- .colseq(length(levs)-1, .ngauges, smooth=0)
      } else if (varname %in% c('T_2M')){
        colours <- .colseq(length(levs)-1, .cordex_t2m, smooth=0)
      } else if (varname %in% c('T_2M_DIFF')){
        colours <- .colseq(length(levs)-1, .cordex_t2m_diff, smooth=0)
      } else {
        colours     <- rbfun(ncols)
      }
    } else {
      colours <- grey((ncols+1):1/(ncols+1))[2:(ncols+1)]
      sea.col <- 'white'
    }
  } else {
    colours <- rep(col, length.out=ncols)
  }


  if (is.null(xlim)) {
    xlim <- range(lon)
  }
  if (is.null(ylim)) {
    ylim <- range(lat)
  }

  if (hires){
    worlddb <- 'worldHires'
    if (any(lon > 180)) worlddb <- 'worldHires2'
  } else {
    worlddb <- 'world'
    if (any(lon > 180)) worlddb <- 'world2'
  }
  if (alt.poli & interior){
    data(polibound)
    world <- polibound
    if (any(lon > 180)){
      world$x[world$x < 0] <- world$x[world$x < 0] + 180
    }
    for (add.name in c('.*Lake.*', '.*Sea.*')){
      world.add <- try(map(worlddb, region=add.name, plot=FALSE, xlim=xlim, ylim=ylim), silent=TRUE)
      if (class(world.add) != 'try-error' & length(world.add) > 0){
        world   <- list(x=c(world$x, NA, world.add$x),
                        y=c(world$y, NA, world.add$y))
      }
    }

  } else {
    world       <- map(worlddb,interior=interior, plot=FALSE,
                       xlim=xlim, ylim=ylim)
  }
  if (!interior){
    # remove Lesotho and add the Lakes and Seas
    for (add.name in c('.*Lake.*', '.*Sea.*', '.*Island.*')){
      world.add <- try(map(worlddb, region=add.name, plot=FALSE, xlim=xlim, ylim=ylim), silent=TRUE)
      if (class(world.add) != 'try-error' & length(world.add) > 0){
        world   <- list(x=c(world$x, NA, world.add$x),
                        y=c(world$y, NA, world.add$y),
                        names=c(world$names, world.add$names))
      }
    }
    world.remove<- map(worlddb, region='Lesotho', plot=FALSE)
    ind.i       <- which(world$x %in% world.remove$x & world$y %in% world.remove$y)
    world$x[ind.i] <- NA
    world$y[ind.i] <- NA
  }
  dx          <- mean(diff(lon))
  dy          <- mean(diff(lat))
  w.i         <- is.na(world$y) | is.na(world$x) | (!is.na(world$x) & !is.na(world$y) &
                                                    world$x >= min(lon - 0.5*dx) & world$x <= max(lon + 0.5*dx) &
                                                    world$y >= min(lat - 0.5*dx) & world$y <= max(lat + 0.5*dx))
  world$x     <- world$x[w.i]
  world$y     <- world$y[w.i]

  data.tmp      <- data
  data.tmp[!lsm]  <- NA
  image(lon, lat, data.tmp, breaks=levs, axes=FALSE, add=add,
        col=colours, main=main, xlab=xlab, ylab=ylab,
        cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
        xlim=xlim, ylim=ylim)

  if (any(!lsm) & !is.null(sea.col)){
    ## sea points with sea.col
    data.tmp  <- data
    data.tmp[lsm]   <- NA
    image(lon, lat, data.tmp, breaks=c(-1e10,1e10),
          col=sea.col, add=TRUE, axes=FALSE, xlab='', ylab='',
          xlim=xlim, ylim=ylim)
  }

  if (exists('alt')){
    if (is.null(alt.lev)) alt.lev <- pretty(alt, 10)
    contour(lon, lat, alt, lev=alt.lev, drawlabels=FALSE, add=TRUE)
  }

  if (exists('signif')){
    #    print('signif')
    if (is.null(signif.lev)) signif.lev <- pretty(signif, 1)
    #    filled.contour(lon, lat, signif, lev=signif.lev, drawlabels=FALSE, add=TRUE, col='green')
    sigcols <- c(rgb(0,0,0,signif.transp,maxColorValue=255),
                 rgb(255,255,255,0,maxColorValue=255))
    image(lon,lat,signif, col=sigcols,breaks=c(-1,0,1), axes=FALSE,add=TRUE,
          xlim=xlim, ylim=ylim)

  }

  if (!is.null(file.small)){

    # make sponge zone transparent if dev == pdf, otherwise white
    if (!is.null(names(dev.cur())) & names(dev.cur()) == 'pdf'){
      rect(min(lon.small), min(lat.small), max(lon.small), max(lat.small),
           border=1, lwd=3, col=rgb(1,0,0,0.7))
    } else {
      rect(min(lon.small), min(lat.small), max(lon.small), max(lat.small),
           border=1, lwd=1, col='white')
    }

    data.small.tmp  <- data.small

    data.small.tmp[!lsm.small] <- NA
    image(lon.small, lat.small, data.small.tmp, breaks=levs,
          col=colours,  add=TRUE, axes=FALSE, xlab='', ylab='')

    if (any(!lsm.small) & !is.null(sea.col)){
      data.small.tmp          <- data.small
      data.small.tmp[lsm.small]   <- NA
      image(lon.small, lat.small, data.small.tmp, breaks=c(-1e10,1e10),
            col=sea.col, add=TRUE, axes=FALSE, xlab='', ylab='')
    }

    if (exists('alt.small')){
      if (is.null(alt.lev)) alt.lev <- pretty(alt.small, 10)
      contour(lon.small, lat.small, alt.small, lev=alt.lev, drawlabels=FALSE, add=TRUE)
    }
  }

  if (!exists('alt')){
    ##lon-lat-lines
    lines(world$x, world$y, lwd=map.lwd)
  }

  if (box.outer){
    box(lwd=1)
  }

  if (grid){

    if (missing(lon.ind)) {
      lon.ind     <- pretty(lon, nlongrid)
      lon.ind     <- lon.ind[lon.ind > min(lon) & lon.ind < max(lon)]
    }
    if (missing(lat.ind)){
      lat.ind     <- pretty(lat, nlatgrid)
      lat.ind     <- lat.ind[lat.ind > min(lat) & lat.ind < max(lat)]
    }

    abline(h=lat.ind, lty=grid.lty)
    abline(v=lon.ind, lty=grid.lty)

    if (grid.txt){
      lon.ind2 <- lon.ind
      if (any(lon > 180)) lon.ind2[lon.ind > 180] <- lon.ind[lon.ind > 180] - 360
      lon.txt <- paste(abs(lon.ind2), '^o', c('~W', '', '~E')[sign(lon.ind2) + 2])
      #lon.txt <- paste(abs(lon.ind2), '^o', c('', '', '')[sign(lon.ind2) + 2])
      lab.w <- strwidth(parse(text=lon.txt), cex=cex.axis)
      lon.at <- lon.ind
      for (i in (min(which(!is.na(lon.at)))+1):max(which(!is.na(lon.at)))){
        lo.i <- max(which(!is.na(lon.at[1:(i-1)])))
        dist <- lon.at[i] - lon.at[lo.i]
        if (dist < 0.6*(lab.w[i] + lab.w[lo.i])) lon.at[i] <- NA
      }
      if ( myaxis == '' | myaxis == 'all') {
        axis(1, at=lon.at, labels=parse(text=lon.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis)
      } else if (myaxis == 'none' ) NA
      #      axis(1, at=lon.at, labels=parse(text=lon.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis)
      lon.at <- lon.ind
      for (i in (min(which(!is.na(lon.at)))+1):max(which(!is.na(lon.at)))){
        lo.i <- max(which(!is.na(lon.at[1:(i-1)])))
        dist <- lon.at[i] - lon.at[lo.i]
        if (dist < 0.6*(lab.w[i] + lab.w[lo.i])) lon.at[i] <- NA
      }
      if ( myaxis == '' | myaxis == 'all' | myaxis == 'topleft' | myaxis == 'topright' | myaxis == 'topleftright' | myaxis == 'toponly' ) {
        axis(3, at=lon.at, labels=parse(text=lon.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis)
      } else if (myaxis == 'none' ) NA
      #axis(3, at=lon.at, labels=parse(text=lon.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis)

      #      lat.txt <- paste(abs(lat.ind), '^o', c('', '', '')[sign(lat.ind) + 2])
      lat.txt <- paste(abs(lat.ind), '^o', c('~S', '', '~N')[sign(lat.ind) + 2])
      lab.w <- strheight(parse(text=lat.txt), cex=cex.axis)
      lat.at <- lat.ind
      for (i in (min(which(!is.na(lat.at)))+1):max(which(!is.na(lat.at)))){
        lo.i <- max(which(!is.na(lat.at[1:(i-1)])))
        dist <- lat.at[i] - lat.at[lo.i]
        if (dist < (lab.w[i] + lab.w[lo.i])) lat.at[i] <- NA
      }
      if ( myaxis == '' | myaxis == 'all' | myaxis == 'topleft' | myaxis == 'topleftright' ) {
        axis(2, at=lat.at, labels=parse(text=lat.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis, las=1)
      } else if (myaxis == 'none' ) NA
      #axis(2, at=lat.at, labels=parse(text=lat.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis, las=1)
      lat.at <- lat.ind
      for (i in (min(which(!is.na(lat.at)))+1):max(which(!is.na(lat.at)))){
        lo.i <- max(which(!is.na(lat.at[1:(i-1)])))
        dist <- lat.at[i] - lat.at[lo.i]
        if (dist < (lab.w[i] + lab.w[lo.i])) lat.at[i] <- NA
      }
      if ( myaxis == '' | myaxis == 'all' | myaxis == 'topright' | myaxis == 'topleftright' ) {
        axis(4, at=lat.at, labels=parse(text=lat.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis, las=1)
      } else if (myaxis == 'none' ) NA
      #axis(4, at=lat.at, labels=parse(text=lat.txt), tick=FALSE, line=-0.5, cex.axis=cex.axis, las=1)
    }
  }

  ## pollon, pollat
  out <- list(col=colours, lev=levs, sea.col=sea.col, flag_values=flag_values,
              flag_meanings=flag_meanings, longname=longname, units=units)
  class(out) <- 'plotmap'
  invisible(out)
}
