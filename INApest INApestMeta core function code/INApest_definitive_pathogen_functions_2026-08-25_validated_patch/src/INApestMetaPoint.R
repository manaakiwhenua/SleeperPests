###############################################################################
### INApestMetaPoint
### Point-based analogue of INApestMeta for small-area incursion modelling
###
### Each biological point represents exactly ONE individual (e.g. invasive tree)
### or ONE colony/nest (e.g. yellow-legged hornet nest).
###
### Main changes relative to INApestMeta:
###   * no spatial nodes;
###   * no SDD connectivity matrix;
###   * propagules disperse from explicit x/y coordinates using kernels;
###   * establishment can be moderated by a native base-R habitat grid;
###   * active propagules can nudge from a provisional destination to the
###     closest location with the highest better habitat inside a search radius;
###   * SEAM-like information spread is represented by distance-based transfer
###     from persistent information sites.
###
### Coordinates, habitat-grid extents, kernel distances, HabitatSearchRadius,
### KRadius and InfoRadius must use compatible linear units. GIS libraries may be
### used upstream to prepare the base-R habitat grid, but are not required here.
###############################################################################

.ipp_clip01 <- function(x) pmin(1, pmax(0, x))

.ipp_validate_probability_schedule <- function(x,
                                               Ntimesteps,
                                               name,
                                               allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(invisible(TRUE))
    stop(name, " cannot be NULL.")
  }
  if (!is.numeric(x) || is.function(x))
    stop(name, " must be a numeric scalar or numeric vector of length Ntimesteps.")
  if (!(length(x) == 1L || length(x) == Ntimesteps))
    stop(name, " must be a scalar or a numeric vector of length Ntimesteps (", Ntimesteps, ").")
  if (any(!is.finite(x))) stop(name, " must contain only finite values.")
  if (any(x < 0 | x > 1)) stop(name, " values must all be between 0 and 1 inclusive.")
  invisible(TRUE)
}

.ipp_resolve <- function(x, points, timestep, perm, Ntimesteps, name) {
  n <- nrow(points)
  if (n == 0L) return(numeric(0))
  if (is.function(x)) {
    out <- x(points = points, timestep = timestep, perm = perm)
  } else if (length(x) == 1L) {
    out <- rep(as.numeric(x), n)
  } else if (is.numeric(x) && length(x) == Ntimesteps) {
    out <- rep(as.numeric(x[timestep]), n)
  } else {
    stop(name, " must be a scalar, a numeric vector of length Ntimesteps, or a function(points, timestep, perm).")
  }
  if (length(out) == 1L) out <- rep(out, n)
  if (length(out) != n) stop(name, " must return one value per supplied point.")
  as.numeric(out)
}

.ipp_probability <- function(mean_value, sd_value, points, timestep, perm, Ntimesteps, name) {
  mu <- .ipp_resolve(mean_value, points, timestep, perm, Ntimesteps, name)
  if (is.null(sd_value)) return(.ipp_clip01(mu))
  sd <- .ipp_resolve(sd_value, points, timestep, perm, Ntimesteps, paste0(name, "SD"))
  .ipp_clip01(stats::rnorm(length(mu), mean = mu, sd = pmax(0, sd)))
}

###############################################################################
### Spatial-grid and habitat helpers
###############################################################################
INApestSpatialGrid <- function(values, xmin, xmax, ymin, ymax, row_origin = c("ymax", "ymin")) {
  row_origin <- match.arg(row_origin)
  d <- dim(values)
  if (is.null(d) || !(length(d) %in% c(2L, 3L))) stop("values must be a numeric matrix or a 3-D numeric array [row, column, timestep].")
  if (!is.numeric(values)) stop("Spatial-grid values must be numeric.")
  if (d[1] < 1L || d[2] < 1L) stop("Spatial grid must contain at least one row and one column.")
  if (length(d) == 3L && d[3] < 1L) stop("A time-varying spatial grid must contain at least one timestep layer.")
  finite_values <- values[!is.na(values)]
  if (length(finite_values) && any(!is.finite(finite_values))) stop("Spatial-grid values must be finite or NA.")
  if (length(finite_values) && any(finite_values < 0 | finite_values > 1)) stop("Spatial-grid values must all be between 0 and 1 inclusive.")
  extent <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  if (!is.numeric(extent) || length(extent) != 4L || any(!is.finite(extent))) stop("xmin, xmax, ymin and ymax must be finite numeric scalars.")
  if (xmax <= xmin || ymax <= ymin) stop("Spatial-grid extent requires xmax > xmin and ymax > ymin.")
  nr <- d[1]; nc <- d[2]
  structure(list(values = values, xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax), nrow = as.integer(nr), ncol = as.integer(nc), nlayer = if (length(d) == 3L) as.integer(d[3]) else 1L, xres = (as.numeric(xmax)-as.numeric(xmin))/nc, yres = (as.numeric(ymax)-as.numeric(ymin))/nr, row_origin = row_origin), class = "INApestSpatialGrid")
}
print.INApestSpatialGrid <- function(x, ...) { cat("INApestSpatialGrid:", x$ncol, "columns x", x$nrow, "rows x", x$nlayer, "layer(s)\n", "extent:", x$xmin, x$xmax, x$ymin, x$ymax, "\n", "resolution:", x$xres, "x", x$yres, "\n", "row 1 origin:", x$row_origin, "\n"); invisible(x) }
INApestHabitatGrid <- function(values, xmin, xmax, ymin, ymax, row_origin = c("ymax", "ymin")) { out <- INApestSpatialGrid(values,xmin,xmax,ymin,ymax,row_origin); class(out) <- c("INApestHabitatGrid","INApestSpatialGrid"); out }
print.INApestHabitatGrid <- function(x, ...) { cat("INApestHabitatGrid:", x$ncol, "columns x", x$nrow, "rows x", x$nlayer, "layer(s)\n", "extent:", x$xmin, x$xmax, x$ymin, x$ymax, "\n", "resolution:", x$xres, "x", x$yres, "\n", "row 1 origin:", x$row_origin, "\n"); invisible(x) }

INApestSpatialSchedule <- function(..., from = 1, to = NULL, grids = NULL) {
  dots <- list(...)
  if (!is.null(grids)) {
    if (length(dots)>0L) stop("Supply spatial grids either directly in ... or through grids=, not both.")
    if (inherits(grids,"INApestSpatialGrid") || inherits(grids,"SpatRaster")) dots <- list(grids) else if (is.list(grids)) dots <- grids else stop("grids= must be an INApestSpatialGrid, a SpatRaster, or a list of grids.")
  }
  if (length(dots)==1L && is.list(dots[[1]]) && !inherits(dots[[1]],"INApestSpatialGrid") && !inherits(dots[[1]],"SpatRaster")) dots <- dots[[1]]
  if (length(dots)<1L) stop("Supply at least one spatial grid to INApestSpatialSchedule().")
  ngrid <- length(dots)
  if (!is.numeric(from) || length(from)!=ngrid || any(!is.finite(from))) stop("from must contain one finite timestep for each supplied grid.")
  if (any(from<1) || any(from!=floor(from))) stop("from timesteps must be positive integers.")
  from <- as.integer(from)
  if (is.unsorted(from, strictly=TRUE)) stop("from timesteps must be strictly increasing.")
  if (is.null(to)) to <- c(from[-1L]-1L,Inf) else {
    if (!is.numeric(to)||length(to)!=ngrid||any(is.na(to))) stop("to must contain one end timestep for each supplied grid.")
    finite_to <- is.finite(to); if (any(to[finite_to]<1)||any(to[finite_to]!=floor(to[finite_to]))) stop("Finite to timesteps must be positive integers.")
    if (any(to<from)) stop("Each to timestep must be greater than or equal to its from timestep.")
    if (ngrid>1L && any(from[-1L] <= to[-ngrid])) stop("Spatial schedule ranges must not overlap.")
  }
  for (i in seq_len(ngrid)) {
    g <- dots[[i]]
    if (inherits(g,"INApestSpatialGrid")) { if (g$nlayer!=1L) stop("Grid ",i," in INApestSpatialSchedule() has ",g$nlayer," layers. Scheduled grids must be static (one layer); use its own schedule periods or supply the 3-D grid directly instead.") }
    else if (inherits(g,"SpatRaster")) { if (requireNamespace("terra",quietly=TRUE) && terra::nlyr(g)!=1L) stop("Scheduled SpatRaster grids must contain one layer each.") }
    else stop("Every scheduled item must be an INApestSpatialGrid (or optionally a one-layer SpatRaster).")
  }
  native <- vapply(dots,inherits,logical(1),what="INApestSpatialGrid")
  if (sum(native)>1L) { ref <- dots[[which(native)[1L]]]; for (i in which(native)[-1L]) { g <- dots[[i]]; same_geometry <- isTRUE(all.equal(c(ref$xmin,ref$xmax,ref$ymin,ref$ymax,ref$nrow,ref$ncol),c(g$xmin,g$xmax,g$ymin,g$ymax,g$nrow,g$ncol),tolerance=0)) && identical(ref$row_origin,g$row_origin); if(!same_geometry) stop("All INApestSpatialGrid objects in one schedule must have the same extent, dimensions and row_origin.") } }
  labels <- names(dots); if (is.null(labels)) labels <- rep("",ngrid); blank <- is.na(labels)|!nzchar(labels); labels[blank] <- paste0("grid",which(blank)); names(dots) <- labels
  structure(list(grids=dots,from=from,to=as.numeric(to),labels=labels),class="INApestSpatialSchedule")
}
print.INApestSpatialSchedule <- function(x, ...) { cat("INApestSpatialSchedule:\n"); for(i in seq_along(x$grids)){ period <- if(is.infinite(x$to[i])) paste0("timestep ",x$from[i]," onward") else if(x$from[i]==x$to[i]) paste0("timestep ",x$from[i]) else paste0("timesteps ",x$from[i],"-",x$to[i]); cat("  ",x$labels[i],": ",period,"\n",sep="")}; invisible(x) }
.ipp_is_spatial_grid <- function(x) inherits(x,"INApestSpatialGrid")
.ipp_is_spatial_schedule <- function(x) inherits(x,"INApestSpatialSchedule")
.ipp_is_habitat_grid <- function(x) .ipp_is_spatial_grid(x)
.ipp_schedule_index <- function(schedule,timestep,name="SpatialSurface") { if(length(timestep)!=1L||!is.finite(timestep)||timestep<1) stop("timestep must be a positive finite scalar when a spatial schedule is queried."); idx <- which(timestep>=schedule$from & timestep<=schedule$to); if(length(idx)==0L) stop(name," has no spatial grid assigned to timestep ",timestep,". Adjust the from/to values in INApestSpatialSchedule()."); if(length(idx)>1L) stop(name," assigns more than one spatial grid to timestep ",timestep,"."); idx }
.ipp_spatial_layer <- function(SpatialSurface,timestep,name="SpatialSurface") {
  if(is.null(SpatialSurface)||is.function(SpatialSurface)) return(SpatialSurface)
  if(.ipp_is_spatial_schedule(SpatialSurface)){i <- .ipp_schedule_index(SpatialSurface,timestep,name); return(.ipp_spatial_layer(SpatialSurface$grids[[i]],timestep,name))}
  if(.ipp_is_spatial_grid(SpatialSurface)){h<-SpatialSurface;if(h$nlayer==1L)return(h);if(timestep>h$nlayer)stop("A time-varying ",name," grid must have at least Ntimesteps layers.");h$values<-matrix(h$values[,,timestep],nrow=h$nrow,ncol=h$ncol);h$nlayer<-1L;return(h)}
  if(inherits(SpatialSurface,"SpatRaster")){if(!requireNamespace("terra",quietly=TRUE))stop("Package 'terra' is required only when ",name," is supplied as a SpatRaster. Use INApestSpatialGrid() for the native base-R path.");if(terra::nlyr(SpatialSurface)==1L)return(SpatialSurface);if(timestep>terra::nlyr(SpatialSurface))stop("A multi-layer ",name," SpatRaster must have at least Ntimesteps layers.");return(SpatialSurface[[timestep]])}
  stop(name," must be NULL, an INApestSpatialGrid (or INApestHabitatGrid), an INApestSpatialSchedule, a function(x, y, timestep), or optionally a terra SpatRaster.")
}
.ipp_validate_spatial_surface <- function(SpatialSurface,Ntimesteps,name){if(is.null(SpatialSurface))return(invisible(TRUE));if(.ipp_is_spatial_schedule(SpatialSurface)){for(tt in seq_len(Ntimesteps))invisible(.ipp_spatial_layer(SpatialSurface,tt,name));return(invisible(TRUE))};invisible(.ipp_spatial_layer(SpatialSurface,Ntimesteps,name));invisible(TRUE)}
.ipp_grid_cell_centres <- function(h){x<-h$xmin+(seq_len(h$ncol)-0.5)*h$xres;y<-if(identical(h$row_origin,"ymax"))h$ymax-(seq_len(h$nrow)-0.5)*h$yres else h$ymin+(seq_len(h$nrow)-0.5)*h$yres;list(x=x,y=y)}
.ipp_spatial_value <- function(x,y,SpatialSurface,timestep,name="SpatialSurface") {
  if(length(x)==0L)return(numeric(0));if(length(x)!=length(y))stop("x and y must have the same length when a spatial surface is queried.");if(is.null(SpatialSurface))return(rep(1,length(x)));h<-.ipp_spatial_layer(SpatialSurface,timestep,name)
  if(is.function(h)){v<-h(x=x,y=y,timestep=timestep);if(length(v)==1L)v<-rep(v,length(x));if(length(v)!=length(x))stop(name," function must return one value per coordinate pair.")}
  else if(.ipp_is_spatial_grid(h)){v<-rep(0,length(x));inside<-is.finite(x)&is.finite(y)&x>=h$xmin&x<=h$xmax&y>=h$ymin&y<=h$ymax;if(any(inside)){col<-floor((x[inside]-h$xmin)/h$xres)+1L;col<-pmin(h$ncol,pmax(1L,col));row<-if(identical(h$row_origin,"ymax"))floor((h$ymax-y[inside])/h$yres)+1L else floor((y[inside]-h$ymin)/h$yres)+1L;row<-pmin(h$nrow,pmax(1L,row));v[inside]<-h$values[cbind(row,col)]}}
  else {z<-terra::extract(h,cbind(x,y));value_name<-names(h)[1];if(!(value_name%in%names(z)))stop("Could not identify the spatial-value column returned by terra::extract.");v<-as.numeric(z[[value_name]])}
  v[is.na(v)]<-0;if(any(!is.finite(v)))stop(name," returned non-finite values.");if(any(v<0|v>1))stop(name," values must all be between 0 and 1 inclusive.");as.numeric(v)
}
.ipp_probability_spatial <- function(mean_value,sd_value,SpatialSurface,points,timestep,perm,Ntimesteps,name){mu<-.ipp_resolve(mean_value,points,timestep,perm,Ntimesteps,name);spatial<-.ipp_spatial_value(points$x,points$y,SpatialSurface,timestep,paste0(name,"Spatial"));mu<-.ipp_clip01(mu*spatial);if(is.null(sd_value))return(mu);sd<-.ipp_resolve(sd_value,points,timestep,perm,Ntimesteps,paste0(name,"SD"));out<-.ipp_clip01(stats::rnorm(length(mu),mean=mu,sd=pmax(0,sd)));out[spatial<=0]<-0;out}
.ipp_habitat_layer <- function(HabitatSuitability,timestep).ipp_spatial_layer(HabitatSuitability,timestep,"HabitatSuitability")
.ipp_habitat_value <- function(x,y,HabitatSuitability,timestep).ipp_spatial_value(x,y,HabitatSuitability,timestep,"HabitatSuitability")
.ipp_habitat_search <- function(x,y,HabitatSuitability,HabitatSearchRadius,timestep,HabitatSearchCandidates=128L){
  n<-length(x);if(n==0L)return(data.frame(x=numeric(0),y=numeric(0),habitat=numeric(0),habitat_nudged=logical(0)));if(length(y)!=n)stop("x and y must have the same length for habitat search.");current<-.ipp_habitat_value(x,y,HabitatSuitability,timestep);out_x<-x;out_y<-y;out_h<-current;nudged<-rep(FALSE,n);if(is.null(HabitatSuitability)||HabitatSearchRadius<=0)return(data.frame(x=out_x,y=out_y,habitat=out_h,habitat_nudged=nudged));h<-.ipp_habitat_layer(HabitatSuitability,timestep)
  if(.ipp_is_habitat_grid(h)){centres<-.ipp_grid_cell_centres(h);tol<-sqrt(.Machine$double.eps);r2<-HabitatSearchRadius^2;for(i in seq_len(n)){candidate_cols<-which(abs(centres$x-x[i])<=HabitatSearchRadius+tol);candidate_rows<-which(abs(centres$y-y[i])<=HabitatSearchRadius+tol);if(!length(candidate_cols)||!length(candidate_rows))next;rr<-rep(candidate_rows,times=length(candidate_cols));cc<-rep(candidate_cols,each=length(candidate_rows));d2<-(centres$x[cc]-x[i])^2+(centres$y[rr]-y[i])^2;in_radius<-d2<=r2+tol;if(!any(in_radius))next;rr<-rr[in_radius];cc<-cc[in_radius];d2<-d2[in_radius];hv<-h$values[cbind(rr,cc)];hv[is.na(hv)]<-0;hv<-.ipp_clip01(hv);best<-max(hv);if(!is.finite(best)||best<=current[i]+tol)next;best_rows<-which(abs(hv-best)<=tol);j<-best_rows[which.min(d2[best_rows])];out_x[i]<-centres$x[cc[j]];out_y[i]<-centres$y[rr[j]];out_h[i]<-hv[j];nudged[i]<-TRUE}}
  else if(inherits(h,"SpatRaster")){value_name<-names(h)[1];for(i in seq_len(n)){p<-terra::vect(data.frame(x=x[i],y=y[i]),geom=c("x","y"),crs=terra::crs(h));search_area<-terra::buffer(p,width=HabitatSearchRadius);z<-terra::extract(h,search_area,cells=TRUE,xy=TRUE,ID=FALSE);if(is.null(z)||nrow(z)==0L||!(value_name%in%names(z)))next;z<-z[is.finite(z[[value_name]]),,drop=FALSE];if(nrow(z)==0L)next;d2<-(z$x-x[i])^2+(z$y-y[i])^2;keep<-d2<=HabitatSearchRadius^2+sqrt(.Machine$double.eps);z<-z[keep,,drop=FALSE];d2<-d2[keep];if(nrow(z)==0L)next;best<-max(z[[value_name]],na.rm=TRUE);if(!is.finite(best)||best<=current[i]+sqrt(.Machine$double.eps))next;best_rows<-which(abs(z[[value_name]]-best)<=sqrt(.Machine$double.eps));j<-best_rows[which.min(d2[best_rows])];out_x[i]<-z$x[j];out_y[i]<-z$y[j];out_h[i]<-.ipp_clip01(z[[value_name]][j]);nudged[i]<-TRUE}}
  else {HabitatSearchCandidates<-max(1L,as.integer(HabitatSearchCandidates));for(i in seq_len(n)){angle<-stats::runif(HabitatSearchCandidates,0,2*pi);dist<-HabitatSearchRadius*sqrt(stats::runif(HabitatSearchCandidates));cx<-c(x[i],x[i]+dist*cos(angle));cy<-c(y[i],y[i]+dist*sin(angle));hv<-.ipp_habitat_value(cx,cy,h,timestep);best<-max(hv);if(best<=current[i]+sqrt(.Machine$double.eps))next;best_rows<-which(abs(hv-best)<=sqrt(.Machine$double.eps));d2<-(cx[best_rows]-x[i])^2+(cy[best_rows]-y[i])^2;j<-best_rows[which.min(d2)];out_x[i]<-cx[j];out_y[i]<-cy[j];out_h[i]<-hv[j];nudged[i]<-TRUE}}
  data.frame(x=out_x,y=out_y,habitat=out_h,habitat_nudged=nudged)
}

.ipp_draw_displacement <- function(kernel,parents,timestep,perm,name){n<-nrow(parents);if(n==0L)return(data.frame(dx=numeric(0),dy=numeric(0)));if(!is.function(kernel))stop(name," must be a function.");z<-kernel(n=n,parents=parents,timestep=timestep,perm=perm);if(is.matrix(z))z<-as.data.frame(z);if(!is.data.frame(z))stop(name," must return a data.frame or matrix with columns dx and dy.");if(!all(c("dx","dy")%in%names(z)))stop(name," output must contain columns dx and dy.");if(nrow(z)!=n)stop(name," must return one displacement per propagule.");if(any(!is.finite(z$dx))||any(!is.finite(z$dy)))stop(name," returned non-finite displacements.");z[,c("dx","dy"),drop=FALSE]}
.ipp_empty_info_sites <- function() data.frame(info_id=integer(0),source_point_id=integer(0),x=numeric(0),y=numeric(0),created_timestep=integer(0),last_known_timestep=integer(0),active=logical(0),stringsAsFactors=FALSE)
.ipp_add_or_refresh_info_sites <- function(info_sites,locations,timestep){if(nrow(locations)==0L)return(info_sites);if(!"id"%in%names(locations))locations$id<-NA_integer_;for(i in seq_len(nrow(locations))){source_id<-as.integer(locations$id[i]);existing<-integer(0);if(!is.na(source_id)&&nrow(info_sites)>0L)existing<-which(!is.na(info_sites$source_point_id)&info_sites$source_point_id==source_id);if(length(existing)>0L){j<-existing[1];info_sites$x[j]<-locations$x[i];info_sites$y[j]<-locations$y[i];info_sites$last_known_timestep[j]<-timestep;info_sites$active[j]<-TRUE}else{new_id<-if(nrow(info_sites)==0L)1L else max(info_sites$info_id)+1L;info_sites<-rbind(info_sites,data.frame(info_id=new_id,source_point_id=source_id,x=locations$x[i],y=locations$y[i],created_timestep=timestep,last_known_timestep=timestep,active=TRUE,stringsAsFactors=FALSE))}};info_sites}
.ipp_transfer_information <- function(points,info_sites,timestep,InfoRadius,InfoTransferProb,InfoKernel,perm,Ntimesteps){
  if(nrow(points)==0L)return(points);src_sites<-info_sites[info_sites$active,,drop=FALSE];src_points<-points[points$have_info,c("id","x","y"),drop=FALSE];if(nrow(src_points)>0L){names(src_points)[names(src_points)=="id"]<-"source_point_id";src_points$info_id<-NA_integer_;src_points$created_timestep<-NA_integer_;src_points$last_known_timestep<-NA_integer_;src_points$active<-TRUE;if(nrow(src_sites)>0L)src_sites<-src_sites[is.na(src_sites$source_point_id)|!(src_sites$source_point_id%in%src_points$source_point_id),,drop=FALSE]}
  src<-rbind(src_sites[,c("info_id","source_point_id","x","y","created_timestep","last_known_timestep","active"),drop=FALSE],if(nrow(src_points)>0L)src_points[,c("info_id","source_point_id","x","y","created_timestep","last_known_timestep","active"),drop=FALSE]else .ipp_empty_info_sites());if(nrow(src)==0L)return(points);targets<-which(!points$have_info);if(length(targets)==0L)return(points);for(j in targets){d<-sqrt((src$x-points$x[j])^2+(src$y-points$y[j])^2);if(!is.null(InfoKernel)){p_each<-InfoKernel(distance=d,source=src,target=points[j,,drop=FALSE],timestep=timestep,perm=perm);if(length(p_each)==1L)p_each<-rep(p_each,length(d));if(length(p_each)!=length(d))stop("InfoKernel must return one probability per source-target distance.");p_each<-.ipp_clip01(p_each)}else{if(InfoRadius<=0)next;p0<-.ipp_resolve(InfoTransferProb,points[j,,drop=FALSE],timestep,perm,Ntimesteps,"InfoTransferProb");p_each<-ifelse(d<=InfoRadius,p0[1],0)};p_transfer<-1-prod(1-p_each);if(stats::rbinom(1L,1L,p_transfer)==1L)points$have_info[j]<-TRUE};points
}
.ipp_apply_local_capacity <- function(candidates,existing_points,LocalK,KRadius,timestep,perm,Ntimesteps){if(nrow(candidates)==0L)return(logical(0));if(!is.function(LocalK)&&length(LocalK)==1L&&is.infinite(LocalK))return(rep(TRUE,nrow(candidates)));if(KRadius<=0)stop("Finite LocalK requires KRadius > 0 because points do not have an intrinsic node area.");k<-.ipp_resolve(LocalK,candidates,timestep,perm,Ntimesteps,"LocalK");k<-pmax(0,floor(k));keep<-rep(FALSE,nrow(candidates));candidate_order<-sample(seq_len(nrow(candidates)));accepted_xy<-existing_points[,c("x","y"),drop=FALSE];for(ii in candidate_order){if(k[ii]<=0)next;n_local<-0L;if(nrow(accepted_xy)>0L){d2<-(accepted_xy$x-candidates$x[ii])^2+(accepted_xy$y-candidates$y[ii])^2;n_local<-sum(d2<=KRadius^2)};if(n_local<k[ii]){keep[ii]<-TRUE;accepted_xy<-rbind(accepted_xy,candidates[ii,c("x","y"),drop=FALSE])}};keep}
INApestPointKernelExponential <- function(mean_distance){force(mean_distance);function(n,parents=NULL,timestep=NULL,perm=NULL){d<-stats::rexp(n,rate=1/mean_distance);a<-stats::runif(n,0,2*pi);data.frame(dx=d*cos(a),dy=d*sin(a))}}
INApestPointKernelLognormal <- function(meanlog,sdlog){force(meanlog);force(sdlog);function(n,parents=NULL,timestep=NULL,perm=NULL){d<-stats::rlnorm(n,meanlog=meanlog,sdlog=sdlog);a<-stats::runif(n,0,2*pi);data.frame(dx=d*cos(a),dy=d*sin(a))}}
INApestPointKernelFixed <- function(distance){force(distance);function(n,parents=NULL,timestep=NULL,perm=NULL){a<-stats::runif(n,0,2*pi);data.frame(dx=distance*cos(a),dy=distance*sin(a))}}
.ipp_event <- function(perm,timestep,event,point_id,parent_id,x,y,detail,provisional_x=NA_real_,provisional_y=NA_real_,habitat=NA_real_,habitat_nudged=NA){n<-length(point_id);data.frame(perm=rep_len(perm,n),timestep=rep_len(timestep,n),event=rep_len(event,n),point_id=point_id,parent_id=rep_len(parent_id,n),x=x,y=y,detail=rep_len(detail,n),provisional_x=rep_len(provisional_x,n),provisional_y=rep_len(provisional_y,n),habitat=rep_len(habitat,n),habitat_nudged=rep_len(habitat_nudged,n),stringsAsFactors=FALSE)}

###############################################################################
### Main function: INApestMetaPoint
###############################################################################
INApestMetaPoint <- function(
  ModelName="INApestMetaPoint",Nperm,Ntimesteps,InitialPoints,Pathogen=NULL,
  Survival=1,PropaguleProduction,PropaguleEstablishment=1,EnvEstabProb=1,
  SDDkernel,LDDkernel=NULL,LDDrate=0,
  HabitatSuitability=NULL,HabitatSearchRadius=0,HabitatSearchCandidates=128,
  LocalK=Inf,KRadius=0,
  DetectionProb=0,DetectionSD=NULL,DetectionSpatial=NULL,
  ManageProb=0,ManageSD=NULL,ManageSpatial=NULL,
  MortalityProb=0,MortalitySD=NULL,MortalitySpatial=NULL,
  FecundityReduction=0,FecundityReductionSD=NULL,FecundityReductionSpatial=NULL,
  SpreadReduction=0,SpreadReductionSD=NULL,SpreadReductionSpatial=NULL,
  SpreadReductionAppliesTo=c("LDD","all"),
  InitialInfo=NULL,InfoRadius=0,InfoTransferProb=0,InfoKernel=NULL,
  InfoRetentionProb=1,InfoPersistenceSteps=NA,ExternalInfoProb=0,OngoingExternalInfo=FALSE,
  ExternalIncursionGenerator=NULL,OngoingExternalInvasion=FALSE,
  OutputDir=NA,SaveResults=FALSE,DoProgress=TRUE,Seed=NULL
){
  SpreadReductionAppliesTo<-match.arg(SpreadReductionAppliesTo);if(!is.null(Seed))set.seed(Seed)
  if(!is.data.frame(InitialPoints))stop("InitialPoints must be a data.frame.");if(!all(c("x","y")%in%names(InitialPoints)))stop("InitialPoints must contain x and y columns.");if(any(!is.finite(InitialPoints$x))||any(!is.finite(InitialPoints$y)))stop("InitialPoints x and y must be finite.");if(!is.null(Pathogen)&&!inherits(Pathogen,"INApestPointPathogenInteraction"))stop("Pathogen must be NULL or created by INApestPointPathogenInteraction().")
  if(Nperm<1||Nperm!=floor(Nperm))stop("Nperm must be a positive integer.");if(Ntimesteps<1||Ntimesteps!=floor(Ntimesteps))stop("Ntimesteps must be a positive integer.");if(!is.function(SDDkernel))stop("SDDkernel must be a dispersal-kernel function.");if(length(LDDrate)!=1L||!is.finite(LDDrate)||LDDrate<0||LDDrate>1)stop("LDDrate must be a scalar between 0 and 1.");if(LDDrate>0&&is.null(LDDkernel))stop("Supply LDDkernel when LDDrate > 0.");if(HabitatSearchRadius<0)stop("HabitatSearchRadius must be >= 0.")
  .ipp_validate_spatial_surface(HabitatSuitability,Ntimesteps,"HabitatSuitability");.ipp_validate_spatial_surface(DetectionSpatial,Ntimesteps,"DetectionSpatial");.ipp_validate_spatial_surface(ManageSpatial,Ntimesteps,"ManageSpatial");.ipp_validate_spatial_surface(MortalitySpatial,Ntimesteps,"MortalitySpatial");.ipp_validate_spatial_surface(FecundityReductionSpatial,Ntimesteps,"FecundityReductionSpatial");.ipp_validate_spatial_surface(SpreadReductionSpatial,Ntimesteps,"SpreadReductionSpatial")
  if(KRadius<0)stop("KRadius must be >= 0.");if(InfoRadius<0)stop("InfoRadius must be >= 0.");.ipp_validate_probability_schedule(FecundityReduction,Ntimesteps,"FecundityReduction");.ipp_validate_probability_schedule(FecundityReductionSD,Ntimesteps,"FecundityReductionSD",allow_null=TRUE);if(OngoingExternalInvasion&&!is.function(ExternalIncursionGenerator))stop("OngoingExternalInvasion = TRUE requires ExternalIncursionGenerator(timestep, perm).")
  if(is.null(ManageSD)&&!is.function(ManageProb))ManageSD<-ManageProb/10;if(is.null(DetectionSD)&&!is.function(DetectionProb))DetectionSD<-DetectionProb/10;if(is.null(MortalitySD)&&!is.function(MortalityProb))MortalitySD<-MortalityProb/10;if(is.null(FecundityReductionSD))FecundityReductionSD<-FecundityReduction/10;.ipp_validate_probability_schedule(FecundityReductionSD,Ntimesteps,"FecundityReductionSD");if(is.null(SpreadReductionSD)&&!is.function(SpreadReduction))SpreadReductionSD<-(1-SpreadReduction)/10
  snapshots<-list();events<-list();final_points<-list();info_site_results<-list();summaries<-list();pathogen_events<-list();snapshot_i<-0L;event_i<-0L;summary_i<-0L;pathogen_event_i<-0L
  for(perm in seq_len(Nperm)){
    points<-data.frame(x=as.numeric(InitialPoints$x),y=as.numeric(InitialPoints$y),stage=if("stage"%in%names(InitialPoints))as.character(InitialPoints$stage)else"default",stringsAsFactors=FALSE);points$id<-seq_len(nrow(points));points$parent_id<-NA_integer_;points$birth_timestep<-0L;points$have_info<-FALSE;points$detected<-FALSE;points$managing<-FALSE;points$last_known_timestep<-NA_integer_;if(!is.null(Pathogen)){field<-Pathogen$PathogenStateField;if(field%in%names(InitialPoints))points[[field]]<-as.character(InitialPoints[[field]]);points<-Pathogen$Initialize(points,perm=perm)};next_id<-nrow(points)+1L
    if(!is.null(InitialInfo)){if(is.logical(InitialInfo)&&length(InitialInfo)==nrow(points))points$have_info<-InitialInfo else if(is.numeric(InitialInfo)){idx<-intersect(as.integer(InitialInfo),seq_len(nrow(points)));points$have_info[idx]<-TRUE}else stop("InitialInfo must be NULL, a logical vector with one value per InitialPoints row, or row indices.")}else if("have_info"%in%names(InitialPoints)){points$have_info<-as.logical(InitialPoints$have_info);points$have_info[is.na(points$have_info)]<-FALSE}
    info_sites<-.ipp_empty_info_sites()
    if(nrow(points)>0L){initial_detection_prob<-.ipp_probability_spatial(DetectionProb,DetectionSD,DetectionSpatial,points,1L,perm,Ntimesteps,"DetectionProb");initially_detected<-stats::rbinom(nrow(points),1L,initial_detection_prob)==1L;points$detected[initially_detected]<-TRUE;points$have_info[initially_detected]<-TRUE;points$last_known_timestep[initially_detected]<-0L;if(any(initially_detected))info_sites<-.ipp_add_or_refresh_info_sites(info_sites,points[initially_detected,c("id","x","y"),drop=FALSE],timestep=0L)}
    if(!is.null(Pathogen)&&isTRUE(Pathogen$Pathogen$DetectionTriggersInfo)&&nrow(points)>0L){pdet<-Pathogen$Detect(points,1L,perm);points$have_info[pdet]<-TRUE;points$last_known_timestep[pdet]<-0L;if(any(pdet))info_sites<-.ipp_add_or_refresh_info_sites(info_sites,points[pdet,c("id","x","y"),drop=FALSE],timestep=0L)}
    for(timestep in seq_len(Ntimesteps)){
      if(DoProgress)cat("\r","Realisation",perm,"Timestep",timestep,"...");n_start<-nrow(points);n_natural_deaths<-0L;n_management_deaths<-0L;n_propagules<-0L;n_established<-0L;n_external<-0L;n_new_detections<-0L;n_managing<-0L;n_pathogen_deaths<-0L;n_new_infections<-0L;n_pathogen_detections<-0L
      if(nrow(points)>0L){manage_prob<-.ipp_probability_spatial(ManageProb,ManageSD,ManageSpatial,points,timestep,perm,Ntimesteps,"ManageProb");points$managing<-stats::rbinom(nrow(points),1L,manage_prob*as.numeric(points$have_info))==1L;n_managing<-sum(points$managing);survival_prob<-.ipp_clip01(.ipp_resolve(Survival,points,timestep,perm,Ntimesteps,"Survival"));survives_natural<-stats::rbinom(nrow(points),1L,survival_prob)==1L;natural_dead<-!survives_natural;n_natural_deaths<-sum(natural_dead);mortality_prob<-.ipp_probability_spatial(MortalityProb,MortalitySD,MortalitySpatial,points,timestep,perm,Ntimesteps,"MortalityProb");management_dead<-rep(FALSE,nrow(points));eligible<-which(survives_natural&points$managing);if(length(eligible)>0L)management_dead[eligible]<-stats::rbinom(length(eligible),1L,mortality_prob[eligible])==1L;n_management_deaths<-sum(management_dead);dead<-natural_dead|management_dead;if(any(dead)){event_i<-event_i+1L;events[[event_i]]<-.ipp_event(perm,timestep,"death",points$id[dead],points$parent_id[dead],points$x[dead],points$y[dead],ifelse(management_dead[dead],"management_mortality","natural_mortality"))};if(any(management_dead))info_sites<-.ipp_add_or_refresh_info_sites(info_sites,points[management_dead,c("id","x","y"),drop=FALSE],timestep);points<-points[!dead,,drop=FALSE]}
      if(nrow(points)>0L){lambda<-pmax(0,.ipp_resolve(PropaguleProduction,points,timestep,perm,Ntimesteps,"PropaguleProduction"));fecundity_reduction<-.ipp_clip01(.ipp_probability_spatial(FecundityReduction,FecundityReductionSD,FecundityReductionSpatial,points,timestep,perm,Ntimesteps,"FecundityReduction"));lambda<-pmax(0,lambda*(1-fecundity_reduction*as.numeric(points$managing)));n_by_parent<-stats::rpois(nrow(points),lambda);n_propagules<-sum(n_by_parent);if(n_propagules>0L){parent_rows<-rep(seq_len(nrow(points)),n_by_parent);propagule_parents<-points[parent_rows,,drop=FALSE];is_ldd<-stats::rbinom(n_propagules,1L,LDDrate)==1L;spread_reduction<-.ipp_probability_spatial(SpreadReduction,SpreadReductionSD,SpreadReductionSpatial,propagule_parents,timestep,perm,Ntimesteps,"SpreadReduction");if(SpreadReductionAppliesTo=="LDD")suppressed<-is_ldd&propagule_parents$managing&(stats::runif(n_propagules)<spread_reduction)else suppressed<-propagule_parents$managing&(stats::runif(n_propagules)<spread_reduction);keep_propagule<-!suppressed;if(any(keep_propagule)){propagule_parents<-propagule_parents[keep_propagule,,drop=FALSE];is_ldd<-is_ldd[keep_propagule];n_dispersing<-nrow(propagule_parents);dx<-numeric(n_dispersing);dy<-numeric(n_dispersing);sdd_idx<-which(!is_ldd);if(length(sdd_idx)>0L){z<-.ipp_draw_displacement(SDDkernel,propagule_parents[sdd_idx,,drop=FALSE],timestep,perm,"SDDkernel");dx[sdd_idx]<-z$dx;dy[sdd_idx]<-z$dy};ldd_idx<-which(is_ldd);if(length(ldd_idx)>0L){z<-.ipp_draw_displacement(LDDkernel,propagule_parents[ldd_idx,,drop=FALSE],timestep,perm,"LDDkernel");dx[ldd_idx]<-z$dx;dy[ldd_idx]<-z$dy};provisional_x<-propagule_parents$x+dx;provisional_y<-propagule_parents$y+dy;destination<-.ipp_habitat_search(provisional_x,provisional_y,HabitatSuitability,HabitatSearchRadius,timestep,HabitatSearchCandidates);candidates<-data.frame(x=destination$x,y=destination$y,stage=propagule_parents$stage,parent_id=propagule_parents$id,dispersal_type=ifelse(is_ldd,"LDD","SDD"),provisional_x=provisional_x,provisional_y=provisional_y,habitat=destination$habitat,habitat_nudged=destination$habitat_nudged,stringsAsFactors=FALSE);propagule_estab<-.ipp_resolve(PropaguleEstablishment,candidates,timestep,perm,Ntimesteps,"PropaguleEstablishment");env_estab<-.ipp_resolve(EnvEstabProb,candidates,timestep,perm,Ntimesteps,"EnvEstabProb");established<-stats::rbinom(nrow(candidates),1L,.ipp_clip01(propagule_estab*env_estab*candidates$habitat))==1L;candidates<-candidates[established,,drop=FALSE];if(nrow(candidates)>0L){capacity_keep<-.ipp_apply_local_capacity(candidates,points,LocalK,KRadius,timestep,perm,Ntimesteps);candidates<-candidates[capacity_keep,,drop=FALSE]};if(nrow(candidates)>0L){n_established<-nrow(candidates);recruits<-data.frame(x=candidates$x,y=candidates$y,stage=candidates$stage,id=seq.int(next_id,length.out=n_established),parent_id=candidates$parent_id,birth_timestep=timestep,have_info=FALSE,detected=FALSE,managing=FALSE,last_known_timestep=NA_integer_,stringsAsFactors=FALSE);if(!is.null(Pathogen))recruits[[Pathogen$PathogenStateField]]<-"S";next_id<-next_id+n_established;event_i<-event_i+1L;events[[event_i]]<-.ipp_event(perm,timestep,"establishment",recruits$id,recruits$parent_id,recruits$x,recruits$y,candidates$dispersal_type,candidates$provisional_x,candidates$provisional_y,candidates$habitat,candidates$habitat_nudged);points<-rbind(points,recruits)}}}}
      if(OngoingExternalInvasion){ext<-ExternalIncursionGenerator(timestep=timestep,perm=perm);if(!is.null(ext)&&nrow(ext)>0L){if(!is.data.frame(ext)||!all(c("x","y")%in%names(ext)))stop("ExternalIncursionGenerator must return NULL/zero rows or a data.frame with x and y.");if(!"stage"%in%names(ext))ext$stage<-"default";n_external<-nrow(ext);ext_points<-data.frame(x=as.numeric(ext$x),y=as.numeric(ext$y),stage=as.character(ext$stage),id=seq.int(next_id,length.out=n_external),parent_id=NA_integer_,birth_timestep=timestep,have_info=FALSE,detected=FALSE,managing=FALSE,last_known_timestep=NA_integer_,stringsAsFactors=FALSE);if(!is.null(Pathogen)){field<-Pathogen$PathogenStateField;if(field%in%names(ext))ext_points[[field]]<-as.character(ext[[field]])else ext_points[[field]]<-"S"};next_id<-next_id+n_external;points<-rbind(points,ext_points);event_i<-event_i+1L;events[[event_i]]<-.ipp_event(perm,timestep,"external_incursion",ext_points$id,NA_integer_,ext_points$x,ext_points$y,"external")}}
      if(!is.null(Pathogen)&&nrow(points)>0L){pathogen_contacts<-Pathogen$Contact(points=points,home_range=NULL,timestep=timestep,perm=perm,context=list(ModelName=ModelName,phase="pathogen"));pz<-Pathogen$Update(points=points,contacts=pathogen_contacts,home_range=NULL,timestep=timestep,perm=perm,context=list(ModelName=ModelName,phase="pathogen"));points<-pz$points;if(!is.null(pz$events)&&nrow(as.data.frame(pz$events))){pe<-as.data.frame(pz$events,stringsAsFactors=FALSE);pe$perm<-perm;pe$timestep<-timestep;pathogen_event_i<-pathogen_event_i+1L;pathogen_events[[pathogen_event_i]]<-pe[,c("perm","timestep",setdiff(names(pe),c("perm","timestep"))),drop=FALSE];if("new_infection"%in%names(pe))n_new_infections<-sum(pe$new_infection%in%TRUE);if("pathogen_detected"%in%names(pe))n_pathogen_detections<-sum(pe$pathogen_detected%in%TRUE)};if(".pathogen_death"%in%names(points)){dead<-points$.pathogen_death%in%TRUE;n_pathogen_deaths<-sum(dead);if(any(dead)){event_i<-event_i+1L;events[[event_i]]<-.ipp_event(perm,timestep,"death",points$id[dead],points$parent_id[dead],points$x[dead],points$y[dead],"pathogen_mortality")};points<-points[!dead,,drop=FALSE];if(".pathogen_death"%in%names(points))points$.pathogen_death<-NULL}}
      if(nrow(points)>0L){informed<-which(points$have_info);if(length(informed)>0L){persistence<-.ipp_resolve(InfoPersistenceSteps,points[informed,,drop=FALSE],timestep,perm,Ntimesteps,"InfoPersistenceSteps");programmed<-!is.na(persistence);if(any(programmed)){last_known<-points$last_known_timestep[informed];age<-timestep-last_known;stop_info<-programmed&(is.na(last_known)|age>=persistence);if(any(stop_info))points$have_info[informed[stop_info]]<-FALSE};retention_group<-informed[is.na(persistence)&points$have_info[informed]];if(length(retention_group)>0L){retention_prob<-.ipp_clip01(.ipp_resolve(InfoRetentionProb,points[retention_group,,drop=FALSE],timestep,perm,Ntimesteps,"InfoRetentionProb"));retain<-stats::rbinom(length(retention_group),1L,retention_prob)==1L;points$have_info[retention_group[!retain]]<-FALSE}}}
      if(nrow(info_sites)>0L){active<-which(info_sites$active);if(length(active)>0L){site_points<-data.frame(x=info_sites$x[active],y=info_sites$y[active]);persistence<-.ipp_resolve(InfoPersistenceSteps,site_points,timestep,perm,Ntimesteps,"InfoPersistenceSteps");programmed<-!is.na(persistence);if(any(programmed)){age<-timestep-info_sites$last_known_timestep[active];stop_sites<-active[programmed&age>=persistence];if(length(stop_sites)>0L)info_sites$active[stop_sites]<-FALSE};retention_group<-active[is.na(persistence)&info_sites$active[active]];if(length(retention_group)>0L){site_points<-data.frame(x=info_sites$x[retention_group],y=info_sites$y[retention_group]);retention_prob<-.ipp_clip01(.ipp_resolve(InfoRetentionProb,site_points,timestep,perm,Ntimesteps,"InfoRetentionProb"));retain<-stats::rbinom(length(retention_group),1L,retention_prob)==1L;info_sites$active[retention_group[!retain]]<-FALSE}}}
      if(nrow(points)>0L){points<-.ipp_transfer_information(points,info_sites,timestep,InfoRadius,InfoTransferProb,InfoKernel,perm,Ntimesteps);if(OngoingExternalInfo){external_info_prob<-.ipp_clip01(.ipp_resolve(ExternalInfoProb,points,timestep,perm,Ntimesteps,"ExternalInfoProb"));externally_informed<-stats::rbinom(nrow(points),1L,external_info_prob)==1L;points$have_info[externally_informed]<-TRUE}}
      if(nrow(points)>0L){detection_prob<-.ipp_probability_spatial(DetectionProb,DetectionSD,DetectionSpatial,points,timestep,perm,Ntimesteps,"DetectionProb");detected_now<-stats::rbinom(nrow(points),1L,detection_prob)==1L;new_detection<-detected_now&!points$detected;n_new_detections<-sum(new_detection);points$detected[detected_now]<-TRUE;points$have_info[detected_now]<-TRUE;points$last_known_timestep[detected_now]<-timestep;if(any(detected_now))info_sites<-.ipp_add_or_refresh_info_sites(info_sites,points[detected_now,c("id","x","y"),drop=FALSE],timestep);if(any(new_detection)){event_i<-event_i+1L;events[[event_i]]<-.ipp_event(perm,timestep,"detection",points$id[new_detection],points$parent_id[new_detection],points$x[new_detection],points$y[new_detection],"detected")}}
      if(nrow(points)>0L){snapshot_i<-snapshot_i+1L;snapshot<-points;snapshot$perm<-perm;snapshot$timestep<-timestep;snapcols<-c("perm","timestep","id","parent_id","x","y","stage","birth_timestep","have_info","detected","managing","last_known_timestep");if(!is.null(Pathogen))snapcols<-c(snapcols,Pathogen$PathogenStateField);snapshots[[snapshot_i]]<-snapshot[,snapcols,drop=FALSE]}
      summary_i<-summary_i+1L;summaries[[summary_i]]<-data.frame(perm=perm,timestep=timestep,n_start=n_start,n_natural_deaths=n_natural_deaths,n_management_deaths=n_management_deaths,n_propagules=n_propagules,n_established=n_established,n_external=n_external,n_end=nrow(points),n_detected=if(nrow(points)>0L)sum(points$detected)else 0L,n_new_detections=n_new_detections,n_pathogen_detections=n_pathogen_detections,n_new_infections=n_new_infections,n_pathogen_deaths=n_pathogen_deaths,n_infectious=if(!is.null(Pathogen)&&nrow(points)>0L)sum(points[[Pathogen$PathogenStateField]]=="I")else 0L,n_have_info=if(nrow(points)>0L)sum(points$have_info)else 0L,n_managing=n_managing,stringsAsFactors=FALSE)
    }
    fp<-points;fp$perm<-rep(perm,nrow(fp));final_points[[perm]]<-fp;ip<-info_sites;ip$perm<-rep(perm,nrow(ip));info_site_results[[perm]]<-ip
  }
  if(DoProgress)cat("\n");PointHistory<-if(length(snapshots)>0L)do.call(rbind,snapshots)else data.frame();EventLog<-if(length(events)>0L)do.call(rbind,events)else data.frame();FinalPoints<-if(length(final_points)>0L)do.call(rbind,final_points)else data.frame();InfoSites<-if(length(info_site_results)>0L)do.call(rbind,info_site_results)else data.frame();Summary<-if(length(summaries)>0L)do.call(rbind,summaries)else data.frame();PathogenEvents<-if(length(pathogen_events)>0L)do.call(rbind,pathogen_events)else data.frame();out<-list(ModelName=ModelName,PointHistory=PointHistory,EventLog=EventLog,FinalPoints=FinalPoints,InfoSites=InfoSites,Summary=Summary,PathogenEvents=PathogenEvents);class(out)<-c("INApestMetaPoint","INApestPoint","list");if(SaveResults){if(is.na(OutputDir))OutputDir<-"";if(nzchar(OutputDir)&&!dir.exists(OutputDir))dir.create(OutputDir,recursive=TRUE);saveRDS(out,file=file.path(OutputDir,paste0(ModelName,"_PointResults.rds")))};out
}
INApestPoint <- INApestMetaPoint
