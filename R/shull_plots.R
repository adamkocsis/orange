# Plot methods need to be defined

#' Plotting of spherical hulls
#'
#' @rdname plot
#' @param x Spherical hull object.
#' @param add Logical parameter specifying whether the previous plot should be overwritten.
#' @param ... Additional plotting parameters as defined in par.
#' @export
plot.shull <- function(x, ...){
	if(x$type=="centroidcircle") plot_centroidcircle(x,...)


}

plot_centroidcircle <- function(x, add=FALSE, ...){
	circ <- sc_shape(x=x$center, r=x$radius)
	if(!add) plot(circ, col=NA)
	# plot the center
	points(x=x$center[1], y=x$center[2], col="red", pch=3)
	# get the shape
	arcs(circ, ...)
}
