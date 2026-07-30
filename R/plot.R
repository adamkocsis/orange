#' Empty world map (Plate Carée)
#'
#' @param xlim Standard \code{xlim} argument of par in degress longitude.
#' @param ylim Standard \code{ylim} argument of par in degrees latitude.
#' @param grat Graticule resolution. Use \code{grat=NULL}, if you don't want any. 
#' @param col Plot background, RGBA possible.
#' @param grat.col Graticule color RGBA possible. 
#' @param border Border color of rectange on the outside of the plot.
#' @param thick The 0 coordinate circles are highlighted with thicker lines, given here
#' @return The function has no return value.
#' @export
#' @examples
#' data(kentsamples)
#' emptymap(l)
#' points(kentsamples$central_s, pch=16, col=1)
#' points(kentsamples$arctic_m, pch=16, col=2)
#' points(kentsamples$antarctic_l, pch=16, col=3)
#' points(kentsamples$dateline_l, pch=16, col=4)
emptymap <- function(xlim=c(-180, 180), ylim=c(-90, 90), grat=c(30, 30), col="#DDDDFE", grat.col="#FFFFFF", lwd=2, thick=lwd*2, border=grat.col){
	# an empty map
	plot(NULL, xlim=xlim, ylim=ylim, asp=1, xlab="", ylab="", axes=FALSE)
	rect(xleft=xlim[1], xright=xlim[2], ytop=ylim[2], ybottom=ylim[1], col=col, border=border, lwd=lwd)

	if(!is.null(thick)){
		segments(x0=xlim[1], y0=0, x1=xlim[2], y1=0, lwd=thick, col=grat.col)
		segments(x0=0, y0=ylim[1], x1=0, y1=ylim[2], lwd=thick, col=grat.col)

	}

	if(!is.null(grat)){
		longgrat <- seq(xlim[1], xlim[2], grat[1])
		segments(x0=longgrat, y0=ylim[1], x1=longgrat, y1=ylim[2], col=grat.col, lwd=lwd)
		latgrat <- seq(ylim[1], ylim[2], grat[2])
		segments(x0=xlim[1], y0=latgrat, x1=xlim[2], y1=latgrat, col=grat.col, lwd=lwd)
	}
	rect(xleft=xlim[1], xright=xlim[2], ytop=ylim[2], ybottom=ylim[1], col=NA, border=border, lwd=lwd)
}
