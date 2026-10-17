# High-level, user interface function for the density-based area estimation.

#' Area estimator based on kernel density estimates
#'
#' Apply a geographic kernel density estimator, calculate its confidence region and calculate its area
#'
#' This is a high level function combining \code{kde} and \code{confregion}.
#'
#' @param x Eiher a 2-column numeric matrix with two columns: longitudes and latitudes, or a \code{data.frame} with these columns.
#' @param s Structure to be occupied. Currently only a terra::SpatRaster, default is \code{rast(res=2)}.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param q Confidence level of the density estimation, the level arguments of \code{confregion}.
#' @export
#' @examples
#' # input coordinate
#' data(kentsamples)
#' coords <- kentsamples$dateline_l[1:50, ]
#' 
#' # area estimate with diagnostic plot
#' densityarea(x=coords, q=0.9, kappa=10, plot=TRUE)
#' 
#' # full return for detailed exploration
#' detailed <- densityarea(x=coords, q=0.9, kappa=10, full=TRUE)
#' detailed
#' @rdname densityarea
setGeneric(
	name="densityarea",
	def=function(x,s, ...){
		standardGeneric("densityarea")
	}
)

# default s:  
#' @rdname densityarea
setMethod(
	"densityarea",
	signature=c(x="ANY", s="missing"),
	function(x, ...){
		densityarea(x=x, s=terra::rast(res=2), ...)
	}
)


#' @rdname densityarea
setMethod(
	"densityarea",
	signature=c(x="matrix", s="SpatRaster"),
	function(x, s, q, kappa, long=NULL, lat=NULL, majorheading=NULL, beta=0, full=FALSE, plot=FALSE){

		# input processing
		y <- x
		if(!is.null(long) & !is.null(lat)) y <- x[,c(long, lat)]
		if(ncol(y)!=2) stop ("You must provide a 2-column matrix!")

		# ensure missing values are omitted
		b1 <- !is.na(y[,1])
		b2 <- !is.na(y[,2])

		if(sum(b1)!=sum(b2)) warning("Some coordinate pairs are missing partially, they are omitted!")
		
		# the actual input	
		y <- y[b1 & b2, , drop=FALSE]

		# put this into the kde
		den <- kde(x=y, s=s, kappa=kappa, majorheading=majorheading, beta=beta) 

		# confidence limiting
		conf <-confregion(den, q)

		if(plot){
			terra::plot(den)
			terra::plot(conf$region, add=TRUE, legend=FALSE, col="#FF000033")
			points(y)
		}
		# output wrangling
		if(full){
			res <- list(
				density=den,
				confidence=conf$region,
				estimate=conf$area
			)
			class(res) <- "orange"
			return(res)
		}else{
			return(conf$area)
		}
	}
)

#' @rdname densityarea 
setMethod(
	"densityarea",
	signature=c(x="data.frame", s="SpatRaster"),
	definition=function(x, s, q, kappa, tax=NULL, long="long", lat="lat", majorheading=NULL, beta=0, full=FALSE, plot=FALSE){

		if(!all(c(long, lat) %in% colnames(x))) stop("The 'long' and 'lat' parameters must be valid column names.")
		y <- x[,c(tax, long, lat)]
		
		# taxon - iteration
		if(!is.null(tax)){
			stop("Not yet implemented!")
		# single taxon/set
		}else{
			# delegate to matrix-method
			res <- densityarea(as.matrix(y), s=s, q=q, kappa=kappa, majorheading=majorheading,
				beta=beta, full=full, plot=plot)

			return(res)
		}
	}
)
