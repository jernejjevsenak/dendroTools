#' @method plot dmrs
#' @export

plot.dmrs <- function(x, ..., type = 2){

if (type == 1){

  plot_out <- x[["plot_extreme"]]

} else if (type == 2){

  plot_out <- x[["plot_heatmap"]]

} else {

  stop(paste0("type should be 1 or 2, but instead it is ", type, "!"))
}

  return(plot_out)

}

