#' @method plot dmrs
#' @export

plot.dmrs <- function(x, ..., type = 1){

if (type == 1){

  plot_out <- x[["plot_extreme"]]

} else if (type == 2){

  plot_out <- x[["plot_heatmap"]]

} else if (type == 3){

  plot_out <- x[["plot_specific"]]

} else {

  stop(paste0("type should be on of 1,2 or 3, but instead it is ", type, "!"))
}

  return(plot_out)

}

