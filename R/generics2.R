#' @method plot dmrs
#' @export

plot.dmrs <- function(object, type = 1, ...){

if (type == 1){

  return(object[["plot_extreme"]])

} else if (type == 2){

  return(object[["plot_heatmap"]])

} else if (type == 3){

  return(object[["plot_specific"]])

} else {

  stop(paste0("type should be on of 1,2 or 3, but instead it is ", type, "!"))
}

}

