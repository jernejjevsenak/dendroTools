#' @method plot dmrs
#' @export

plot.dmrs <- function(x, ..., type = 2){

if (type == 1){

  plot_out <- x[["plot_extreme"]]

} else if (type == 2){

  plot_out <- x[["plot_heatmap"]]

} else if (type == 3){

  # heatmap of the partial correlations with the control variable
  plot_out <- x[["plot_heatmap_control"]]

  if (is.null(plot_out)) {
    stop(paste0("type = 3 is only available for results of ",
                "daily_response_seascorr() and monthly_response_seascorr() ",
                "with pcor_control = TRUE."))
  }

} else {

  stop(paste0("type should be 1, 2 or 3, but instead it is ", type, "!"))
}

  return(plot_out)

}

