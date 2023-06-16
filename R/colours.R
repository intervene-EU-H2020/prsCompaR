get_discrete_colours <- function(exclude_multiprs=FALSE) {
  discrete_cols <- RColorBrewer::brewer.pal(n = 9,'Set1')[c(1:5,8,9)]

  if (!exclude_multiprs){
    return(c(discrete_cols, 'black'))
  } else {
    return(discrete_cols)
  }
}
