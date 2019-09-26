cut_region <- function(Coord, Lval, Rval){
#   select.fg <- (Coord[,1] <= (Lval + 4 * Rval)) &
#     (Coord[,2] <= (2 * Rval))
  select.fg <- rep(TRUE, nrow(Coord))
  return(select.fg)
}