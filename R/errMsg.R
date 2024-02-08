## Error message for classes not implemented...
.errNotImp <- function(x){
    cl <- class(x)
    message("Method not implemented for objects of class ", cl)
}