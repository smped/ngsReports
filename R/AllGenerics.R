## Set the main Generics
setGeneric("getSummary", function(object){standardGeneric("getSummary")})

setGeneric("getModule",function(object, module){standardGeneric("getModule")})

setGeneric("version", function(object){standardGeneric("version")})

setGeneric("getColours", function(object){standardGeneric("getColours")})

setGeneric("setColours", function(object, PASS, WARN, FAIL, MAX){
    standardGeneric("setColours")
})

setGeneric("setAlpha", function(object, alpha){standardGeneric("setAlpha")})
