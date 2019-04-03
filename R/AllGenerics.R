## Set the main Generics
setGeneric("getSummary", function(object){standardGeneric("getSummary")})

setGeneric("getModule",function(object, module){standardGeneric("getModule")})

setGeneric("fqcVersion", function(object){standardGeneric("fqcVersion")})

setGeneric("fqName", function(object){standardGeneric("fqName")})

setGeneric("getColours", function(object){standardGeneric("getColours")})

setGeneric("setColours", function(object, PASS, WARN, FAIL, MAX){
    standardGeneric("setColours")
})

setGeneric("setAlpha", function(object, alpha){standardGeneric("setAlpha")})