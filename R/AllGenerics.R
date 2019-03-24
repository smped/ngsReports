## Set the main Generics
setGeneric(
    "getSummary",
    function(object){standardGeneric("getSummary")}
)

setGeneric(
    "getModule",
    function(object, module){standardGeneric("getModule")}
)

setGeneric(
    "Total_Deduplicated_Percentage",
    function(object){standardGeneric("Total_Deduplicated_Percentage")}
)
setGeneric("Version", function(object){standardGeneric("Version")})

setGeneric("getColours", function(object){standardGeneric("getColours")})

setGeneric(
    "setColours",
    function(object, PASS, WARN, FAIL, MAX){
        standardGeneric("setColours")
    }
)

setGeneric("setAlpha", function(object, alpha){standardGeneric("setAlpha")})
