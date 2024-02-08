
test_that("fqcVersion gives correct output",{
    expect_equal(fqcVersion(fdl[[1]]), "0.11.2")
    expect_equal(fqcVersion(fdl)$version, rep("0.11.2", length(fdl)))
})

test_that("fqcVersion handles unimplements classes",{
    expect_message(
        fqcVersion(c()), "Method not implemented for objects of class NULL"
    )
})
