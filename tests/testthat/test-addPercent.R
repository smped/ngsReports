context("Check correct behaviour for .addPercent()")

test_that("Percent sign is added",{
    expect_equal(.addPercent("a"), "a%")
})

test_that("Message is returned by factors",{
    expect_message(.addPercent(as.factor("a")))
})
