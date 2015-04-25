# ----------------------------------------------------------------------
# Unit tests
# ----------------------------------------------------------------------

# unit tests for the model selection functions

source("../../R/ModelSelection.R")

# ----
# AIC weights function
# ----

test_that("AIC weights", {

    library(assertthat)

    # basic results, expected input, example from the article
    AIC_results <- c(204, 202, 206, 206, 214)
    results <- AICw(AIC_results)

    expect_that( results, is_a("list") )
    expect_that( results$AICweights, is_a("numeric") )
    expect_that( results$AICdelta, is_a("numeric") )
    expect_that( length(results), equals(2) )
    expect_that( length(results$AICweights), equals(length(AIC_results)) )
    expect_that( length(results$AICdelta), equals(length(AIC_results)) )
    expect_that( all(results$AICweights > 0), is_true() )
    expect_that( all(results$AICdelta >= 0), is_true() )
    
    print("AIC weights - all tests passed")
})


test_that("AIC weights comparison", {

    library(assertthat)

    # basic results, expected input, example from the article
    AIC_results <- c(204, 202, 206, 206, 214)
    AICweights <- AICw(AIC_results)$AICweights
    results <- compare_AICw(AICweights, NA)

    expect_is( results, "data.frame" )
    expect_true( all(results > 0) )
    expect_equal( min(results), min(AICweights)/max(AICweights) )
    expect_equal( max(results), max(AICweights)/min(AICweights) )
    
    # with model names
    model_names <- c("A1", "A2", "B1", "B2", "C")
    results <- compare_AICw(AICweights, model_names)

    expect_is( results, "data.frame" )
    expect_identical( colnames(results), model_names )
    expect_identical( rownames(results), model_names )
    expect_true( all(results > 0) )
    expect_equal( min(results), min(AICweights)/max(AICweights) )
    expect_equal( max(results), max(AICweights)/min(AICweights) )

    # some inputs that should fail
    expect_error( compare_AICw(NA, NA) )
    expect_error( compare_AICw(matrix(2, 2, 2), NA) )
    expect_error( compare_AICw(c(NA, 1, 1), NA) )
    expect_error( compare_AICw(c(-2, 2), NA) )
    expect_error( compare_AICw(c(2, 2), model_names=c(1, NA)) )
    expect_error( compare_AICw(c(2, 2), matrix(2, 2, 2)) )
    expect_error( compare_AICw(c(2, 2), list(1, NA)) )
    expect_error( compare_AICw(c(2, 2), c("1", "NA", "3")) )

    print("AIC weights comparison - all tests passed")
})
