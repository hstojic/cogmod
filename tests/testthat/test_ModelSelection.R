# ----------------------------------------------------------------------
# Unit tests
# ----------------------------------------------------------------------

# unit tests for the model selection functions

#source("../../R/ModelSelection.R")
#library(assertthat)
#library(testthat)

context("Model selection")


# ----
# AIC weights function
# ----

test_that("AIC weights", {

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



# ----
# BIC weights function
# ----

test_that("BIC weights", {

    # basic results, expected input, example from the article
    BIC_results <- c(211, 212.4, 216.4, 219.9, 227.9)
    results <- BICw(BIC_results)

    expect_that( results, is_a("list") )
    expect_that( results$BICweights, is_a("numeric") )
    expect_that( results$BICdelta, is_a("numeric") )
    expect_that( length(results), equals(2) )
    expect_that( length(results$BICweights), equals(length(BIC_results)) )
    expect_that( length(results$BICdelta), equals(length(BIC_results)) )
    expect_that( all(results$BICweights > 0), is_true() )
    expect_that( all(results$BICdelta >= 0), is_true() )
    
    print("BIC weights - all tests passed")
})



test_that("BIC weights comparison", {

    # basic results, expected input, example from the article
    BIC_results <- c(211, 212.4, 216.4, 219.9, 227.9)
    BICweights <- BICw(BIC_results)$BICweights
    results <- compare_BICw(BICweights, NA)

    expect_is( results, "data.frame" )
    expect_true( all(results > 0) )
    expect_equal( min(results), min(BICweights)/max(BICweights) )
    expect_equal( max(results), max(BICweights)/min(BICweights) )
    
    # with model names
    model_names <- c("A1", "A2", "B1", "B2", "C")
    results <- compare_BICw(BICweights, model_names)

    expect_is( results, "data.frame" )
    expect_identical( colnames(results), model_names )
    expect_identical( rownames(results), model_names )
    expect_true( all(results > 0) )
    expect_equal( min(results), min(BICweights)/max(BICweights) )
    expect_equal( max(results), max(BICweights)/min(BICweights) )

    # some inputs that should fail
    expect_error( compare_BICw(NA, NA) )
    expect_error( compare_BICw(matrix(2, 2, 2), NA) )
    expect_error( compare_BICw(c(NA, 1, 1), NA) )
    expect_error( compare_BICw(c(-2, 2), NA) )
    expect_error( compare_BICw(c(2, 2), model_names=c(1, NA)) )
    expect_error( compare_BICw(c(2, 2), matrix(2, 2, 2)) )
    expect_error( compare_BICw(c(2, 2), list(1, NA)) )
    expect_error( compare_BICw(c(2, 2), c("1", "NA", "3")) )

    print("BIC weights comparison - all tests passed")
})



# ----
# Modeling table
# ----

test_that("BIC/AIC modeling table", {

    # basic results, expected input, example from the article
    logLik <- c(-100, -98, -100, -99, -103)
    model_names <- c("A1", "A2", "B1", "B2", "C")
    noPar <- c(2, 3, 3, 4, 4)
    noObs <- rep(240, 5)
    results <- MStable(model_names, logLik, noPar, noObs, FALSE)

    expect_that( results, is_a("data.frame") )
    
    print("BIC/AIC modeling table - all tests passed")
})
