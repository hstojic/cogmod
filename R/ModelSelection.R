
# ----------------------------------------------------------------------
# Compute AIC
# ----------------------------------------------------------------------

#' Compute Akaike information criterion of the model. 
#' 
#' This is a standard implementation where log likelihood is multiplied by -2 and penalty term is added in the form of 2 times the number of free parameters in the model. Optionally, finite sample correction can be used. See Burnham and Anderson (2002) for details.  
#'
#' @param logLikelihood Either a log of the maximum likelihood of the model, or a vector of logged probabilities for each trial. Note that it should not be a negative of the log likelihood.
#' @param noPar Number of free parameters in the model.
#' @param correction TRUE if finite sample correction is needed, by default FALSE. Recommended to use if number of observations divided by number of free parameters is less than 40.
#' @return The AIC value of the model in a form of a scalar. 
#' @export
#' @import assertthat
#' @examples
#' # 100 artificial trials with probability for each trial
#' set.seed(1234)
#' maxLikelihood <- runif(100)
#' # you can use either a final log of the maximum likelihood
#' logLikelihood1 <- sum(log(maxLikelihood))
#' # or logged probabilites of each trial
#' logLikelihood2 <- log(maxLikelihood)
#' # computing the AIC value, both give the same value
#' AIC(logLikelihood1, 3)
#' AIC(logLikelihood2, 3)
#' # if finite sample correction is used, then we need to supply number of observations as well
#' AIC(logLikelihood2, 3, length(logLikelihood2))

AIC <- function(logLikelihood, noPar, noObs = NA, correction = FALSE) {
    # basic checks of the inputs
    not_empty(logLikelihood)
    not_empty(noPar)
    is.count(noPar)  # number of par an integer
    assert_that(noPar >= 0)  # non-negative integer
    assert_that(  # a numeric vector
        is.numeric(logLikelihood), 
        !is.matrix(logLikelihood),
        !is.list(logLikelihood),
        !is.data.frame(logLikelihood)
    )
    assert_that(all(logLikelihood <= 0))  # if false, probably not logged
    if (correction) {
        is.count(noObs)  # number of par an integer
        assert_that(noObs >= 0)
    }

    # computing AIC
    if(!correction) {
        AIC <- -2*sum(logLikelihood) + 2*noPar
    } else {
        AIC <- -2*sum(logLikelihood) + 2*noPar + 
                (2*noPar*(noPar+1)) / (noObs - noPar - 1)
    }
    return(AIC)
}


# ----------------------------------------------------------------------
# Compute AIC weights
# ----------------------------------------------------------------------

#' Compute AIC weights for a set of models. 
#' 
#' See Wagenmakers and Farrel (2004) for details.
#'
#' @param AIC_results A numeric vector of AIC values.
#' @return A list of length two. First element is a numeric vector with relative likelihoods and second element is a numeric vector with AIC weights. Vectors are of the same length as inputs, order is preserved. 
#' @export
#' @import assertthat
#' @examples
#' # example from Wagenmakers and Farrel (2004) article
#' AIC_results <- c(204, 202, 206, 206, 214)
#' # computing AIC weights
#' AICw(AIC_results)

AICw <- function(AIC_results) {
    # basic checks of inputs
    not_empty(AIC_results)
    assert_that(  # a numeric vector
        is.numeric(AIC_results), 
        !is.matrix(AIC_results),
        !is.list(AIC_results),
        !is.data.frame(AIC_results)
    )
    assert_that(all(AIC_results > 0)) 

    # first computing the differences between the model with the smallest AIC and all other models
    best_model <- which.min(AIC_results)
    AIC_delta <- AIC_results - AIC_results[best_model]
    
    # second, we compute relative model likelihoods
    rel_likelihood <- exp(-0.5*AIC_delta)
    weights <- rel_likelihood / sum(rel_likelihood)

    # returning results
    return(list(AICdelta = AIC_delta, AICweights = weights))
}


# ----------------------------------------------------------------------
# Compare AIC weights
# ----------------------------------------------------------------------

#' Compare AIC weights for a set of models. 
#' 
#' With AIC weights instead of raw AIC values one can get the sense of magnitude how much more likely (in terms of Kullback-Leibler discrepancy) one model is than some other model that has also been fitted. To get such a likelihood we divide the AIC weight of one model with the AIC weight of the other fitted model. See Wagenmakers and Farrel (2004) for details.
#'
#' @param AICweights A numeric vector of AIC weigths.
#' @param model_names Optional character vector with names of the models.
#' @return A data frame where each cell represents a likelihood of the model denoted by the row name with respect to the model denoted by the cell's column name. Order of the models from the input to the function is preserved. 
#' @export
#' @import assertthat
#' @examples
#' # example from Wagenmakers and Farrel (2004) article
#' AIC <- c(204, 202, 206, 206, 214)
#' AICweights <- AICw(AIC)$AICweights
#' # comparing AIC weights
#' compare_AICw(AICweights)

compare_AICw <- function(AICweights, model_names = NA) {

    # basic checks of inputs
    not_empty(AICweights)
    assert_that(  # a numeric vector
        is.numeric(AICweights), 
        !is.matrix(AICweights),
        !is.list(AICweights),
        !is.data.frame(AICweights)
    )
    assert_that(all(AICweights >= 0)) 
    if (!is.na(model_names)[1]) {
        assert_that(  # a numeric vector
            is.character(model_names), 
            !is.matrix(model_names),
            !is.list(model_names),
            !is.data.frame(model_names),
            length(AICweights) == length(model_names)
        )
    }

    # computing comparisons
    comparison <- as.data.frame(outer(AICweights, 1/AICweights))
    if (!is.na(model_names)[1]) {
        colnames(comparison) <- model_names
        rownames(comparison) <- model_names
    } else {
        colnames(comparison) <- paste(1:length(AICweights))
        rownames(comparison) <- paste(1:length(AICweights))
    }
    return(comparison)
}



# ----------------------------------------------------------------------
# Compute BIC
# ----------------------------------------------------------------------

#' Compute Bayesian information criterion of the model. 
#' 
#' This is a standard implementation where log likelihood is multiplied by -2 and penalty term is added in the form of number of parameters times log of number of observations. See Schwarz (1978) for details.
#'
#' @param logLikelihood Either a natural logarithm of the maximum likelihood of the model, or a vector of logged probabilities for each trial. Note that it should not be a negative of the log likelihood.
#' @param noPar Number of free parameters in the model.
#' @param noObs Number of observations used to obtain the maximum likelihood of the model.
#' @return The BIC value of the model in a form of a scalar. 
#' @export
#' @import assertthat
#' @examples
#' # 100 artificial trials with probability for each trial
#' set.seed(1234)
#' maxLikelihood <- runif(100)
#' # you can use either a final log of the maximum likelihood
#' logLikelihood1 <- sum(log(maxLikelihood))
#' # or logged probabilites of each trial
#' logLikelihood2 <- log(maxLikelihood)
#' # computing the BIC value, both give the same value
#' BIC(logLikelihood1, 3, 100)
#' BIC(logLikelihood2, 3, length(logLikelihood2))

BIC <- function(logLikelihood, noPar, noObs) {
    # basic checks of the inputs
    not_empty(logLikelihood)
    not_empty(noPar)
    not_empty(noObs)
    is.count(noPar)  # number of par an integer
    assert_that(noPar >= 0)  # non-negative integer
    is.count(noObs)  # number of observations an integer
    assert_that(noObs > 0)  # non-negative integer
    assert_that(  # a numeric vector
        is.numeric(logLikelihood), 
        !is.matrix(logLikelihood),
        !is.list(logLikelihood),
        !is.data.frame(logLikelihood)
    )
    assert_that(all(logLikelihood <= 0))  # if false, probably not logged

    # computing BIC
    BIC <- -2*sum(logLikelihood) + noPar*log(noObs)
    return(BIC)
}


# ----------------------------------------------------------------------
# Compute BIC weights
# ----------------------------------------------------------------------

#' Compute BIC weights for a set of models. 
#' 
#' See Wagenmakers and Farrel (2004) for details.
#'
#' @param BIC_results A numeric vector of BIC values.
#' @return A list of length two. First element is a numeric vector with relative likelihoods and second element is a numeric vector with BIC weights. Vectors are of the same length as inputs, order is preserved. 
#' @export
#' @import assertthat
#' @examples
#' # example from Wagenmakers and Farrel (2004) article
#' BIC_results <- c(204, 202, 206, 206, 214)
#' # computing BIC weights
#' BICw(BIC_results)


BICw <- function(BIC_results) {
    
    # basic checks of inputs
    not_empty(BIC_results)
    assert_that(  # a numeric vector
        is.numeric(BIC_results), 
        !is.matrix(BIC_results),
        !is.list(BIC_results),
        !is.data.frame(BIC_results)
    )
    assert_that(all(BIC_results > 0)) 

    # first computing the differences between the model with the smallest BIC and all other models
    best_model <- which.min(BIC_results)
    BIC_delta <- BIC_results - BIC_results[best_model]
    
    # second, we compute relative model likelihoods
    rel_likelihood <- exp(-0.5*BIC_delta)
    weights <- rel_likelihood / sum(rel_likelihood)
    
    # returning results
    return(list(BICdelta = BIC_delta, BICweights = weights))
}



# ----------------------------------------------------------------------
# Compare BIC weights
# ----------------------------------------------------------------------

#' Compare BIC weights for a set of models. 
#' 
#' With BIC weights instead of raw BIC values one can get the sense of magnitude how much more likely (in terms of Kullback-Leibler discrepancy) one model is than some other model that has also been fitted. To get such a likelihood we divide the BIC weight of one model with the BIC weight of the other fitted model. See Wagenmakers and Farrel (2004) for details.
#'
#' @param BICweights A numeric vector of BIC weigths.
#' @param model_names Optional character vector with names of the models.
#' @return A data frame where each cell represents a likelihood of the model denoted by the row name with respect to the model denoted by the cell's column name. Order of the models from the input to the function is preserved. 
#' @export
#' @import assertthat
#' @examples
#' # example from Wagenmakers and Farrel (2004) article
#' BIC <- c(204, 202, 206, 206, 214)
#' BICweights <- BICw(BIC)$BICweights
#' # comparing BIC weights
#' compare_BICw(BICweights)

compare_BICw <- function(BICweights, model_names = NA) {

    # basic checks of inputs
    not_empty(BICweights)
    assert_that(  # a numeric vector
        is.numeric(BICweights), 
        !is.matrix(BICweights),
        !is.list(BICweights),
        !is.data.frame(BICweights)
    )
    assert_that(all(BICweights >= 0)) 
    if (!is.na(model_names)[1]) {
        assert_that(  # a numeric vector
            is.character(model_names), 
            !is.matrix(model_names),
            !is.list(model_names),
            !is.data.frame(model_names),
            length(BICweights) == length(model_names)
        )
    }

    # computing comparisons
    comparison <- as.data.frame(outer(BICweights, 1/BICweights))
    if (!is.na(model_names)[1]) {
        colnames(comparison) <- model_names
        rownames(comparison) <- model_names
    } else {
        colnames(comparison) <- paste(1:length(BICweights))
        rownames(comparison) <- paste(1:length(BICweights))
    }
    return(comparison)
}



# ----------------------------------------------------------------------
# Print a nice (latex) table with modeling results
# ----------------------------------------------------------------------

#' All the information for model selection in one data frame. 
#' 
#' Information important for model selection in one place - number of parameters, log likelihood, AIC, BIC, AIC and BIC delta, AIC and BIC weights.
#'
#' @param model_names A character vector with names of the models.
#' @param logLik A numeric vector with log likelihoods for each model.
#' @param noPar A numeric vector with number of parameters for each model.
#' @param noObs A numeric vector with number of observations for each model.
#' @param latex A logical argument, TRUE if latex formatted table is needed, FALSE by default.
#' @return A data frame with the information on the models that is recommended for the model selection. 
#' @export
#' @import xtable
#' @examples
#' # example from Wagenmakers and Farrel (2004) article
#' logLik <- c(-100, -98, -100, -99, -103)
#' model_names <- c("A1", "A2", "B1", "B2", "C")
#' noPar <- c(2, 3, 3, 4, 4)
#' noObs <- rep(240, 5)
#' MStable(model_names, logLik, noPar, noObs)

MStable <- function(model_names, logLik, noPar, noObs, latex=FALSE) {

    # computing the AIC and BIC results
    AICresults <- BICresults <- rep(NA, length(model_names))
    for (model in 1:length(model_names)) {
        print(model)
        BICresults[model] <- BIC(logLik[model], noPar[model], noObs[model])
        AICresults[model] <- AIC(logLik[model], noPar[model])
    }
    # putting everything in a data frame
    results <- data.frame(Model = model_names,
                          noPar,
                          logLik,
                          AIC  = AICresults,
                          AICd = AICw(AICresults)$AICdelta, 
                          AICw = AICw(AICresults)$AICweights, 
                          BIC  = BICresults,
                          BICd = BICw(BICresults)$BICdelta,
                          BICw = BICw(BICresults)$BICweights
                        )

    # we return either an R object or latex object
    if (latex) {
        colnames(results) <- c("Model", 
                               "$No.Par_i$", 
                               "$log(L_i)$",
                               "$AIC_i$",
                               "$\\Delta_i(AIC)$",
                               "$w_i(AIC)$",
                               "$BIC_i$",
                               "$\\Delta_i(BIC)$",
                               "$w_i(BIC)$")
        return(print(xtable(results), 
                     include.rownames = F, 
                     sanitize.text.function=function(x){x}) )
    } else {
        return(results)
    }
}

