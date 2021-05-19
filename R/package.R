#' Length selected growth models.
#'
#' Fit growth models to length-at-age fisheries data, correcting for
#' gear selectivity, length-stratified sampling, and ageing errors.
#'
#' Provides functionality to fit von Bertalanffy, Gompertz, logistic,
#' or Schnute-Richards growth models to length-at-age fisheries data
#' by maximum likelihood. The package uses the method of Candy (2007)
#' extended to account for ageing error, alongside fishery selectivity
#' and length-stratified initial sampling.
#'
#' Data required are based on the common practice of length stratified
#' sampling fish that are randomly sampled from a sample of fish
#' landed in a fishing operation. Measured fish selected from a
#' fishery stock based on length-stratified initial sampling
#' (i.e. random or otherwise selection of particular fish that are
#' measured, and contribute to the pool of fish that have the
#' potential to be aged). Aged fish are selected randomly for ageing
#' from the length-stratified sample. Selectivity is based on
#' fishery-specific selectivity.
#'
#' The ageing error matrix is calculated prior to the use of lgsm (via Punt et
#' al. (XXXX) or similar methods). The rows of this matrix correspond to
#' 'actual' age and the columns the 'measured' age, and each column is a
#' vector of conditional probabilities that an individual is of an actual' age
#' given their 'measured' age.
#'
#'
#' @name lsgm-package
#' @docType package
#' @author Genevieve Phillips
NULL
