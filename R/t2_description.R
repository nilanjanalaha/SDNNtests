#' Artificial data resembling the aggregated IgG binding immune
#' response data from HVTN 100 (Laha et al (2020)).
#'
#' The aggregated IgG binding immune response from the HVTN 100 trial has a skewed
#' symmetric histogram (see Laha et al., 2020). Therefore, we use the Gamma distribution to generate an
#' artificial dataset  resembling this variable. This dataset  has only
#' one array t2 of size 180 (the sample size of the aggregated  HVTN 100 immune-response data). t2  is sampled from
#'  Gamma distribution with shape parameter 3.6083 and rate parameter 0.7058, which are the parametric MLEs
#'  of the shape and the rate when a gamma distribution is fitted to  the
#' aggregated immune response data from HVTN 100.
#'
#'
#' @docType data
#'
#' @usage SDNNtests::t2
#'
#' @format An object of class \code{"numeric"}.
#'
#' @keywords datasets
#'
#' @references Laha, N., Moodie, Z., Huang, Y., and Luedtke, A. (2020).
#' \emph{ Improved inference for vaccine-induced immune responses
#'        via shape-constrained methods}. Submitted.
#'
#' @examples
#' t2
#' hist(t2); boxplot(t2)
"t2"
