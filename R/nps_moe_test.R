#' Calculate whether NPS score has significantly changed between 2 samples.
#'
#' @description A margin of error is calculated for each sample, from the number of
#' promoters, neutrals (i.e. passives) and detractors. The standard error of
#' their difference is estimated using the Pythagorean formula, and the absolute
#' difference of the two samples is compared to this multiplied by the critical
#' value (aka z*-value).
#'
#' The return value is in (-1, 0, +1), according to whether a significant decrease
#' is found, no significant change, or a significant increase, respectively. If
#' the total for a sample is 0, then 0 is returned.
#'
#' Formula is based on the one found in this blog post:
#' (https://www.genroe.com/blog/how-to-calculate-margin-of-error-and-other-stats-
#' for-nps/5994).
#'
#' @param p_0 Number of Promoters in latest sample
#' @param n_0 Number of Neutrals in latest sample
#' @param d_0 Number of Detractors in latest sample
#' @param p_1 Number of Promoters in oldest sample
#' @param n_1 Number of Neutrals in oldest sample
#' @param d_1 Number of Detractors in oldest sample
#' @param z_val Critical value multiplier; 1.96 by default for a 95% confidence
#' interval. See [this table]
#' (http://www.ltcconline.net/greenl/courses/201/estimation/smallConfLevelTable.htm)
#' for further values of z_val for common confidence intervals.
#'
#' @return A value in (-1, 0, +1); see notes above.
#'
#' @examples
#' # Test with a 99% confidence interval
#' \dontrun{
#' nps_moe_test(123, 456, 789, 321, 654, 987, z_val = 2.58)
#' }

nps_moe_test <- function(p_0, n_0, d_0,
                         p_1, n_1, d_1,
                         z_val = 1.96) {
  if (NA %in% c(p_0, n_0, d_0, p_1, n_1, d_1)) {
    return(0)
  }

  t_0 <- p_0 + n_0 + d_0
  if (t_0 == 0) {
    return(0)
  }
  nps_0 <- (p_0 - d_0) / t_0
  t_1 <- p_1 + n_1 + d_1
  if (t_1 == 0) {
    return(0)
  }
  nps_1 <- (p_1 - d_1) / t_1

  var_0 <- ((1 - nps_0)^2 * p_0 + nps_0^2 * n_0 + (-1 - nps_0)^2 * d_0) / t_0
  var_1 <- ((1 - nps_1)^2 * p_1 + nps_1^2 * n_1 + (-1 - nps_1)^2 * d_1) / t_1

  se_0 <- sqrt(var_0 / t_0)
  se_1 <- sqrt(var_1 / t_1)

  if (abs(nps_0 - nps_1) > z_val * sqrt(se_0^2 + se_1^2)) {
    if (nps_0 > nps_1) {
      return(1)
    }
    return(-1)
  } else {
    0
  }
}
