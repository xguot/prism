#' Route Longitudinal Missing Data
#'
#' @param data A data frame.
#' @param threshold_n Sample size threshold.
#' @param threshold_miss Missingness rate threshold.
#' @param threshold_skew Skewness threshold (p-value for Shapiro-Wilk).
#'
#' @return A character string indicating the selected method.
#' @export
route_missing_data <- function(data, 
                               threshold_n = 100, 
                               threshold_miss = 0.2, 
                               threshold_skew = 0.05) {
        n <- nrow(data)
        miss_rate <- mean(is.na(data))
        
        # Check normality for the first non-missing column as a proxy for skew
        # In a real scenario, we'd check all relevant columns.
        skew_p <- shapiro.test(na.omit(data[[1]]))$p.value

        if (n < threshold_n || miss_rate > threshold_miss || skew_p < threshold_skew) {
                return("FIML/RMB")
        }

        return("Lagrange-constrained Random Forest")
}
