#' Plot the prevalence estimates
#'
#' @param fit A rater fit object.
#'
#' @return A plot of the prevalence estimates extracted from the fit.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text coord_cartesian labs
#'     theme_bw
#' @importFrom rlang .data
#'
#' @noRd
#'
plot_pi <- function(fit) {
  pi <- pi_point_estimate(fit)
  pi_ci <- as.data.frame(posterior_interval(fit,pars="pi"))
  plot_data <- data.frame(cat = as.factor(1:length(pi)),
                          pi = pi,
                          round_pi = round(pi, 2),
                          upper=pi_ci$"95%",
                          lower=pi_ci$"5%")

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$cat, y = .data$pi)) +
    ggplot2::geom_point(stat="identity")+
    ggplot2::geom_errorbar (mapping = ggplot2::aes(ymin=.data$lower,ymax=.data$upper)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$round_pi), hjust = -.5) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Category",
                  y = "Prevalence prob.") +
    ggplot2::theme_bw() +
    NULL

  plot
}

#' Plot the posterior distribution of prevalence probabilities
#'
#' @param fit A rater fit object
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text coord_cartesian labs
#'     theme_bw
#'
#' @importFrom tidybayes stat_eye
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr everything
#' @return A plot of posterior distributions of prevalence probabilities of latent classes
#' @export
#'
#'
plot_pi_dist <- function(fit) {
  pi_data <- as.data.frame(posterior_samples(fit,pars="pi"))
  pi_data <- tidyr::pivot_longer(pi_data,cols=dplyr::everything(),values_to="pi",names_to="cat",names_prefix="pi.", names_ptypes=factor())
  #pi_mean <- dplyr::group_by(pi_data,cat)
  #pi_mean<- dplyr::summarise( pi_mean, cat=unique(cat), label=round(mean(pi),2), pi_mean=mean(pi) )
  plot <- ggplot2::ggplot(pi_data, ggplot2::aes(x = .data$cat, y = .data$pi)) +
    tidybayes::stat_eye(fill = "steelblue") +
    #ggplot2::geom_text(data=pi_mean,mapping=ggplot2::aes(label = .data$label,y=.data$pi_mean,x=.data$cat), hjust = -.9) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Category",
                  y = "Prevalence prob.") +
    ggplot2::theme_bw() +
    NULL

  plot
}


#' Plot the rater accuracy estimates
#'
#' @param fit rater fit object
#' @param which which raters to plot
#'
#' @return Plot of the rate accuracy estimates
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text facet_wrap labs guides
#'      scale_fill_gradient theme_bw theme element_rect element_blank
#' @importFrom rlang .data
#'
#' @noRd
#'
plot_theta <- function(fit, which = NULL) {
  theta <- theta_point_estimate(fit, which = which)

  # theta will always have dim[[2]] and it will always be == K
  K <- dim(theta)[[2]]

  # would be great if we could treat in arrays and matrices the 'same'
  if (length(dim(theta)) > 2) {
    J <- dim(theta)[[1]]
    value <- unlist(lapply(1:J, function(x) as.vector(theta[x, , ])))
  } else {
    J <- 1
    value <- as.vector(theta)
  }
  which <- if (is.null(which)) 1:J else which

  plot_data <- data.frame(
                  x = factor(rep(rep(1:K, each = K), J), levels = 1:K),
                  y = factor(rep(rep(1:K, K), J), levels = K:1),
                  rater = rep(which, each = K^2),
                  value = value,
                  round_value = round(value, 2))
  rownames(plot_data) <- NULL

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y)) +
   ggplot2::geom_tile(ggplot2::aes(fill = .data$value), col = "black") +
   ggplot2::geom_text(ggplot2::aes(label = .data$round_value)) +
   ggplot2::facet_wrap(~ rater) +
   # TODO add way to change defaults
   ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
   ggplot2::labs(y = "True label",
                 x = "Assigned label") +
   ggplot2::guides(fill = FALSE) +
   ggplot2::theme_bw() +
   ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                  panel.grid.major = ggplot2::element_blank(),
                  panel.grid.minor = ggplot2::element_blank(),
                  panel.border     = ggplot2::element_blank()) +
   NULL

  plot
}

#' Plot the latent class estimates of a rater fit.
#'
#' @param x numeric matrix object
#' @param ... Other arguments
#'
#' @return Plot of the rate accuracy estimates
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text labs theme_bw theme
#'     scale_fill_gradient guides element_blank
#' @importFrom rlang .data
#'
#' @noRd
#'
plot_class_probabilities <- function(fit) {

  x <- class_probabilities(fit)

  # We could validate more stringently here if required
  if (!is.numeric(x)) {
    stop("Can only plot numeric matrices.", call. = FALSE)
  }

  I <- nrow(x)
  K <- ncol(x)

  plot_data <- data.frame(x = factor(rep(1:K, each = I), levels = 1:K),
                          y = factor(rep(1:I, K), levels = I:1),
                          prob = as.vector(x),
                          round_prob = round(as.vector(x), 2))

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$prob), colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = .data$round_prob)) +
    ggplot2::labs(x = "Latent Class",
                  y = "Item") +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::guides(fill = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border     = ggplot2::element_blank()) +
    NULL

  plot
}
