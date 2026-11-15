library(rpsftmPlus)
library(survival)

# With no censoring or switching, an exponential model equates arm means.
trial <- data.frame(t = c(1, 2, 3, 20, 4, 5, 6, 7), event = 1,
                    arm = rep(0:1, each = 4), rx = rep(0:1, each = 4))
fit <- rpsftm(Surv(t, event) ~ rand(arm, rx), data = trial,
              test = survreg, dist = "exponential", alpha = .1,
              low_psi = -4, hi_psi = 4)
refit <- rpsftm.refit(fit, k = 1)
expected <- log(mean(trial$t[trial$arm == 0]) / mean(trial$t[trial$arm == 1]))
stopifnot(inherits(refit, "survreg"),
          abs(refit$psi - expected) < 1e-5,
          isTRUE(all.equal(refit$CI, fit$CI)), refit$fit$conf.int == .9)

# For switchers, change only the requested modifier; keep test and ... .
trial <- data.frame(t = c(1:8, 2 * (1:8)), event = 1,
                    arm = rep(0:1, each = 8), rx = rep(c(.25, 1), each = 8))
fit <- rpsftm(Surv(t, event) ~ rand(arm, rx), data = trial,
              test = survreg, dist = "exponential", alpha = .1,
              low_psi = -4, hi_psi = 4)
for (k in c(1, .8)) {
  trial$modifier <- ifelse(trial$arm == 1, 1, k)
  expected <- rpsftm(Surv(t, event) ~ rand(arm, rx), data = trial,
                     test = survreg, dist = "exponential", alpha = .1,
                     low_psi = -4, hi_psi = 4, treat_modifier = modifier)
  actual <- rpsftm.refit(fit, k)
  stopifnot(inherits(actual, "survreg"),
            isTRUE(all.equal(actual$psi, expected$psi)),
            isTRUE(all.equal(actual$CI, expected$CI)),
            isTRUE(all.equal(actual$eval_z, expected$eval_z)),
            isTRUE(all.equal(actual$Sstar, expected$Sstar)))
}

# A separately supplied formula may not retain the caller's option bindings.
# Keep that unsupported case explicit rather than use same-named outer values.
external.formula <- Surv(t, event) ~ rand(arm, rx)
chosen.alpha <- .05
make.fit <- function() {
  chosen.alpha <- .1
  rpsftm(external.formula, data = trial, alpha = chosen.alpha,
          low_psi = -4, hi_psi = 4)
}
fit <- make.fit()
error <- try(rpsftm.refit(fit, 1), silent = TRUE)
stopifnot(inherits(error, "try-error"),
          grepl("formula supplied directly", as.character(error)))
