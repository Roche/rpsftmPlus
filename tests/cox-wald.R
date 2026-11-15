library(rpsftmPlus)
library(survival)

trial <- data.frame(t = c(1:8, 1.4 * (1:8)), event = 1,
                    arm = rep(0:1, each = 8), rx = rep(c(.75, 1), each = 8))
fit <- rpsftm(Surv(t, event) ~ rand(arm, rx), data = trial,
              low_psi = -5, hi_psi = 3)
cox <- rpsftm.coxph(fit)
s <- summary(cox)
# For the single treatment coefficient, Wald = beta^2 / adjusted variance.
expected <- as.numeric(coef(cox)^2 / vcov(cox))
expected.p <- pchisq(expected, df = 1, lower.tail = FALSE)
stopifnot(isTRUE(all.equal(as.numeric(cox$wald.test), expected)),
          isTRUE(all.equal(unname(s$waldtest["pvalue"]), expected.p)),
          isTRUE(all.equal(unname(s$coefficients[1, "Pr(>|z|)"]), expected.p)),
          expected.p > .05,
          s$conf.int[1, "lower .95"] < 1, s$conf.int[1, "upper .95"] > 1)
