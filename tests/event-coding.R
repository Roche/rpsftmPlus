library(rpsftmPlus)
library(survival)

# The same patients, represented with the two supported numeric status codes.
trial01 <- data.frame(
  t = c(2, 3, 4, 6, 8, 10, 3, 5, 6, 8, 10, 12),
  event = c(1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0),
  arm = rep(0:1, each = 6),
  rx = c(.5, 0, .5, 0, .5, 0, rep(1, 6)))
trial12 <- transform(trial01, event = event + 1)
fit01 <- rpsftm(Surv(t, event) ~ rand(arm, rx), data = trial01,
               low_psi = -3, hi_psi = 3)
fit12 <- rpsftm(Surv(t, event) ~ rand(arm, rx), data = trial12,
               low_psi = -3, hi_psi = 3)
stopifnot(identical(fit01$Sstar, fit12$Sstar))

counts <- rpsftm.atrisk(fit12, eval.times = 20)
stopifnot(all(counts$Count[counts$Status == "Event"] == 4),
          all(counts$Count[counts$Status == "Censored"] == 2),
          isTRUE(all.equal(counts, rpsftm.atrisk(fit01, eval.times = 20))))
cox01 <- rpsftm.coxph(fit01)
cox12 <- rpsftm.coxph(fit12)
stopifnot(cox12$n == 12, cox12$nevent == 8,
          isTRUE(all.equal(coef(cox01), coef(cox12))),
          isTRUE(all.equal(vcov(cox01), vcov(cox12))))
