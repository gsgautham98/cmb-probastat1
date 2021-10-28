# Task 1, 8, 9
s <- rnorm(10)
t <- rnorm(10)
print(t.test(s, t))

p <- rpois(10, 1)
q <- rpois(10, 1)
print(t.test(p, q))

p.large <- rpois(1000, 1)
q.large <- rpois(1000, 1)
print(t.test(p.large, q.large))

# Task 2, 8, 9
results.normal <- replicate(10000, t.test(rnorm(10), rnorm(10))$p.value)
results.poisson <- replicate(10000, t.test(rpois(10, 1), rpois(10, 1))$p.value)
results.poisson.large <- replicate(10000, t.test(rpois(1000, 1), rpois(1000, 1))$p.value)
hist(results.normal, xlab = "p-values", ylab = "Frequency", main = "Histogram of p-values (Normal distribution)")
hist(results.poisson, xlab = "p-values", ylab = "Frequency", main = "Histogram of p-values (Poisson distribution)")
hist(results.poisson.large, xlab = "p-values", ylab = "Frequency", main = "Histogram of p-values (Poisson distribution - large sample)")

# Task 3, 8, 9
test.custom.means <- function(n, x, v, f) {
  pvalues <- c()
  for (i in 1:1000) {
    if (f == 0) {
      s <- rnorm(n, -x, v)
      t <- rnorm(n, x, v)
    }
    else {
      s <- rpois(n, 1) - x
      t <- rpois(n, 1) + x
    }
    pvalues <- append(pvalues, t.test(s, t)$p.value)
  }
  return(pvalues)
}
hist(test.custom.means(10, 0.25, 1, 0), xlab = "p-values", ylab = "Frequency", main = "Histogram of p-values (Normal dist - diff means)")
hist(test.custom.means(10, 0.25, 1, 1), xlab = "p-values", ylab = "Frequency", main = "Histogram of p-values (Poisson dist - diff means)")
hist(test.custom.means(1000, 0.05, 1, 1), xlab = "p-values", ylab = "Frequency", main = "Histogram of p-values (Poisson dist - diff means & larger sample)")

# Task 4, 8
rpowers.normal <- c()
rpowers.poisson <- c()
tpowers <- c()
delta.list <- seq(0, 2, 0.01)
for (d in delta.list) {
  pvalues.normal <- replicate(1000, t.test(rnorm(10, -d/2), rnorm(10, d/2))$p.value)
  pvalues.poisson <- replicate(1000, t.test(rpois(10,1) - d/2, rpois(10, 1) + d/2)$p.value)
  tpowers <- append(tpowers, power.t.test(10, delta = d)$power)
  rpowers.normal <- append(rpowers.normal, sum(pvalues.normal < 0.05) / 1000)
  rpowers.poisson <- append(rpowers.poisson, sum(pvalues.poisson < 0.05) / 1000)
}
plot(delta.list, tpowers, type = "l", col = "blue", xlab = "delta", ylab = "Power", main = "Power comparison between t-tests on different distributions")
lines(delta.list, rpowers.poisson, col = "green")
lines(delta.list, rpowers.normal, col = "red")
legend(x = "topleft", legend = c("Theoretical", "Poisson", "Normal"), cex = 0.8, lty = c(1, 1, 1), col = c("blue", "green", "red"), bty = "n")

# Task 5, 9
for (sample.size in c(10, 100, 1000)) {
  v.equal.power <- c()
  v.unequal.power <- c()
  for (d in delta.list) {
    pvalues.equal <- c()
    pvalues.unequal <- c()
    for (i in 1:1000) {
      s <- rnorm(sample.size, -d/2)
      t <- rnorm(sample.size, d/2)
      pvalues.equal <- append(pvalues.equal, t.test(s, t, var.equal = TRUE)$p.value)
      pvalues.unequal <- append(pvalues.unequal, t.test(s, t)$p.value)
    }
    v.equal.power <- append(v.equal.power, sum(pvalues.equal < 0.05) / 1000)
    v.unequal.power <- append(v.unequal.power, sum(pvalues.unequal < 0.05) / 1000)
  }
  plot(delta.list, v.equal.power, type = "l", col = "blue", xlab = "delta", ylab = "Power", main = paste0("Plot with sample size ", sample.size))
  points(delta.list, v.unequal.power, col = "red", pch = 20, cex = 0.5)
  legend(x = "bottomright", legend = c("Equal variances t-test", "Unequal variances t-test"), cex = 0.8, lty = c(1, 2), col = c("blue", "red"), bty = "n")
}

# Task 6
v.seq <- seq(0.05, 4, 0.05)
f.pos.eq.list <- c()
f.pos.uneq.list <- c()
for (v in v.seq) {
  f.pos.equal <- replicate(1000, t.test(rnorm(100, 0, 0.5), rnorm(100, 0, v), var.equal = TRUE)$p.value) < 0.05
  f.pos.unequal <- replicate(1000, t.test(rnorm(100, 0, 0.5), rnorm(100, 0, v))$p.value) < 0.05
  f.pos.eq.list <- append(f.pos.eq.list, mean(f.pos.equal))
  f.pos.uneq.list <- append(f.pos.uneq.list, mean(f.pos.unequal))
}
f.pos.eq.rate <- f.pos.eq.list / seq_along(v.seq)
f.pos.uneq.rate <- f.pos.uneq.list / seq_along(v.seq)
plot(v.seq, f.pos.eq.rate, type = "l", col = "blue", xlab = "Variance", ylab = "False positive rate", main = "Evolution of false positive rate for a normally-distributed sample")
lines(v.seq, f.pos.uneq.rate, col = "red")
legend(x = "topright", c("Equal variances t-test", "Unequal variances t-test"), cex = 0.8, lty = c(1, 1), col = c("blue", "red"), bty = "n")

# Task 7, 8, 9
n.list <- seq(100, 5100, 250)
t.powers.normal <- c()
w.powers.normal <- c()
t.powers.poisson <- c()
w.powers.poisson <- c()
for (n in n.list) {
  t.powers.normal <- append(t.powers.normal, mean(replicate(1000, t.test(rnorm(n), rnorm(n), var.equal = TRUE)$p.value) < 0.05))
  w.powers.normal <- append(w.powers.normal, mean(replicate(1000, wilcox.test(rnorm(n), rnorm(n))$p.value) < 0.05))
  t.powers.poisson <- append(t.powers.poisson, mean(replicate(1000, t.test(rpois(n, 1), rpois(n, 1))$p.value) < 0.05))
  w.powers.poisson <- append(w.powers.poisson, mean(replicate(1000, wilcox.test(rpois(n, 1), rpois(n, 1))$p.value) < 0.05))
}
plot(n.list, t.powers.normal, main = "Comparison of t-test and Mann-Whitney U test - normal dist", xlab = "Sample size", ylab = "Power", col = "blue", type = "l")
lines(n.list, w.powers.normal)
plot(n.list, t.powers.poisson, main = "Comparison of t-test and Mann-Whitney U test - Poisson dist", xlab = "Sample size", ylab = "Power", col = "blue", type = "l")
lines(n.list, w.powers.poisson)