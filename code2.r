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

# Task 4
tpowers <- c()
rpowers <- c()
delta.list <- seq(0, 2, 0.01)
for (d in delta.list) {
  pvalues <- replicate(1000, t.test(rnorm(10, -d/2), rnorm(10, d/2))$p.value)
  tpowers <- append(tpowers, power.t.test(10, delta = d)$power)
  rpowers <- append(rpowers, sum(pvalues < 0.05) / 1000)
}
plot(delta.list, rpowers, type = "l", col = "blue")
lines(delta.list, tpowers, col = "red")

# Task 5
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
  plot(delta.list, v.equal.power, type = "l", col = "blue", main = paste0("Plot with sample size ", sample.size))
  points(delta.list, v.unequal.power, col = "red", pch = 20, cex = 0.5)
}