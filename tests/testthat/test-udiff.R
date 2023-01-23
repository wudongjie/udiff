test_that("multiplication works", {
  usdata_ex <- read.csv(list.files(system.file('extdata', package = 'udiff'), full.names = T)[1])
  #usdata_ex$y <- as.factor(usdata_ex$y)
  #usdata_ex$x <- as.factor(usdata_ex$x)
  #usdata_ex$g <- as.factor(usdata_ex$g)
  browser()
  formula1 <- y ~ x + g
  w <- runif(39294, 0.01, 0.99)
  result <- udiff(formula1, usdata_ex, weights=w)
  browser()
  print(summary(result))
  #fmm2 <- em::em(result, latent=2, algo="cem")
  #print(summary(fmm2))
})
