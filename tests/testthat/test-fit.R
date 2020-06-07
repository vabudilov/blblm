test_that("fitting test", {
  set.seed(3)
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  res<-coef(fit)
  itwas<-c(51.91244091, -8.96848121, -0.13196459,  0.03239431)
 # expect_equal(sum(res), sum(itwas))
#  expect_equal(res[1], itwas[1])
#  expect_equal(sum(res[1]), 51.91244091)
  expect_equal(sum(itwas-res), 0)

  set.seed(345)
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, ifgeneral = TRUE)
  res<-coef(fit)
  itwas<-c( 56.63892608, -10.17525530,  -0.17983333,   0.04461437 )
  expect_equal(sum(itwas-res), 0)

    })
