context("Data generators")

test_that("Generate mixture data (binary)", {

	dat <- gen_mixture()
	plot(dat)
	expect_s3_class(dat, "iprobitData")

})

test_that("Generate mixture data (multinomial)", {

  dat <- gen_mixture(m = 5, proportion = c(0.4, 0.3, 0.2, 0.05, 0.05))
  plot(dat)
  expect_s3_class(dat, "iprobitData")

})

test_that("Generate circle data (binary)", {

  dat <- gen_circle()
  plot(dat)
  expect_s3_class(dat, "iprobitData")

})

test_that("Generate circle data (multinomial)", {

  dat <- gen_circle(m = 5)
  plot(dat)
  expect_s3_class(dat, "iprobitData")

})

test_that("Generate spiral data (binary)", {

  dat <- gen_spiral()
  plot(dat)
  expect_s3_class(dat, "iprobitData")

})

test_that("Generate spiral data (multinomial)", {

  dat <- gen_spiral(m = 5, cycles = 5, sd = 0.1)
  plot(dat)
  expect_s3_class(dat, "iprobitData")

})