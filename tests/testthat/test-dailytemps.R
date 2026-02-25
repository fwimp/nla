test_that("dt_modifiedlp works", {
  lpout <- dt_modifiedlp(10, tbar = seq(20,30,1), clamp = TRUE)
  expect_shape(lpout, dim = c(24,11))
  expect_true(is.numeric(lpout))
})

test_that("dt_sin works", {
	sinout <- dt_sin(10, tbar = seq(20,30,1), clamp = TRUE)
	expect_shape(sinout, dim = c(24,11))
	expect_true(is.numeric(sinout))
})

test_that("dt_tri works", {
	triout <- dt_tri(10, tbar = seq(20,30,1), clamp = TRUE)
	expect_shape(triout, dim = c(24,11))
	expect_true(is.numeric(triout))
})

test_that("dt_modifiedlp clamps appropriately", {
	lpout_unclamped <- dt_modifiedlp(10, tbar = seq(20,30,1), clamp = FALSE)
	lpout <- dt_modifiedlp(10, tbar = seq(20,30,1), clamp = TRUE)
	expect_lte(max(lpout), max(lpout_unclamped))
	expect_gte(min(lpout), min(lpout_unclamped))
})

test_that("dt_sin clamps appropriately", {
	sinout_unclamped <- dt_sin(10, tbar = seq(20,30,1), clamp = FALSE)
	sinout <- dt_sin(10, tbar = seq(20,30,1), clamp = TRUE)
	expect_lte(max(sinout), max(sinout_unclamped))
	expect_gte(min(sinout), min(sinout_unclamped))
})

test_that("dt_tri clamps appropriately", {
	triout_unclamped <- dt_tri(10, tbar = seq(20,30,1), clamp = FALSE)
	triout <- dt_tri(10, tbar = seq(20,30,1), clamp = TRUE)
	expect_lte(max(triout), max(triout_unclamped))
	expect_gte(min(triout), min(triout_unclamped))
})
