
test_that("compute_NPV_works", {

testdata = c(2, 4)
short = compute_NPV(100, time=testdata[1], 0.01)
long = compute_NPV(100, time = testdata[2], 0.01)
expect_true(short > long)
expect_equal(compute_NPV(value=0,time=100,discount=0.1), 0)
expect_equal(compute_NPV(value=100, time=0, discount=0.1), 100)
})

test_that("compute_NPV_haserrorchecking", {
  # throw an error if time is negative 
  expect_error(compute_NPV(value=100, time=-1, discount=-0.1))  
  # throw an error if discount is negative
  expect_error(compute_NPV(value=100, time=1, discount=-0.1))  
  # throw a warning if value is negative
  expect_warning(compute_NPV(value=-100, time=1, discount=0.1))
})