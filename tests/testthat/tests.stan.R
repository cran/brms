test_that("Test that stan.reffects contains the correct strings", {
  expect_match(stan.reffects(list(c("Intercept","PROD"), "site"), f = c("Intercept","Prod"))[5], 
               "cor_site\\[1,2\\]", all = FALSE)
  expect_equal(stan.reffects(list(c("Intercept"), "site"), f = c("Intercept","Prod"))[5], 
               "  sd_site_Intercept <- sd_site; \n", all = FALSE)
})

test_that("Test that stan.model accepts supported links", {
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "sratio", 
                        link="probit_approx"), "Phi_approx")
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "cumulative", 
                        link="probit"), "Phi")
  expect_match(stan.model(rating ~ treat + period + carry, data = inhaler, family = "poisson", 
                        link="log"), "log")
})

test_that("Test that stan.prior accepts supported prior families", {
  expect_equal(stan.prior("b_x1", prior = list(b = "uniform(0,10)")), 
               "  b ~ uniform(0,10); \n")
  expect_equal(stan.prior(c("b_x1","b_x2"), prior = list(b = "uniform(0,10)", 
               b_x1 = "normal(0,1)"), ind = 1:2), 
               c("  b[1] ~ normal(0,1); \n", "  b[2] ~ uniform(0,10); \n"))
})

test_that("Test that stan.prior returns the correct indices", {
  expect_equal(stan.prior("sd_Intercept"), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan.prior("sd_Intercept", ind = "k"), 
               "  sd ~ cauchy(0,5); \n")
  expect_equal(stan.prior("sd_Intercept", ind = "k", prior = list(sd_Intercept = "normal(0,1)")), 
               "  sd[k] ~ normal(0,1); \n")
  expect_equal(stan.prior(c("sd_x1","sd_x2"), ind = 1:2, prior = list(sd_x1 = "normal(0,1)")),
               c("  sd[1] ~ normal(0,1); \n","  sd[2] ~ cauchy(0,5); \n"))                                                       
})

#test_that("Test that brm.stan does not duplicate the model", {
#  expect_match()
#})


