# Regression models in Stan
# 
# @inheritParams brm
# @return A character string containing the model in stan language
# @examples
# \dontrun{
# stan.model <- brm.stan(y ~ log_Age + log_Base4 * Trt, data = epilepsy, 
#                   family = "poisson", link = "log", prior = list(b_ = "normal(0,5)"))
# }
stan.model <- function(formula, data = NULL, family = "gaussian", link = "identity",
                       prior = list(), partial = NULL, threshold = "flexible",
                       predict = FALSE, save.model = NULL, ...) {
  ef <- extract.effects(formula = formula, partial = partial)  
  data <- model.frame(ef$all, data = data, drop.unused.levels = TRUE)

  is.lin <- family %in% c("gaussian", "student", "cauchy")
  is.ord <- family %in% c("cumulative", "cratio", "sratio", "acat") 
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  is.count <- family %in% c("poisson", "negbinomial", "geometric")
    
  if (family == "categorical") {
    X <- data.frame()
    Xp <- brm.model.matrix(ef$fixed, data, rm.int = is.ord)
  }
  else {
    X <- brm.model.matrix(ef$fixed, data, rm.int = is.ord)
    Xp <- brm.model.matrix(partial, data, rm.int = TRUE)
  }  
  f <- colnames(X)
  p <- colnames(Xp)
  Z <- lapply(ef$random, brm.model.matrix, data = data)
  r <- lapply(Z,colnames)
  n <- ifelse(is(ef$add, "formula") & (is.ord | family %in% c("binomial", "categorical")), "[n]", "")
  
  type <- ifelse(is.lin | is.skew, "real", "int")
  data.text <- paste0(
    "data { \n", paste0(
    "  int<lower=1> N; \n",
    "  ",type," Y[N]; \n"), 
    if (ncol(X)) paste0("  int<lower=1> K; \n", "  matrix[N,K] X; \n"),
    if (length(p)) paste0("  int<lower=1> Kp; \n", "  matrix[N,Kp] Xp; \n"),  
    if (is.lin & is(ef$add, "formula"))
      "  real<lower=0> sigma[N]; \n"
    else if (is.lin & is(ef$add2, "formula"))
      "  vector<lower=0>[N] inv_weights; \n"
    else if (is.ord | is.element(family, c("binomial", "categorical"))) 
      paste0("  int max_obs",toupper(n),"; \n")
    else if (is.skew & is(ef$add, "formula"))
      "  vector[N] cens; \n",
    paste0(sapply(ef$group, function(g) { 
      r <- r[[match(g,ef$group)]]
      paste0(  "  int<lower=1> ",g,"[N]; \n",
      "  int<lower=1> N_",g,"; \n",
      "  int<lower=1> K_",g,"; \n",
      if (length(r) == 1) paste0(
      "  real Z_",g,"[N]; \n")
      else if (length(r) > 1) paste0(
      "  row_vector[K_",g,"] Z_",g,"[N]; \n"))
    }), collapse = ""), 
    "} \n")
  
  max_obs <- ifelse(rep(is(ef$add, "formula") & (is.ord | family == "categorical"), 3), 
    c("MAX_obs", "  int MAX_obs; \n", "  MAX_obs <- max(max_obs); \n"), c("max_obs", rep("", 2)))
  zero <- ifelse(rep(family == "categorical", 2), c("  row_vector[1] zero; \n", "  zero[1] <- 0; \n"), "")
  trans.data.text <- paste0(
    "transformed data { \n",
      max_obs[2], zero[1], max_obs[3], zero[2],
    "} \n")
  
  ord.intercept <- ifelse(is.ord, paste0("  ", ifelse(family == "cumulative",
                   "ordered", "vector"), "[",max_obs[1],"-1] b_Intercept; \n"), "")
  par.text <- paste0(
    "parameters { \n",
    if (length(f)) "  vector[K] b; \n",
    if (length(p)) paste0("  matrix[Kp,",max_obs[1],"-1] bp; \n"),
    if (threshold == "flexible") ord.intercept
    else if (threshold == "equidistant") paste0(
      "  real b_Intercept1; \n",  
      "  real", ifelse(family == "cumulative","<lower=0>", "")," delta; \n"),
    paste0(sapply(ef$group, function(g) {
      r <- r[[match(g, ef$group)]]
      paste0(ifelse(length(r) == 1, paste0(
      "  real r_",g,"[N_",g,"]; \n",
      "  real<lower=0> sd_",g,"; \n"), paste0(
      "  vector[K_",g,"] r_",g,"[N_",g,"]; \n",
      "  vector<lower=0>[K_",g,"] sd_",g,"; \n",
      "  cholesky_factor_corr[K_",g,"] L_",g,"; \n")))}), 
      collapse = ""),
    if (is.lin & !is(ef$add, "formula")) 
      "  real<lower=0> sigma; \n",
    if (family == "student") 
      "  real<lower=1> nu; \n",
    if (family %in% c("gamma", "weibull", "negbinomial")) 
      "  real<lower=1e-10> shape; \n",
    "} \n")
  
  ilink <- c(identity = "", log = "exp", inverse = "inv", sqrt = "square", logit = "inv_logit", 
            probit = "Phi", probit_approx = "Phi_approx", cloglog = "inv_cloglog")[[link]]
  simplify <- family %in% c("binomial","bernoulli","cumulative") & link == "logit" & 
              !(predict & family == "cumulative")
  vectorize <- c(!length(ef$random), 
                 is.element(family,c("binomial","bernoulli")) & link=="logit" | 
                 is.lin & link != "inverse" | is.count & link !="sqrt")
  if (family == "categorical") fe.only <- f
  else fe.only <- setdiff(f, unlist(lapply(r, intersect, y = f)))
  eta.fe <- paste0("  eta <- ", ifelse(length(f), "X * b", "rep_vector(0,N)"), "; \n",
                   ifelse(length(p), "  etap <- Xp * bp; \n", ""))
  eta.re <- sapply(ef$group, function(g) paste0("Z_",g,"[n]*r_",g,"[",g,"[n]]"))
  if (length(eta.re)) 
    eta.re <- paste0("    eta[n] <- eta[n] + ", paste0(eta.re, collapse = " + "), "; \n")
  etap.init <- ifelse(length(p), paste0("  matrix[N,",max_obs[1],"-1] etap; \n"), "")

  reffects <- unlist(lapply(mapply(list, r, ef$group, SIMPLIFY = FALSE), stan.reffects, f = f, 
                                 family = family, prior = prior))
  if (length(reffects)) reffects <- sapply(1:5, function(x) 
    paste0(reffects[seq(x, length(reffects), 5)], collapse = ""))   
  ord <- stan.ord(family, ilink = ilink, partial = length(p), max_obs = max_obs[1], n = n, 
                  predict = predict) 
  
  llh <- stan.llh(family, link = link, add = is(ef$add, "formula"), 
                 add2 = is(ef$add2, "formula"))
  llh.pred <- stan.llh(family, link = link, predict = TRUE, add = is(ef$add, "formula"),  
                 add2 = is(ef$add2, "formula")) 
  if (is.skew & is(ef$add, "formula"))
    cens <- c("if (cens[n] == 0) ", paste0("    else ",
      stan.llh(family, link = link, add = is(ef$add, "formula"), cens = TRUE), "\n"))
  else cens <- rep("", 2)
                      
  priors <- paste0(
    if (length(f)) paste0(stan.prior(paste0("b_",f), prior = prior, ind = 1:length(f)), collapse = ""),
    if (is.ord & threshold == "flexible") 
      stan.prior("b_Intercept", prior = prior, add.type = "Intercept")
    else if (is.ord & threshold == "equidistant") 
      paste0(stan.prior("b_Intercept1", prior = prior, add.type = "Intercept1"),
             stan.prior("delta", prior = prior)),
    if (length(p)) paste0(stan.prior(paste0("b_",p), prior = prior, 
                    ind = 1:length(p), partial = TRUE), collapse = ""), 
    if (is.element(family,c("gamma", "weibull"))) stan.prior("shape", prior = prior),
    if (family == "student") stan.prior("nu", prior = prior),
    if (is.lin & !is(ef$add, "formula")) stan.prior("sigma", prior = prior), reffects[1])
  loop <- !vectorize[1] | is.ord & !(is.ord & simplify & !is(ef$add, "formula"))
  
  model <- paste0(data.text, 
  trans.data.text, par.text,
  "transformed parameters { \n",
    "  vector[N] eta; \n", etap.init, ord[1], 
    if (threshold == "equidistant") ord.intercept,
    reffects[2],
    eta.fe, if (threshold == "equidistant") 
      paste0("  for (k in 1:(",max_obs[1],"-1)) { \n",
      "    b_Intercept[k] <- b_Intercept1 + (k-1.0)*delta; \n  } \n"),
    if(loop) "  for (n in 1:N) { \n",
    eta.re, ord[2], if (loop)"  } \n",
    reffects[3], 
  "} \n",
  "model { \n",
    priors, 
    ifelse(vectorize[2],paste0("  Y ~ ", llh),
      paste0("  for(n in 1:N) { \n    ",
        cens[1],"Y[n] ~ ",llh, cens[2],"  } \n")), 
  "} \n",
  "generated quantities { \n",
    if (length(f)) paste0(
    "  real b_",f,"; \n", collapse =""),
    if (length(p)) paste0(
    "  vector[",max_obs[1],"-1] b_",p,"; \n", collapse =""),
    reffects[4],
    if (predict) paste0(
    "  ",type," Y_pred[N]; \n", 
    "  for (n in 1:N) { \n", 
    "    Y_pred[n] <- ",llh.pred,
    "  } \n"),
    if (length(f)) paste0(
    "  b_",f," <- b[",1:length(f),"]; \n", collapse = ""),
    if (length(p)) paste0(
    "  b_",p," <- to_vector(bp[",1:length(p),"]); \n", collapse = ""),
    reffects[5], 
  "} \n")
  
  if (is.character(save.model)) {
    sink(save.model)
    cat(model)
    sink()
  }
  model
}

# Random effects in Stan 
# 
# @return A vector of strings containing the random effects in stan language
stan.reffects <- function(rg, f, family = "gaussian", prior = list()) {
  r <- rg[[1]]
  g <- rg[[2]]
  is.ord <- is.element(family, c("cumulative", "cratio", "sratio", "acat")) 
  out <- rep("", 5)
  out[1] <- paste0(stan.prior(paste0("sd_",g,"_",r), add.type = g, prior = prior ,
                     ind = ifelse(length(r) == 1, "", list(1:length(r)))[[1]]))
  if (length(r) == 1) {
    m <- ifelse(r == "Intercept" & is.ord | family == "categorical" | !is.element(r, f), "0", 
                paste0("b[", which(r==f), "]"))
    out[1] <- paste0(out[1],"  r_",g," ~ normal(",m,",sd_",g,"); \n")
    out[4] <- paste0("  real<lower=0> sd_",g,"_",r,"; \n")
    out[5] <- paste0("  sd_",g,"_",r," <- sd_",g,"; \n")
  }  
  else if(length(r) > 1) {
    out[1] <- paste0(out[1], stan.prior(paste0("L_",g), prior = prior, add.type = g),
                     "  r_",g," ~ multi_normal_cholesky(mu_",g,",diag_pre_multiply(sd_",g,",L_",g,")); \n")
    out[2] <- paste0("  vector[K_",g,"] mu_",g,"; \n")
    out[3] <- paste0(sapply(1:length(r), function(i) {
      if (is.element(r[i],f) & family != "categorical") paste0("  mu_",g,"[",i,"] <- b[",which(r[i]==f),"]; \n") 
      else paste0("  mu_",g,"[",i,"] <- 0; \n")}), collapse = "")
    out[4] <- paste0(paste0("  real<lower=0> sd_",g,"_",r,"; \n", collapse = ""),
                     "  corr_matrix[K_",g,"] cor_",g,"; \n",
                     paste0(unlist(lapply(2:length(r),function(i) lapply(1:(i-1), function(j)
                       paste0("  real<lower=-1,upper=1> cor_",g,"_",r[j],"_",r[i],"; \n")))),
                       collapse = ""))
    out[5] <- paste0(paste0("  sd_",g,"_",r," <- sd_",g,"[",1:length(r),"]; \n", collapse = ""),
                     "  cor_",g," <- multiply_lower_tri_self_transpose(L_",g,"); \n",
                     paste0(unlist(lapply(2:length(r),function(i) lapply(1:(i-1), function(j)
                       paste0("  cor_",g,"_",r[j],"_",r[i]," <- cor_",g,"[",j,",",i,"]; \n")))),
                       collapse = ""))       
  }
  out
}

# Ordinal effects in Stan
# 
# @return A vector of strings containing the ordinal effects in stan language
stan.ord <- function(family, ilink = "inv_logit", partial = FALSE, max_obs = "max_obs", 
                         n = "", predict = FALSE) {
  simplify <- is.element(family, c("cumulative", "categorical")) & ilink == "inv_logit" & n != "[n]" & !predict
  if (!is.element(family, c("cumulative", "cratio", "sratio", "acat", "categorical")) | simplify) 
    return(rep("", 3)) 
  th <- function(k) {
    sign <- ifelse(is.element(family, c("cumulative", "sratio"))," - ", " + ")
    ptl <- ifelse(partial, paste0(sign, "etap[n,k]"), "") 
    if (sign == " - ") out <- paste0("b_Intercept[",k,"]", ptl, " - eta[n]")
    else out <- paste0("eta[n]", ptl, " - b_Intercept[",k,"]")
  }  
  add.loop <- ifelse(n == "[n]", paste0("    for (k in (max_obs[n]+1):MAX_obs) p[n,k] <- 0.0; \n"), "")
  hd <- ifelse(rep(n == "[n]" & is.element(family, c("sratio","cratio")), 2), 
               c("head(", paste0(",max_obs",n,"-1)")), "")
  sc <- ifelse(family=="sratio","1-","")
  ord <- paste0("  vector[",max_obs,"] p[N]; \n", 
                if (!is.element(family, c("cumulative", "categorical"))) 
                  paste0("  vector[",max_obs,"-1] q[N]; \n"))
  if (family == "categorical" & ilink == "inv_logit") ord[2] <- paste0(
    "    p[n,1] <- 1.0; \n",
    "    for (k in 2:max_obs",n,") { \n",
    "      p[n,k] <- exp(eta[n,k-1]); \n",
    "    } \n", add.loop,
    "    p[n] <- p[n]/sum(p[n]); \n")
  else if (family == "cumulative") ord[2] <- paste0(
    "    p[n,1] <- ",ilink,"(",th(1),"); \n",
    "    for (k in 2:(max_obs",n,"-1)) { \n", 
    "      p[n,k] <- ",ilink,"(",th("k"),") - ",ilink,"(",th("k-1"),"); \n", 
    "    } \n",
    "    p[n,max_obs",n,"] <- 1 - ",ilink,"(",th(paste0("max_obs",n,"-1")),"); \n", 
    add.loop)
  else if (is.element(family,c("sratio", "cratio"))) ord[2] <- paste0(
    "    for (k in 1:(max_obs",n,"-1)) { \n",
    "      q[n,k] <- ",sc, ilink,"(",th("k"),"); \n",
    "      p[n,k] <- 1-q[n,k]; \n",
    "      for (kk in 1:(k-1)) p[n,k] <- p[n,k] * q[n,kk]; \n", 
    "    } \n",
    "    p[n,max_obs",n,"] <- prod(",hd[1],"q[n]",hd[2],"); \n",
    add.loop)
  else if (family == "acat") 
    if (ilink == "inv_logit") ord[2] <- paste0(
      "    p[n,1] <- 1.0; \n",
      "    for (k in 1:(max_obs",n,"-1)) { \n",
      "      q[n,k] <- ",th("k"),"; \n",
      "      p[n,k+1] <- q[n,1]; \n",
      "      for (kk in 2:k) p[n,k+1] <- p[n,k+1] + q[n,kk]; \n",
      "      p[n,k+1] <- exp(p[n,k+1]); \n",
      "    } \n", add.loop,
      "    p[n] <- p[n]/sum(p[n]); \n")
  else ord[2] <- paste0(                   
    "    for (k in 1:(max_obs",n,"-1)) \n",
    "      q[n,k] <- ",ilink,"(",th("k"),"); \n",
    "    for (k in 1:max_obs",n,") { \n",     
    "      p[n,k] <- 1.0; \n",
    "      for (kk in 1:(k-1)) p[n,k] <- p[n,k] * q[n,kk]; \n",
    "      for (kk in k:(max_obs",n,"-1)) p[n,k] <- p[n,k] * (1-q[n,kk]); \n",      
    "    } \n", add.loop,
    "    p[n] <- p[n]/sum(p[n]); \n")
  ord
}

# Priors in Stan
# 
# Define priors for parameters in Stan language
# 
# @param par A vector of parameter names
# @param prior A named list of strings containing the priors for parameters in \code{par}.
# @param ind An optional index to allow for different priors for different elements of a parameter vector (see 'Examples')
# @param s An integer >= 0 defining the number of spaces in front of the output string
# @param ... Other potential arguments
# @inheritParams brm
# 
# @return A vector of character strings each defining a line of code in stan language that
#   defines the prior of a parameter in \code{par}. If a parameter has has no corresponding prior in \code{prior} 
#   and also no internal default in \code{stan.prior} (see 'Details'), an empty string is returned.
#      
# @examples 
# stan.prior(c("b_x1","sd_Site","sd_obs"), 
#           prior = list(b_ = "normal(0,1)", sd = "cauchy(0,2.5)"))
#  
# # returns a cauchy prior for sd_Site and a uniform prior for sd_obs                
# stan.prior(c("sd_Site","sd_obs"), 
#           prior = list(sd = "cauchy(0,5)", sd_obs = "uniform(0,100)"))            
# 
# # returns a uniform prior for the first element of b                      
# stan.prior("b", prior = list(beta = "uniform(0,10)"), ind = "[1]")
# 
# # returns the default prior for nu
# stan.prior("nu")
stan.prior = function(par, prior = list(), add.type = NULL, ind = rep("", length(par)), 
                      partial = FALSE, s = 2) { 
  if (length(par) != length(ind)) 
    stop("Arguments par and ind must have the same length")
  if (length(par) > 1 & all(ind == ""))
    warning("Indices are all zero")
  type <- unique(unlist(lapply(par, function(par) 
    unlist(regmatches(par, gregexpr("[^_]*", par)))[1])))
  if (length(type) > 1) stop("Only one parameter type can be handled at once")
  if (!is.null(add.type)) {
    type <- c(type, paste0(type,"_",add.type))
    if (!all(grepl(type[2], par))) 
      stop("Additional parameter type not present in all parameters")
  }  
  default.prior <- list(b = "", bp = "", sigma = "", delta = "",
    L = "lkj_corr_cholesky(1)", sd = "cauchy(0,5)", nu = "uniform(1,60)", 
    shape = "gamma(0.01,0.01)") 
  if (!is.null(prior[[type[2]]])) base.prior <- prior[[type[2]]]
  else if (!is.null(prior[[type[1]]])) base.prior <- prior[[type[1]]]
  else base.prior <- default.prior[[type[1]]]
  
  if (type[1] == "b" & partial) type[1] <- "bp"
  type <- ifelse(is.na(type[2]), type[1], type[2])  
  if (any(par %in% names(prior)))
    out <- sapply(1:length(par), function(i, par, ind) {
      if (ind[i] != "") ind[i] <- paste0("[",ind[i],"]")
      if (!par[i] %in% names(prior))
        prior[[par[i]]] <- base.prior
      if (paste0(prior[[par[i]]],"") != "")  
        return(paste0(paste0(rep(" ", s), collapse = ""), 
               type, ind[i], " ~ ", prior[[par[i]]], "; \n"))
      else return("") }, 
    par = par, ind = ind)
  else if (base.prior != "")
    out <- paste0(paste0(rep(" ", s), collapse = ""), 
      ifelse(identical(type,"bp"), "to_vector(bp)", type), " ~ ", base.prior, "; \n")
  else out <- ""
  out
}

# Likelihoods in stan language
# 
# Define the likelihood of the dependent variable in stan language
# 
# @inheritParams brm
# @param add A flag inicating if the model contains additional information of the response variable
#   (e.g., standard errors in a gaussian linear model)
# @param add2 A flag indicating if the response variable should have unequal weights 
#   Only used if \code{family} is either \code{"gaussian", "student"}, or \code{"cauchy"}.
#    
# @return A character string defining a line of code in stan language 
#   that contains the likelihood of the dependent variable. 
# @examples 
# \dontrun{
# stan.llh(family = "gaussian")
# stan.llh(family = "cumulative", link = "logit")
# }
stan.llh <- function(family, link = "identity", predict = FALSE, add = FALSE,
                    add2 = FALSE, cens = FALSE, engine = "stan") {
  is.ord <- is.element(family, c("cumulative", "cratio", "sratio", "acat"))
  n <- ifelse(predict | is.element(link, c("inverse","sqrt")), "[n]", "")
  ilink <- c(identity = "", log = "exp", inverse = "inv", sqrt = "square", logit = "inv_logit", 
             probit = "Phi", probit_approx = "Phi_approx", cloglog = "inv_cloglog")[[link]]
  simplify <- !predict & (link=="logit" & (is.element(family,c("binomial", "bernoulli"))
                            | is.element(family, c("cumulative", "categorical")) & !add) 
                            | family %in% c("poisson","negbinomial", "geometric") & link == "log") 
  lin.args <- paste0(ilink,"(eta",n,"),sigma", ifelse(add, n, ifelse(add2,"*inv_weights","")))
  if (simplify) 
    llh <- list(poisson = c("poisson_log", "(eta);"), 
            negbinomial = c("neg_binomial_2_log", "(eta,shape);"),
            geometric = c("neg_binomial_2_log", "(eta,1);"),
            cumulative = c("ordered_logistic", "(eta[n],b_Intercept);"),
            categorical = c("categorical_logit", "(to_vector(append_col(zero, eta[n] + etap[n])));"), 
            binomial=c("binomial_logit", "(max_obs,eta);"), 
            bernoulli=c("bernoulli_logit", "(eta);"))[[family]]
  else llh <- list(gaussian = c("normal", paste0("(",lin.args,");")),
               student = c("student_t", paste0("(nu,",lin.args,");")),
               cauchy = c("cauchy", paste0("(",lin.args,");")),
               poisson = c("poisson", paste0("(",ilink,"(eta",n,"));")),
               negbinomial = c("neg_binomial_2", paste0("(",ilink,"(eta",n,"),shape);")),
               geometric = c("neg_binomial_2", paste0("(",ilink,"(eta",n,"),1);")),
               binomial = c("binomial",paste0("(max_obs",ifelse(add,"[n]",""),",",ilink,"(eta[n]));")),
               bernoulli = c("bernoulli",paste0("(",ilink,"(eta[n]));")), 
               gamma = c("gamma", paste0("(shape, shape*inv(",ilink,"(eta[n])));")), 
               exponential = c("exponential", paste0("(",ilink,"(-eta[n]));")),
               weibull = c("weibull", paste0("(shape, inv(",ilink,"(-eta[n]/shape)));")), 
               categorical = c("categorical","(p[n]);"))[[ifelse(is.ord, "categorical", family)]]
  llh <- paste0(llh[1],ifelse(predict, "_rng",""),llh[2],"\n")
  if (cens & family %in% c("gamma", "exponential", "weibull")) {
    llh <- list(gamma = paste0("gamma_cdf_log(Y[n], shape, shape*inv(",ilink,"(eta[n])))"),
                exponential = paste0("exponential_cdf_log(Y[n],",ilink,"(-eta[n]))"),
                weibull = paste0("weibull_cdf_log(Y[n], shape, inv(",ilink,"(-eta[n]/shape)))"))[family]
    llh <- paste0("increment_log_prob(",llh,");")
  }  
  llh
}
