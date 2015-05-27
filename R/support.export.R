#' Parameters of interest for \code{brms} models
#' 
#' @inheritParams brm
#' @param ranef logical; indicating if random effects estimates should be returned
#' 
#' @return A vector of character strings specifying parameters of interest for models produced by the \code{brms} package.
#' @examples 
#' brm.pars(rating ~ treat + period + carry + (1|subject),
#'          data = inhaler, family = "cumulative")
#
#  brm.pars(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit),
#         data = epilepsy, family = c("poisson", "log"))
#'          
#' @export
brm.pars = function(formula, data = NULL, family = "gaussian", partial = NULL,
                    threshold = "flexible", predict = FALSE, ranef = TRUE, 
                    engine = "stan", ...) {
  dots <- list(...)  
  ef <- extract.effects(formula = formula, partial = partial) 
  data <- model.frame(ef$all, data = data, drop.unused.levels = TRUE)
  
  family <- family[1]
  if (is.element(engine,c("stan","jags"))) stan <- engine == "stan"
  else stop("engine must be either stan or jags")
  is.lin <- family %in% c("gaussian", "student", "cauchy")
  is.ord <- family  %in% c("cumulative","cratio","sratio","acat")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  if (!(is.lin | is.ord | is.skew | family %in% 
      c("poisson", "negbinomial", "geometric", "binomial","bernoulli", "categorical")))
    stop(paste(family,"is not a valid family"))
  
  f <- colnames(brm.model.matrix(ef$fixed,data, rm.int = is.ord))
  r <- lapply(lapply(ef$random, brm.model.matrix, data=data, rm.int = is.ord & !stan), colnames)
  p <- colnames(brm.model.matrix(partial, data, rm.int = TRUE))
  out = NULL
  if (is.ord & threshold == "flexible") out <- c(out, "b_Intercept")
  if (is.ord & threshold == "equidistant") out <- c(out, "b_Intercept1", "delta")
  if (length(f)) out <- c(out, paste0("b_",f))
  if (is.ord & length(p)) out <- c(out, paste0("b_",p))
  if (length(ef$group) & engine == "jags") 
    out <- c(out, paste0("V_",ef$group), paste0("VI_",ef$group))
  else if (length(ef$group) & engine == "stan") {
    out <- c(out, unlist(lapply(1:length(ef$group), function(i) 
      paste0("sd_",ef$g[[i]],"_",r[[i]]))))
    out <- c(out, unlist(lapply(1:length(ef$group), function(i)
      if (length(r[[i]])>1) paste0("cor_",ef$g[[i]],"_", unlist(lapply(2:length(r[[i]]), function(j) 
        lapply(1:(j-1), function(k) paste0(r[[i]][k],"_",r[[i]][j]))))))))
    if (ranef) out <- c(out, paste0("r_",ef$group))
  }  
  if (is.lin & !is(ef$add,"formula")) out <- c(out,"sigma")
  if (family == "student") out <- c(out,"nu")
  if (family %in% c("gamma","weibull","negbinomial")) out <- c(out,"shape")
  if (predict) out <- c(out,"Y_pred")
  return(out)
}

#' Extract required data for \code{brms} models
#'
#' @inheritParams brm
#' @param ... Further arguments for testing purposes only
#' 
#' @return A named list of objects containing the required data to fit a \code{brms} model 
#' 
#' @examples
#' data1 <- brm.data(rating ~ treat + period + carry + (1|subject), 
#'          data = inhaler, family = "cumulative")
#' names(data1)
#' 
#' data2 <- brm.data(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'          data = epilepsy, family = c("poisson", "log"))
#' names(data2)
#'          
#' @export
brm.data <- function(formula, data = NULL, family = c("gaussian", "identity"), prior = list(),
                     partial = NULL, engine = "stan", ...) {
  dots <- list(...)  
  ef <- extract.effects(formula = formula, partial = partial) 
  data <- model.frame(ef$all, data = data, drop.unused.levels = TRUE)
  for (g in ef$group) data[[g]] <- as.numeric(as.factor(data[[g]]))
  
  if (is.element(engine,c("stan", "jags"))) stan <- engine == "stan"
  else stop("engine must be either stan or jags")
  
  link <- brm.link(family)
  family <- family[1]
  is.lin <- family %in% c("gaussian", "student", "cauchy")
  is.ord <- family  %in% c("cumulative","cratio","sratio","acat")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  if (!(is.lin | is.ord | is.skew | family %in% 
      c("poisson", "negbinomial", "geometric", "binomial","bernoulli", "categorical")))
    stop(paste(family,"is not a valid family"))
  
  supl.data <- list(N = nrow(data), Y = model.response(data))
  X <- brm.model.matrix(ef$fixed, data, rm.int = is.ord)
  
  if (is.ord | family == "categorical") {
    if (is.factor(supl.data$Y)) supl.data$Y <- as.numeric(supl.data$Y)
    else supl.data$Y <- supl.data$Y - min(supl.data$Y) + 1
  }
  else if (is.factor(supl.data$Y)) 
    stop(paste("family", family, "expects numeric response variable")) 
  if (is.lin) {
    if (is(ef$add, "formula") & length(all.vars(ef$add)) == 1)
      supl.data <- c(supl.data,list(sigma = brm.model.matrix(ef$add, data, rm.int = TRUE)[,1]))
    else if (is(ef$add2, "formula")) {
      inv_weights <- 1/brm.model.matrix(ef$add2, data, rm.int = TRUE)[,1]
      inv_weights <- supl.data$N*inv_weights/sum(inv_weights)
      supl.data <- c(supl.data, list(inv_weights = inv_weights))
    }
  }  
  else if (is.ord | is.element(family, c("binomial", "categorical"))) {
    if (!length(ef$add)) supl.data$max_obs <- max(supl.data$Y)
    else if (is.numeric(ef$add)) supl.data$max_obs <- ef$add
    else if (is(ef$add, "formula") & length(all.vars(ef$add)) == 1) 
      supl.data$max_obs <- brm.model.matrix(ef$add, data, rm.int = TRUE)[,1]
    else stop("Response part of formula is invalid")
    
    if (is(partial,"formula")) {
      if (family %in% c("sratio","cratio","acat")) {
        Xp <- brm.model.matrix(partial, data, rm.int = TRUE)
        supl.data <- c(supl.data, list(Kp = ncol(Xp), Xp = Xp))
        fp <- intersect(colnames(X), colnames(Xp))
        if (length(fp))
          stop(paste("Variables cannot be modeled as fixed and partial effects at the same time.",
                     "Error occured for variables:", paste(fp, collapse = ", ")))
      } 
      else stop("partial is only meaningful for families 'sratio', 'cratio', and 'acat'")  
    } 
    two.cat <- (family == "categorical" | is.ord) & max(supl.data$max_obs) == 2
    if (two.cat) supl.data$Y <- supl.data$Y - 1
    family <- ifelse(family == "binomial" & max(supl.data$Y) == 1 | two.cat, 
                     "bernoulli", family)
  } 
  else if (family %in% c("gamma","exponential","weibull")) {
    if (is(ef$add,"formula") & length(all.vars(ef$add)) == 1) {
      cens <- brm.model.matrix(ef$add, data, rm.int = TRUE)[,1]
      supl.data <- c(supl.data,list(cens = ifelse(cens, 1, 0)))
    }
  }
  
  if (length(ef$random)) {
    Z <- lapply(ef$random, brm.model.matrix, data = data, rm.int = is.ord & !stan)
    r <- lapply(Z, colnames)
    if (family != "categorical")
      to.zero <- unlist(lapply(unlist(lapply(r, intersect, y = colnames(X))), 
                               function(x) which(x == colnames(X)))) 
    else to.zero <- NULL
    X[,to.zero] <- 0
    ncolZ <- lapply(Z, ncol)
    expr <- expression(get(ef$group[[i]], data), length(unique(get(ef$group[[i]], data))), 
                       ncolZ[[i]], Z[[i]])
    for (i in 1:length(ef$random)) {
      name <- paste0(c("", "N_", "K_", "Z_"), ef$group[[i]])
      if (ncolZ[[i]] == 1 & stan) Z[[i]] <- as.vector(Z[[i]])
      for ( j in 1:length(name)) supl.data <- c(supl.data, setNames(list(eval(expr[j])), name[j]))
      if (is.null(dots$Sigma[[ef$group[[i]]]])) mat <- diag(1,ncolZ[[i]])
      else mat <- dots$Sigma[[ef$group[[i]]]]
      if (ncolZ[[i]] > 1) supl.data <- c(supl.data, setNames(list(mat), paste0("Sigma_", ef$group[[i]]))) 
    } 
  }
  if (family == "categorical") {
    supl.data <- c(supl.data, list(Kp = ncol(X), Xp = X))
    X <- data.frame()
  }
  if (stan & ncol(X) > 0) supl.data <- c(supl.data, list(K = ncol(X), X = X))
  else if (!stan & ncol(X) > 0) supl.data <- c(supl.data, list(X = X))
  supl.data
}  
