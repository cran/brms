# This file contains functions dealing with the extended 
# lme4-like formula syntax to specify group-level terms

#' Set up basic grouping terms in \pkg{brms}
#' 
#' Function used to set up a basic grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#' \code{gr} is called implicitly inside the package
#' and there is usually no need to call it directly.
#' 
#' @param ... One or more terms containing grouping factors.
#' @param by An optional factor variable, specifying sub-populations
#'   of the groups. For each level of the \code{by} variable, 
#'   a separate variance-covariance matrix will be fitted. 
#'   Levels of the grouping factor must be nested in levels 
#'   of the \code{by} variable.
#' @param dist Name of the distribution of the group-level effects.
#'   Currently \code{"gaussian"} is the only option.
#' 
#' @seealso \code{\link{brmsformula}}
#' 
#' @examples 
#' \dontrun{
#' # model using basic lme4-style formula
#' fit1 <- brm(count ~ Trt + (1|patient), data = epilepsy)
#' summary(fit1)
#' 
#' # equivalent model using 'gr' which is called anyway internally
#' fit2 <- brm(count ~ Trt + (1|gr(patient)), data = epilepsy)
#' summary(fit2)
#' 
#' # include Trt as a by variable
#' fit3 <- brm(count ~ Trt + (1|gr(patient, by = Trt)), data = epilepsy)
#' summary(fit3)
#' }
#' 
#' @export
gr <- function(..., by = NULL, dist = "gaussian") {
  label <- deparse(match.call())
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) > 1L) {
    stop2("Grouping structure 'gr' expects only a single grouping term")
  }
  stopif_illegal_group(groups[1])
  by <- substitute(by)
  if (!is.null(by)) {
    by <- all.vars(by)
    if (length(by) != 1L) {
      stop2("Argument 'by' must contain exactly one variable.")
    }
  } else {
    by <- ""
  }
  dist <- match.arg(dist, c("gaussian", "student"))
  allvars <- str2formula(c(groups, by))
  nlist(groups, allvars, label, by, dist, type = "")
}

#' Set up multi-membership grouping terms in \pkg{brms}
#' 
#' Function to set up a multi-membership grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#'
#' @inheritParams gr
#' @param weights A matrix specifying the weights of each member.
#'  It should have as many columns as grouping terms specified in \code{...}.
#'  If \code{NULL} (the default), equally weights are used. 
#' @param scale Logical; if \code{TRUE} (the default), 
#'  weights are standardized in order to sum to one per row.
#'  If negative weights are specified, \code{scale} needs
#'  to be set to \code{FALSE}.
#'  
#' @seealso \code{\link{brmsformula}}, \code{\link{mmc}}
#'  
#' @examples 
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(
#'  y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100),
#'  g1 = sample(1:10, 100, TRUE), g2 = sample(1:10, 100, TRUE)
#' )
#' 
#' # multi-membership model with two members per group and equal weights
#' fit1 <- brm(y ~ x1 + (1|mm(g1, g2)), data = dat)
#' summary(fit1)
#' 
#' # weight the first member two times for than the second member
#' dat$w1 <- rep(2, 100)
#' dat$w2 <- rep(1, 100)
#' fit2 <- brm(y ~ x1 + (1|mm(g1, g2, weights = cbind(w1, w2))), data = dat)
#' summary(fit2)
#' 
#' # multi-membership model with level specific covariate values
#' dat$xc <- (dat$x1 + dat$x2) / 2
#' fit3 <- brm(y ~ xc + (1 + mmc(x1, x2) | mm(g1, g2)), data = dat)
#' summary(fit3)
#' }
#'   
#' @export
mm <- function(..., weights = NULL, scale = TRUE, dist = "gaussian") {
  label <- deparse(match.call())
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) < 2) {
    stop2("Multi-membership terms require at least two grouping variables.")
  }
  for (i in seq_along(groups)) {
    stopif_illegal_group(groups[i])
  }
  dist <- match.arg(dist, c("gaussian", "student"))
  scale <- as_one_logical(scale)
  weights <- substitute(weights)
  weightvars <- all.vars(weights)
  allvars <- str2formula(c(groups, weightvars))
  if (!is.null(weights)) {
    weights <- str2formula(deparse_no_string(weights))
    attr(weights, "scale") <- scale
    weightvars <- str2formula(weightvars)
  }
  nlist(
    groups, weights, weightvars, allvars, 
    label, by = "", dist, type = "mm"
  )
}

#' Multi-Membership Covariates
#' 
#' Specify covarariates that vary over different levels 
#' of multi-membership grouping factors thus requiring
#' special treatment. This function is almost solely useful,
#' when called in combination with \code{\link{mm}}. 
#' Outside of multi-membership terms it will behave
#' very much like \code{\link{cbind}}.
#' 
#' @param ... One or more terms containing covariates 
#'   corresponding to the grouping levels specified in \code{\link{mm}}.
#' 
#' @return A matrix with covariates as columns.
#' 
#' @seealso \code{\link{mm}}
#' 
#' @examples 
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(
#'   y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100),
#'   g1 = sample(1:10, 100, TRUE), g2 = sample(1:10, 100, TRUE)
#' )
#' 
#' # multi-membership model with level specific covariate values
#' dat$xc <- (dat$x1 + dat$x2) / 2 
#' fit <- brm(y ~ xc + (1 + mmc(x1, x2) | mm(g1, g2)), data = dat)
#' summary(fit)
#' }
#' 
#' @export
mmc <- function(...) {
  dots <- list(...)
  if (any(ulapply(dots, is_like_factor))) {
    stop2("'mmc' requires numeric variables.")
  }
  out <- cbind(...)
  colnames(out) <- paste0("?", colnames(out))
  out
}

# check if the group part of a group-level term is invalid
# @param group the group part of a group-level term
illegal_group_expr <- function(group) {
  group <- as_one_character(group)
  valid_expr <- ":|([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*"
  rsv_signs <- c("+", "-", "*", "/", "|", "::")
  nzchar(gsub(valid_expr, "", group)) ||
    any(ulapply(rsv_signs, grepl, x = group, fixed = TRUE))
}

stopif_illegal_group <- function(group) {
  if (illegal_group_expr(group)) {
    stop2(
      "Illegal grouping term '", group, "'. It may contain ",
      "only variable names combined by the symbol ':'"
    )
  }
  invisible(NULL)
}

re_lhs <- function(re_terms) {
  get_matches("^[^\\|]*", re_terms) 
}

re_mid <- function(re_terms) {
  get_matches("\\|([^\\|]*\\||)", re_terms)
}

re_rhs <- function(re_terms) {
  sub("^\\|", "", get_matches("\\|[^\\|]*$", re_terms))
}

# extract the three parts of group-level terms
# @param re_terms character vector of RE terms in lme4 syntax
# @return a data.frame with one row per group-level term
re_parts <- function(re_terms) {
  lhs <- re_lhs(re_terms)
  mid <- re_mid(re_terms)
  rhs <- re_rhs(re_terms)
  out <- nlist(lhs, mid, rhs)
  if (any(lengths(out) != length(re_terms))) {
    stop2("Invalid syntax used in group-level terms.")
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# split nested group-level terms and check for special effects terms
# @param re_terms character vector of RE terms in lme4 syntax
split_re_terms <- function(re_terms) {
  if (!length(re_terms)) {
    return(re_terms)
  }
  stopifnot(is.character(re_terms))
  
  # split after grouping factor terms
  re_parts <- re_parts(re_terms)
  new_re_terms <- vector("list", length(re_terms))
  for (i in seq_along(re_terms)) {
    new_re_rhs <- terms(formula(paste0("~", re_parts$rhs[i])))
    new_re_rhs <- attr(new_re_rhs, "term.labels")
    new_re_rhs <- ifelse(
      !grepl("^(gr|mm)\\(", new_re_rhs), 
      paste0("gr(", new_re_rhs, ")"), new_re_rhs
    )
    new_re_terms[[i]] <- paste0(
      re_parts$lhs[i], re_parts$mid[i], new_re_rhs
    )
  }
  re_terms <- unlist(new_re_terms)
  
  # split after coefficient types
  re_parts <- re_parts(re_terms)
  new_re_terms <- type <- vector("list", length(re_terms))
  for (i in seq_along(re_terms)) {
    lhs_form <- formula(paste("~", re_parts$lhs[i]))
    lhs_all_terms <- all_terms(lhs_form)
    # otherwise varying intercepts cannot be handled reliably
    is_cs_term <- grepl_expr(regex_sp("cs"), lhs_all_terms)
    if (any(is_cs_term) && !all(is_cs_term)) {
      stop2("Please specify category specific effects ",
            "in separate group-level terms.")
    }
    new_lhs <- NULL
    # prepare effects of special terms
    valid_types <- c("sp", "cs", "mmc")
    invalid_types <- c("sm", "gp")
    for (t in c(valid_types, invalid_types)) {
      lhs_tform <- do_call(paste0("parse_", t), list(lhs_form))
      if (is.formula(lhs_tform)) {
        if (t %in% invalid_types) {
          stop2("Cannot handle splines or GPs in group-level terms.")
        }
        new_lhs <- c(new_lhs, formula2str(lhs_tform, rm = 1))
        type[[i]] <- c(type[[i]], t)
      }
    }
    # prepare effects of basic terms
    fe_form <- parse_fe(lhs_form)
    fe_terms <- all_terms(fe_form)
    has_intercept <- attr(terms(fe_form), "intercept")
    if (length(fe_terms) || has_intercept && !"cs" %in% type[[i]]) {
      new_lhs <- c(new_lhs, formula2str(fe_form, rm = 1))
      type[[i]] <- c(type[[i]], "")
    }
    if (length(new_lhs) > 1 && re_parts$mid[i] != "||") {
      id <- gsub("\\|", "", re_parts$mid[i])
      if (!nzchar(id)) {
        # ID is required to model coefficients as correlated 
        # if multiple types are provided within the same term
        id <- collapse(sample(0:9, 10, TRUE))
        re_parts$mid[i] <- paste0("|", id, "|")
      }
    }
    new_re_terms[[i]] <- paste0(new_lhs, re_parts$mid[i], re_parts$rhs[i])
    new_re_terms[[i]] <- new_re_terms[[i]][order(type[[i]])]
    type[[i]] <- sort(type[[i]])
  }
  re_terms <- unlist(new_re_terms)
  structure(re_terms, type = unlist(type))
}

# extract group-level terms from a formula of character vector
# @param x formula or character vector
# @param formula return a formula rather than a character string?
# @param brackets include group-level terms in brackets?
get_re_terms <- function(x, formula = FALSE, brackets = TRUE) {
  if (is.formula(x)) {
    x <- all_terms(x)
  }
  re_pos <- grepl("\\|", x)
  out <- x[re_pos]
  if (brackets && length(out)) {
    out <- paste0("(", out, ")")
  } 
  if (formula) {
    if (length(out)) {
      out <- formula(paste("~ 1", collapse("+", out)))
    } else {
      out <- ~ 1
    }
  }
  out
}

# validate the re_formula argument
# @inheritParams extract_draws.brmsfit
# @param formula: formula to match re_formula with
# @return updated re_formula containing only terms existent in formula
check_re_formula <- function(re_formula, formula) {
  old_re_formula <- get_re_terms(formula, formula = TRUE)
  if (is.null(re_formula)) {
    re_formula <- old_re_formula
  } else if (SW(anyNA(re_formula))) {
    re_formula <- ~ 1
  } else {
    re_formula <- get_re_terms(as.formula(re_formula), formula = TRUE)
    new <- parse_bf(re_formula, check_response = FALSE)$dpars$mu$re
    old <- parse_bf(old_re_formula, check_response = FALSE)$dpars$mu$re
    if (nrow(new)) {
      new_terms <- lapply(new$form, terms)
      found <- rep(FALSE, nrow(new))
      for (i in 1:nrow(new)) {
        group <- new$group[[i]]
        old_terms <- lapply(old$form[old$group == group], terms)
        j <- 1
        while (!found[i] && j <= length(old_terms)) {
          found[i] <- isTRUE(
            all(attr(new_terms[[i]], "term.labels") %in% 
                  attr(old_terms[[j]], "term.labels")) &&
              attr(new_terms[[i]], "intercept") <=
              attr(old_terms[[j]], "intercept")
          )
          j <- j + 1
        }
      }  
      new <- new[found, ]
      if (nrow(new)) {
        forms <- ulapply(new$form, formula2str, rm = 1)
        groups <- ulapply(new$gcall, "[[", "label")
        re_terms <- paste("(", forms, "|", groups, ")")
        re_formula <- formula(paste("~", paste(re_terms, collapse = "+")))
      } else {
        re_formula <- ~ 1
      }
    } else {
      re_formula <- ~ 1
    }
  }
  re_formula
}

# remove existing group-level terms in formula and
# add valid group-level terms of re_formula
update_re_terms <- function(formula, re_formula) {
  UseMethod("update_re_terms")
}

#' @export
update_re_terms.mvbrmsformula <- function(formula, re_formula) {
  formula$forms <- lapply(formula$forms, update_re_terms, re_formula)
  formula
}

#' @export
update_re_terms.brmsformula <- function(formula, re_formula) {
  formula$formula <- update_re_terms(formula$formula, re_formula)
  formula$pforms <- lapply(formula$pforms, update_re_terms, re_formula)
  formula
}

#' @export
update_re_terms.formula <- function(formula, re_formula = NULL) {
  if (is.null(re_formula) || get_nl(formula)) {
    return(formula)
  }
  re_formula <- check_re_formula(re_formula, formula)
  new_formula <- formula2str(formula)
  old_re_terms <- get_re_terms(formula, brackets = FALSE)
  if (length(old_re_terms)) {
    # remove old group-level terms
    rm_terms <- c(
      paste0("+(", old_re_terms, ")"),
      paste0("(", old_re_terms, ")"),
      old_re_terms
    )
    new_formula <- rename(new_formula, rm_terms, "")
    if (grepl("~\\+*$", new_formula)) {
      # lhs only formulas are syntactically invalid
      # also check for trailing '+' signs (#769)
      new_formula <- paste(new_formula, "1")
    }
  }
  # add new group-level terms
  new_re_terms <- get_re_terms(re_formula)
  new_formula <- paste(c(new_formula, new_re_terms), collapse = "+")
  new_formula <- formula(new_formula)
  attributes(new_formula) <- attributes(formula)
  new_formula
}

# extract group-level terms
get_re <- function(x, ...) {
  UseMethod("get_re")
}

#' @export
get_re.default <- function(x, ...) {
  NULL
}

# get group-level information in a data.frame
# @param bterms object of class 'brmsterms'
# @param all logical; include ranefs of additional parameters?
#' @export
get_re.brmsterms <- function(x, all = TRUE, ...) {
  if (all) {
    re <- named_list(c(names(x$dpars), names(x$nlpars)))
    for (dp in names(x$dpars)) {
      re[[dp]] <- get_re(x$dpars[[dp]])
    }
    for (nlp in names(x$nlpars)) {
      re[[nlp]] <- get_re(x$nlpars[[nlp]])
    }
    re <- do_call(rbind, re)
  } else {
    x$dpars[["mu"]]$nlpars <- NULL
    re <- get_re(x$dpars[["mu"]])
  }
  re
}

#' @export
get_re.mvbrmsterms <- function(x, ...) {
  do_call(rbind, lapply(x$terms, get_re, ...))
}

#' @export
get_re.btl <- function(x, ...) {
  stopifnot(is.data.frame(x$re))
  px <- check_prefix(x)
  re <- x$re
  re$resp <- rep(px$resp, nrow(re)) 
  re$dpar <- rep(px$dpar, nrow(re))
  re$nlpar <- rep(px$nlpar, nrow(re)) 
  re
}

# gather information on group-level effects
# @param bterms object of class brmsterms
# @param data data.frame containing all model variables
# @param all include REs of all parameters?
# @param old_levels optional original levels of the grouping factors
# @return a tidy data.frame with the following columns:
#   id: ID of the group-level effect 
#   group: name of the grouping factor
#   gn: number of the grouping term within the respective formula
#   coef: name of the group-level effect
#   cn: number of the effect within the ID
#   resp: name of the response variable
#   dpar: name of the distributional parameter
#   nlpar: name of the non-linear parameter
#   cor: are correlations modeled for this effect?
#   ggn: global number of the grouping factor
#   type: special effects type; can be 'sp' or 'cs'
#   gcall: output of functions 'gr' or 'mm'
#   form: formula used to compute the effects
tidy_ranef <- function(bterms, data, all = TRUE, old_levels = NULL) {
  data <- combine_groups(data, get_group_vars(bterms))
  re <- get_re(bterms, all = all)
  ranef <- vector("list", nrow(re))
  used_ids <- new_ids <- NULL
  id_groups <- list()
  j <- 1
  for (i in seq_rows(re)) {
    if (!nzchar(re$type[i])) {
      coef <- colnames(get_model_matrix(re$form[[i]], data)) 
    } else if (re$type[i] == "sp") {
      coef <- tidy_spef(re$form[[i]], data)$coef
    } else if (re$type[i] == "mmc") {
      coef <- rename(all_terms(re$form[[i]]))
    } else if (re$type[i] == "cs") {
      resp <- re$resp[i]
      if (nzchar(resp)) {
        stopifnot(is.mvbrmsterms(bterms))
        nthres <- max(get_thres(bterms$terms[[resp]]))
      } else {
        stopifnot(is.brmsterms(bterms))
        nthres <- max(get_thres(bterms))
      }
      indices <- paste0("[", seq_len(nthres), "]")
      coef <- colnames(get_model_matrix(re$form[[i]], data = data))
      coef <- as.vector(t(outer(coef, indices, paste0)))
    }
    avoid_dpars(coef, bterms = bterms)
    rdat <- data.frame(
      id = re$id[[i]],
      group = re$group[[i]],
      gn = re$gn[[i]],
      gtype = re$gtype[[i]],
      coef = coef, 
      cn = NA,
      resp = re$resp[[i]],
      dpar = re$dpar[[i]],
      nlpar = re$nlpar[[i]],
      ggn = NA,
      cor = re$cor[[i]],
      type = re$type[[i]],
      by = re$gcall[[i]]$by,
      dist = re$gcall[[i]]$dist,
      stringsAsFactors = FALSE
    )
    bylevels <- NULL
    if (nzchar(rdat$by[1])) {
      bylevels <- rm_wsp(levels(factor(get(rdat$by[1], data))))
    }
    rdat$bylevels <- repl(bylevels, nrow(rdat))
    rdat$form <- repl(re$form[[i]], nrow(rdat))
    rdat$gcall <- repl(re$gcall[[i]], nrow(rdat)) 
    # prepare group-level IDs
    id <- re$id[[i]]
    if (is.na(id)) {
      rdat$id <- j
      j <- j + 1
    } else {
      if (id %in% used_ids) {
        k <- match(id, used_ids)
        rdat$id <- new_ids[k]
        new_id_groups <- c(re$group[[i]], re$gcall[[i]]$groups)
        if (!identical(new_id_groups, id_groups[[k]])) {
          stop2("Can only combine group-level terms of the ",
                "same grouping factors.")
        }
      } else {
        used_ids <- c(used_ids, id)
        k <- length(used_ids)
        rdat$id <- new_ids[k] <- j
        id_groups[[k]] <- c(re$group[[i]], re$gcall[[i]]$groups)
        j <- j + 1
      }
    }
    ranef[[i]] <- rdat 
  }
  ranef <- do_call(rbind, c(list(empty_ranef()), ranef))
  # check for overlap between different group types
  rsv_groups <- ranef[nzchar(ranef$gtype), "group"]
  other_groups <- ranef[!nzchar(ranef$gtype), "group"]
  inv_groups <- intersect(rsv_groups, other_groups)
  if (length(inv_groups)) {
    inv_groups <- paste0("'", inv_groups, "'", collapse = ", ")
    stop2("Grouping factor names ", inv_groups, " are resevered.")
  }
  # check for duplicated and thus not identified effects
  dup <- duplicated(ranef[, c("group", "coef", vars_prefix())])
  if (any(dup)) {
    dr <- ranef[which(dup)[1], ]
    stop2(
      "Duplicated group-level effects are not allowed.\n",
      "Occured for effect '", dr$coef, "' of group '", dr$group, "'."
    )
  }
  if (nrow(ranef)) {
    for (id in unique(ranef$id)) {
      ranef$cn[ranef$id == id] <- seq_len(sum(ranef$id == id))
    }
    ranef$ggn <- match(ranef$group, unique(ranef$group))
    if (is.null(old_levels)) {
      rsub <- ranef[!duplicated(ranef$group), ]
      levels <- named_list(rsub$group)
      for (i in seq_along(levels)) {
        # combine levels of all grouping factors within one grouping term
        levels[[i]] <- unique(ulapply(
          rsub$gcall[[i]]$groups, 
          function(g) levels(factor(get(g, data)))
        ))
        # store information of corresponding by levels
        if (nzchar(rsub$by[i])) {
          stopifnot(!nzchar(rsub$type[i]))
          by <- rsub$by[i]
          bylevels <- rsub$bylevels[[i]]
          g <- rsub$gcall[[i]]$groups
          J <- match(get(g, data), levels[[i]])
          df <- unique(data.frame(J, by = rm_wsp(get(by, data))))
          if (nrow(df) > length(unique(J))) {
            stop2("Some levels of '", g, "' correspond ", 
                  "to multiple levels of '", by, "'.")
          }
          df <- df[order(df$J), ]
          by_per_level <- bylevels[match(df$by, bylevels)]
          attr(levels[[i]], "by") <- by_per_level
        }
      }
      attr(ranef, "levels") <- levels 
    } else {
      # for newdata numeration has to depend on the original levels
      attr(ranef, "levels") <- old_levels
    }
  }
  structure(ranef, class = c("ranef_frame", "data.frame"))
}

empty_ranef <- function() {
  structure(
    data.frame(
      id = numeric(0), group = character(0), gn = numeric(0),
      coef = character(0), cn = numeric(0), resp = character(0),
      dpar = character(0), nlpar = character(0), ggn = numeric(0),
      cor = logical(0), type = character(0), form = character(0), 
      stringsAsFactors = FALSE
    ),
    class = c("ranef_frame", "data.frame")
  )
}

is.ranef_frame <- function(x) {
  inherits(x, "ranef_frame")
}

# extract names of all grouping variables
get_group_vars <- function(x, ...) {
  UseMethod("get_group_vars") 
}

#' @export
get_group_vars.brmsfit <- function(x, ...) {
  get_group_vars(x$formula, ...)
}

#' @export
get_group_vars.default <- function(x, ...) {
  get_group_vars(parse_bf(x), ...)
}

#' @export
get_group_vars.brmsterms <- function(x, ...) {
  .get_group_vars(x, ...)
}

#' @export
get_group_vars.mvbrmsterms <- function(x, ...) {
  .get_group_vars(x, ...)
}

.get_group_vars <- function(x, ...) {
  out <- c(get_re_groups(x), get_me_groups(x), get_ac_groups(x))
  out <- out[nzchar(out)]
  if (length(out)) {
    c(out) <- unlist(strsplit(out, ":"))
    out <- sort(unique(out))
  }
  out
}

# get names of grouping variables of re terms
get_re_groups <- function(x, ...) {
  ulapply(get_re(x)$gcall, "[[", "groups")
}

# extract information about groups with a certain distribution
get_dist_groups <- function(ranef, dist) {
  out <- subset2(ranef, dist = dist)
  out[!duplicated(out$group), c("group", "ggn", "id")]
}

# extract list of levels with one element per grouping factor
# @param ... objects with a level attribute
get_levels <- function(...) {
  dots <- list(...)
  out <- vector("list", length(dots))
  for (i in seq_along(out)) {
    levels <- attr(dots[[i]], "levels", exact = TRUE)
    if (is.list(levels)) {
      stopifnot(!is.null(names(levels)))
      out[[i]] <- as.list(levels)
    } else if (!is.null(levels)) {
      stopifnot(isTRUE(nzchar(names(dots)[i])))
      out[[i]] <- setNames(list(levels), names(dots)[[i]])
    }
  }
  out <- unlist(out, recursive = FALSE)
  out[!duplicated(names(out))]
}

# extract names of group-level effects
# @param ranef output of tidy_ranef()
# @param group optinal name of a grouping factor for
#   which to extract effect names
# @param bylevels optional names of 'by' levels for 
#    which to extract effect names
# @return a vector of character strings
get_rnames <- function(ranef, group = NULL, bylevels = NULL) {
  stopifnot(is.data.frame(ranef))
  if (!is.null(group)) {
    group <- as_one_character(group)
    ranef <- subset2(ranef, group = group)
  }
  stopifnot(length(unique(ranef$group)) == 1L)
  out <- paste0(usc(combine_prefix(ranef), "suffix"), ranef$coef)
  if (isTRUE(nzchar(ranef$by[1]))) {
    if (!is.null(bylevels)) {
      stopifnot(all(bylevels %in% ranef$bylevels[[1]]))
    } else {
      bylevels <- ranef$bylevels[[1]]
    }
    bylabels <- paste0(ranef$by[1], bylevels)
    out <- outer(out, bylabels, paste, sep = ":")
  }
  out
}