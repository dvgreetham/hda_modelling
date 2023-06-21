library(survival)
library(icenReg)
library(dplyr)
library(ggplot2)
library(magrittr)


prepare_factors <- function(df) {
    df$jobclas <- factor(df$jobclas,
                         levels = c("full-time", "part-time", "self-employed", "not working")
    )
    df$gender <- factor(df$gender,
                        levels = c("male", "female")
    )
    df$res <- factor(df$res,
                     levels = c("rural", "village", "town", "city")
    )
    df$maristat <- factor(df$maristat,
                          levels = c("single", "married", "partnership", "divorced", "widowed", "separated")
    )
    df$race <- factor(df$race,
                      levels = c("Causasian", "American Black", "Hispanic", "Native American", "Asian", "Mixed", "Other")
    )
    
    if (is.factor(df$region)) {
        # pass
    } else if (is.character(df$region)) {
        df$region <- factor(df$region,
                            levels = c("Europe", "Northern America", "Australasia", "Latin America")
        )
    } else {
        df$region <- factor(df$region,
                            levels = c("0", "1", "2"),
                            labels = c("Europe", "Northern America", "Other")
        )
    }
    
    df
}


prepare_survival <- function(df) {
    # Duration of stage is in the interval duration_min to duration_max
    df$SurvObj <- with(df, Surv(duration_min, duration_max, type = "interval2"))
    
    df
}


### Extract and plot NPMLE data from icenReg fits
npmle_table <- function (x, ...) {
    UseMethod("npmle_table", x)
}


npmle_table.sp_curves <- function (x, group=NULL) {
    dplyr::bind_cols(x$Tbull_ints, x$S_curves)
}


npmle_table.icenReg_fit <- function (x, group=NULL) {
    getSCurves(x) %>% npmle_table
}


npmle_table.ic_npList <- function (x, group="group") {
    dplyr::bind_rows(
        lapply(x$fitList, npmle_table),
        .id=group
    )
}


## Survival plot for non-parametric fit to interval censored data
# WIP: alternative to survminer::ggsurvplot for icenReg fit
plot_np <- function(np_fit, group="group") {
    data <- npmle_table(np_fit, group=group)
    
    if (group %in% names(data))
        col = sym(group)
    else
        col = NULL
    
    # TODO: data needs to be padded so that complete boxes are drawn
    # See lines.sp_curves
    
    # Derived columns to support plotting
    data$baseline_prev <- c(data$baseline[1], head(data$baseline, -1))
    
    ggplot() +
        geom_step(data=data, mapping=aes(x=lower, y=baseline, col=!!col)) +
        geom_step(data=data, mapping=aes(x=upper, y=baseline, col=!!col)) +
        geom_rect(data=data, aes(xmin=lower, xmax=upper, ymin=baseline, ymax=baseline_prev, fill=!!col), alpha=0.5) +
        xlab("Time") +
        ylab("Survival probability")
}


## Forest plot
# This is a modified version of survminer::ggforest
# that supports icenReg ic_ph fits.
plot_forest <- function (x, ...) {
    UseMethod("plot_forest", x)
}


plot_forest.coxph <- function(model, data = NULL, ...) {
    stopifnot(inherits(model, "coxph"))
    
    # get data and variables/terms from cox model
    if(is.null(data)){
        data <- eval(model$call$data)
    }
    terms <- names(attr(model$terms, "dataClasses")[-1])
    
    # use broom to get some required statistics
    coef <- as.data.frame(broom::tidy(model, conf.int = TRUE))
    gmodel <- broom::glance(model)
    
    comment <- paste0("# Events: ", gmodel$nevent, "; Global p-value (Log-Rank): ",
                      format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", round(gmodel$AIC,2),
                      "; Concordance Index: ", round(gmodel$concordance,2))
    
    plot_forest_data(terms, data, coef, comment=comment, ...)
}


plot_forest.ic_ph <- function(model, data = NULL, ...) {
    stopifnot(inherits(model, "ic_ph"))
    
    # get data and variables/terms from cox model
    if (missing(data)) {
        data <- model$getRawData()
        if(is.null(data))
            stop('Could not find data from fit. Original model must be built with data argument (rather than variables found in the Global Environment) supplied to be retreivable')
    }
    terms <- names(model$xlevels)
    
    # get some required statistics
    coef <- data.frame(summary(model)$summaryParameters) %>%
        dplyr::rename(estimate = Estimate, p.value = p, std.error = Std.Error, statistic = z.value) %>%
        dplyr::select(-Exp.Est.)
    coef$conf.low <- apply(model$bsMat, 2, quantile, prob=0.025)
    coef$conf.high <- apply(model$bsMat, 2, quantile, prob=0.975)
    coef <- coef %>% tibble::rownames_to_column("term")
    
    # TODO: Add useful comments for model fit properties
    comment <- NULL
    
    plot_forest_data(terms, data, coef, comment=comment, ...)
}


plot_forest_data <- function(terms, data, coef, comment = NULL,
                             main = "Hazard ratio", cpositions = c(0.02, 0.22, 0.4),
                             fontsize = 0.7, refLabel = "reference", noDigits = 2
) {
    conf.high <- conf.low <- estimate <- NULL
    
    # extract statistics for every variable
    allTerms <- lapply(seq_along(terms), function(i){
        var <- terms[i]
        values <- data[[var]] %>% unlist()
        if (is.factor(values) || is.character(values)) {
            adf <- as.data.frame(table(values))
            cbind(var = var, adf, pos = 1:nrow(adf))
        }
        else if (is.numeric(values)) {
            data.frame(var = var, Var1 = "", Freq = nrow(data),
                       pos = 1)
        }
        else {
            vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
            data.frame(var = vars, Var1 = "", Freq = nrow(data),
                       pos = seq_along(vars))
        }
    })
    allTermsDF <- do.call(rbind, allTerms)
    colnames(allTermsDF) <- c("var", "level", "N", "pos")
    inds <- apply(allTermsDF[,1:2], 1, paste0, collapse="")
    
    # use broom again to get remaining required statistics
    rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
    toShow <- cbind(allTermsDF, coef[inds,])[,c("var", "level", "N", "p.value", "estimate", "conf.low", "conf.high", "pos")]
    toShowExp <- toShow[,5:7]
    toShowExp[is.na(toShowExp)] <- 0
    toShowExp <- format(exp(toShowExp), digits=noDigits)
    toShowExpClean <- data.frame(toShow,
                                 pvalue = signif(toShow[,4],noDigits+1),
                                 toShowExp)
    toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, noDigits+1), " ",
                                   ifelse(toShowExpClean$p.value < 0.05, "*",""),
                                   ifelse(toShowExpClean$p.value < 0.01, "*",""),
                                   ifelse(toShowExpClean$p.value < 0.001, "*",""))
    toShowExpClean$ci <- paste0("(",toShowExpClean[,"conf.low.1"]," - ",toShowExpClean[,"conf.high.1"],")")
    toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
    toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
    toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
    toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
    toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
    toShowExpClean$var = as.character(toShowExpClean$var)
    toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
    # make label strings:
    toShowExpClean$N <- paste0("(N=",toShowExpClean$N,")")
    
    #flip order
    toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
    
    rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = TRUE)
    breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
    rangeplot <- rangeb
    # make plot twice as wide as needed to create space for annotations
    rangeplot[1] <- rangeplot[1] - diff(rangeb)
    # increase white space on right for p-vals:
    rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
    
    width <- diff(rangeplot)
    # y-coordinates for labels:
    y_variable <- rangeplot[1] +  cpositions[1] * width
    y_nlevel <- rangeplot[1]  +  cpositions[2] * width
    y_cistring <- rangeplot[1]  +  cpositions[3] * width
    y_stars <- rangeb[2]
    x_annotate <- seq_len(nrow(toShowExpClean))
    
    # geom_text fontsize is in mm (https://github.com/tidyverse/ggplot2/issues/1828)
    annot_size_mm <- fontsize *
        as.numeric(grid::convertX(unit(theme_get()$text$size, "pt"), "mm"))
    
    p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
        geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                      ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                      fill = ordered(seq_along(var) %% 2 + 1))) +
        scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
        geom_point(pch = 15, size = 4) +
        geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), width = 0.15) +
        geom_hline(yintercept = 1, linetype = 3) +
        coord_flip(ylim = exp(rangeplot)) +
        ggtitle(main) +
        scale_y_log10(
            name = "",
            labels = sprintf("%g", breaks),
            expand = c(0.02, 0.02),
            breaks = breaks) +
        theme_light() +
        theme(panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              legend.position = "none",
              panel.border=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        xlab("") +
        annotate(geom = "text", x = x_annotate, y = exp(y_variable),
                 label = toShowExpClean$var, fontface = "bold", hjust = 0,
                 size = annot_size_mm) +
        annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0,
                 label = toShowExpClean$level, vjust = -0.1, size = annot_size_mm) +
        annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
                 label = toShowExpClean$N, fontface = "italic", hjust = 0,
                 vjust = ifelse(toShowExpClean$level == "", .5, 1.1),
                 size = annot_size_mm) +
        annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
                 label = toShowExpClean$estimate.1, size = annot_size_mm,
                 vjust = ifelse(toShowExpClean$estimate.1 == "reference", .5, -0.1)) +
        annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
                 label = toShowExpClean$ci, size = annot_size_mm,
                 vjust = 1.1,  fontface = "italic") +
        annotate(geom = "text", x = x_annotate, y = exp(y_stars),
                 label = toShowExpClean$stars, size = annot_size_mm,
                 hjust = -0.2,  fontface = "italic")
    
    if (!is.null(comment))
        p <- p +
        annotate(geom = "text", x = 0.5, y = exp(y_variable),
                 label = comment,
                 size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")
    
    # switch off clipping for p-vals, bottom annotation:
    gt <- ggplot_gtable(ggplot_build(p))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    # grid.draw(gt)
    # invisible(p)
    ggpubr::as_ggplot(gt)
}


find_turnbull_intervals <- function(sample) {
    # Find Turnbull intervals
    left_edges <- unique(sample$duration_min)
    right_edges <- unique(sample$duration_max)
    df_edges <- data.frame(edges = c(left_edges, right_edges),
                           edge_type = c(rep('left', length(left_edges)),
                                         rep('right', length(right_edges))))
    
    df_edges <- df_edges %>% arrange(edges)
    
    turnbull_left <- numeric(length(left_edges)) + NA
    turnbull_right <- numeric(length(right_edges)) + NA
    
    last_edge_type <- 'left'
    idx_turnbull <- 1
    for (i_edge in seq_len(nrow(df_edges))) {
        if (df_edges[i_edge, 'edge_type'] == 'left') {
            if (last_edge_type == 'left') {
                # Overwrite
                turnbull_left[idx_turnbull] <- df_edges[i_edge, 'edges']
            } else {
                # this is a new interval
                idx_turnbull <- idx_turnbull + 1
                turnbull_left[idx_turnbull] <- df_edges[i_edge, 'edges']
                last_edge_type <- 'left'
            }
        } else {  # edge_type is 'right'
            if (last_edge_type == 'left') {
                turnbull_right[idx_turnbull] <- df_edges[i_edge, 'edges']
                last_edge_type <- 'right'
            }
            # else do nothing -- keep the leftmost right edge
        }
    }
    
    turnbull_left <- turnbull_left[!is.na(turnbull_left)]
    turnbull_right <- turnbull_right[!is.na(turnbull_right)]
    
    return(list(turnbull_left = turnbull_left, turnbull_right = turnbull_right))
}


calculate_turnbull_estimator <- function(time_vec, sample, s_init=NULL,
                                         max_iters=1000, tol=1e-6) {
    
    sample_size <- nrow(sample)
    
    turnbull_intervals <- find_turnbull_intervals(sample)
    turnbull_left <- turnbull_intervals$turnbull_left
    turnbull_right <- turnbull_intervals$turnbull_right
    
    if (is.null(s_init)) {
        s_init <- numeric(length(turnbull_left)) + 1/length(turnbull_left)
    }
    
    s_old <- s_init
    alpha <- calculate_alpha_matrix(sample, turnbull_left, turnbull_right)
    
    for (i_iter in seq_len(max_iters)) {
        v <- alpha %*% s_old
        mu <- s_old * (t(alpha) %*% (1/v))
        s_new <- mu / sum(mu)
        
        s_max_diff <- max(abs(s_new - s_old))
        if (s_max_diff < tol) {
            print(paste0("Number of iterations: ", i_iter+1))
            break
        }
        
        if (i_iter == max_iters) {
            print("Maximum number of iterations reached")
        }
        
        s_old <- s_new
    }
    
    # Construct survival function
    surv_func_upper <- numeric(length(time_vec))
    surv_func_lower <- numeric(length(time_vec))
    surv_func_upper[time_vec <= turnbull_left[1]] <- 1
    surv_func_lower[time_vec <= turnbull_left[1]] <- 1
    surv_func_upper[time_vec > turnbull_right[length(turnbull_right)]] <- 0
    
    for (j_interval in seq_along(turnbull_left)) {
        if (j_interval < length(turnbull_left)) {
            idx_to_right_of_turnbull_interval <- time_vec > turnbull_right[j_interval] &
                time_vec <= turnbull_left[j_interval+1]
        } else {
            idx_to_right_of_turnbull_interval <- time_vec > turnbull_right[j_interval]
        }
        
        value_to_right_of_interval <- 1 - sum(s_new[1:j_interval])
        surv_func_upper[idx_to_right_of_turnbull_interval] <- value_to_right_of_interval
        surv_func_lower[idx_to_right_of_turnbull_interval] <- value_to_right_of_interval
        
        # Values within intervals are arbitrary
        idx_in_turnbull_interval = time_vec > turnbull_left[j_interval] &
            time_vec <= turnbull_right[j_interval]
        if (j_interval == 1) {
            surv_func_upper[idx_in_turnbull_interval] <- 1
        } else {
            surv_func_upper[idx_in_turnbull_interval] <- 1 - sum(s_new[1:(j_interval-1)])  # value to left of interval
        }
        
        surv_func_lower[idx_in_turnbull_interval] <- value_to_right_of_interval
        
    }
    
    surv_func_var <- calculate_survival_variance(sample, turnbull_left, turnbull_right,
                                                 s_new, time_vec, alpha)
    
    return(list(surv_func_lower = surv_func_lower,
                surv_func_upper = surv_func_upper,
                turnbull_estimates = s_new,
                surv_func_var = surv_func_var))
}


calculate_survival_variance <- function(sample, turnbull_left, turnbull_right,
                                        turnbull_estimates, time_vec, alpha) {
    
    surv_func_var <- numeric(length(time_vec))
    sample_size <- nrow(sample)
    
    # Calculate variance estimates
    fisher_infmn <- calculate_turnbull_fisher_information(sample, turnbull_left,
                                                          turnbull_right,
                                                          turnbull_estimates,
                                                          alpha)
    inv_fisher_infmn <- solve(fisher_infmn)
    
    for (j_interval in seq_along(turnbull_left)) {
        if (j_interval < length(turnbull_left)) {
            idx_to_right_of_turnbull_interval <- time_vec > turnbull_right[j_interval] &
                time_vec <= turnbull_left[j_interval+1]
        } else {
            idx_to_right_of_turnbull_interval <- time_vec > turnbull_right[j_interval]
        }
        
        elem_vec <- c(rep(1, j_interval),
                      rep(0, length(turnbull_estimates) - j_interval))
        surv_func_var[idx_to_right_of_turnbull_interval] <- t(elem_vec) %*%
            (inv_fisher_infmn %*% elem_vec) / sample_size
        
        idx_in_turnbull_interval = time_vec > turnbull_left[j_interval] &
            time_vec <= turnbull_right[j_interval]
        
        surv_func_var[idx_in_turnbull_interval] <- NA
    }
    
    return(surv_func_var)
}


calculate_alpha_matrix <- function(sample, turnbull_left, turnbull_right) {
    # Construct the $\alpha$ matrix, which is defined by
    #
    # $$
    # \alpha_{ij} = I\{(p_j, q_j] \subseteq (L_i, U_i]\},
    # $$
    #
    # where $I$ denotes an indicator function, and $(L_i, U_i]$ is the interval in which the event time for subject $i$ is known to lie.
    #
    # This indicates whether it's possible that subject $i$ had an event in the interval $(t_{j-1}, t_j]$.
    
    sample_size <- nrow(sample)
    alpha <- matrix(0, nrow = sample_size, ncol = length(turnbull_left))
    
    for (i_subj in seq_len(sample_size)) {
        for (j_interval in seq_along(turnbull_left)) {
            alpha[i_subj, j_interval] <- (turnbull_left[j_interval] >=
                                              sample[i_subj, 'duration_min']
                                          & turnbull_right[j_interval] <=
                                              sample[i_subj, 'duration_max'])
        }
    }
    
    return(alpha)
}


calculate_turnbull_fisher_information <- function(sample, turnbull_left,
                                                  turnbull_right,
                                                  turnbull_estimates, alpha=NULL) {
    
    if (is.null(alpha)) {
        alpha <- calculate_alpha_matrix(sample, turnbull_left, turnbull_right)
    }
    
    fisher_infmn <- matrix(0, nrow = ncol(alpha), ncol = ncol(alpha))
    
    denom <- (alpha %*% turnbull_estimates)^2
    for (k_row in seq_len(nrow(fisher_infmn))) {
        for (l_col in seq_len(ncol(fisher_infmn))) {
            fisher_infmn[k_row, l_col] <- sum(alpha[, k_row] * alpha[, l_col] / denom)
        }
    }
    
    return(fisher_infmn)
}


run_single_imputation <- function(n_mc_runs, time_vec, sample,
                                  use_exponential_greenwood_CIs = FALSE) {
    all_surv_funcs <- matrix(0, nrow = n_mc_runs, ncol = length(time_vec))
    all_sd <- matrix(0, nrow = n_mc_runs, ncol = length(time_vec))
    idx_observed <- is.finite(sample[['duration_max']])
    n_observed <- sum(idx_observed)
    
    for (i_run in seq_len(n_mc_runs)) {
        # random samples for each interval
        sample_imputed <- sample
        
        sample_imputed[!idx_observed, 'duration_true'] <-
            sample_imputed[!idx_observed, 'duration_min']
        
        high <- sample_imputed[idx_observed, 'duration_max']
        low <- sample_imputed[idx_observed, 'duration_min']
        sample_imputed[idx_observed, 'duration_true'] <-
            runif(n = n_observed) * (high - low) + low
        
        # Kaplan-Meier fit
        sample_imputed$SurvObj <- Surv(sample_imputed$duration_true,
                                       is.finite(sample_imputed$duration_max))
        km_fit <- survfit(SurvObj ~ 1, data = sample_imputed)
        km_fit <- survfit0(km_fit)
        
        # Fill in survival function matrix
        for (i_time in seq_along(km_fit$time)) {
            if (i_time < length(km_fit$time)) {
                idx_assign <- time_vec >= km_fit$time[i_time] & time_vec < km_fit$time[i_time+1]
            } else {
                idx_assign <- time_vec >= km_fit$time[i_time]
            }
            
            all_surv_funcs[i_run, idx_assign] <- km_fit$surv[i_time]
            alpha <- 0.05  # Default
            normal_quantile <- qnorm(alpha/2, lower.tail = FALSE)
            S_hat <- km_fit$surv[i_time]
            
            if (use_exponential_greenwood_CIs) {
                # Go from upper 97.5% point of standard Greenwood formula to standard
                # deviation for Greenwood exponential formula
                all_sd[i_run, idx_assign] <- (km_fit$upper[i_time] - S_hat) /
                    normal_quantile / S_hat / log(S_hat)
            } else {
                # standard greenwood
                all_sd[i_run, idx_assign] <- (km_fit$upper[i_time] - S_hat) /
                    normal_quantile
            }
        }
    }
    
    return(list(all_surv_funcs = all_surv_funcs, all_sd = all_sd))
}


get_rc_km_curve_from_ic_data <- function(sample, time_vec) {
    sample_KM_true <- sample
    sample_KM_true[is.infinite(sample_KM_true$duration_max), 'duration_true'] <-
        sample_KM_true[is.infinite(sample_KM_true$duration_max), 'duration_min']
    sample_KM_true$SurvObj <- Surv(sample_KM_true$duration_true,
                                   is.finite(sample_KM_true$duration_max))
    km_fit_sample_true <- survfit(SurvObj ~ 1, data = sample_KM_true, conf.int = FALSE)
    km_true <- numeric(length(time_vec))
    drop_times <- km_fit_sample_true$time
    surv_sample_values <- km_fit_sample_true$surv
    surv_sample <- numeric(length(time_vec))
    
    for (i_drop_time in seq_along(drop_times)) {
        if (i_drop_time == 1) {
            surv_sample[time_vec < drop_times[i_drop_time]] <- 1
        }
        
        if (i_drop_time < length(drop_times)) {
            surv_sample[time_vec >= drop_times[i_drop_time] &
                            time_vec < drop_times[i_drop_time+1]] <- surv_sample_values[i_drop_time]
        } else {
            surv_sample[time_vec >= drop_times[i_drop_time]] <- surv_sample_values[i_drop_time]
        }
    }
    
    return(surv_sample)
}


calculate_interval_censored_CIs_from_single_imputation <- function(time_vec,
                                                                   single_imp_out,
                                                                   alpha) {
    vars_imp <- apply(single_imp_out$all_surv_funcs, 2, var)
    
    # For adapted Greenwood formula:
    # variance due to imputation + estimated variance due to sampling
    total_variance_est <- vars_imp + colMeans(single_imp_out$all_sd^2)
    
    # Construct adapted Greenwood's confidence interval:
    ci_width <- qnorm(alpha/2, lower.tail = FALSE) * sqrt(total_variance_est)
    S_hat <- colMeans(single_imp_out$all_surv_funcs)
    ci_lower <- pmax(0, S_hat - ci_width)
    ci_upper <- pmin(1, S_hat + ci_width)
    
    # Fill in missing CI values with their last value or zero
    for (i_time in seq_along(time_vec)) {
        if (is.na(ci_lower[i_time])) {
            ci_lower[i_time] <- 0
        }
        if (is.na(ci_upper[i_time])) {
            ci_upper[i_time] <- ci_upper[i_time-1]
        }
    }
    
    # Upper limit should never increase, since every member of the ensemble
    # must be a valid survival curve
    for (i_time in 2:length(time_vec)) {
        ci_upper[i_time] <- min(ci_upper[i_time], ci_upper[i_time-1])
    }
    
    return(data.frame(S_hat = S_hat,
                      ci_lower = ci_lower,
                      ci_upper = ci_upper,
                      time = time_vec))
}


run_grouped_single_imputation <- function(sample, n_mc_runs, time_vec,
                                          grouping_var, alpha,
                                          use_exponential_greenwood_CIs = FALSE) {
    level_names <- levels(sample[[grouping_var]])
    
    df_CIs_all <- foreach(level_name = level_names,
                          .packages = "survival",
                          .combine = rbind) %do% {
                              single_imp_out <- run_single_imputation(
                                  n_mc_runs,
                                  time_vec,
                                  sample[sample[[grouping_var]] == level_name, ],
                                  use_exponential_greenwood_CIs = use_exponential_greenwood_CIs
                              )
                              
                              df_CIs <- calculate_interval_censored_CIs_from_single_imputation(time_vec,
                                                                                               single_imp_out,
                                                                                               alpha)
                              df_CIs[[grouping_var]] <- level_name
                              df_CIs
                          }
    
    return(df_CIs_all)
}


plot_single_imputation_results <- function(df_CIs, grouping_var, alpha = 0.4) {
    g <- ggplot(df_CIs, aes(x = time)) +
        geom_line(aes(y = S_hat, colour = .data[[grouping_var]]), size = 1) +
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = .data[[grouping_var]]),
                    alpha = alpha) +
        theme_classic(base_size = 25) +
        labs(x = "Time (days)", y = "Survival probability")
    
    return(g)
}


prepare_survival_data <- function(therapy_data, last_visits) {
    therapy_data %<>%
        group_by(subjid) %>%
        mutate(cmendy = ifelse(is.na(cmendy),
                               last_visits$last_visit[last_visits$subjid == first(subjid)],
                               cmendy),
               surv_time = cmendy - cmstdy,
               event = 1 - cmenrf) %>%
        ungroup()
    
    return(therapy_data)
}


collect_surv_outputs <- function(surv_obj) {
    # Helper function for summarising survival fits from within
    # a dplyr pipeline
    surv_fit <- survfit(surv_obj ~ 1)
    tibble(time = surv_fit$time, surv = surv_fit$surv,
           lower = surv_fit$lower, upper = surv_fit$upper)
}