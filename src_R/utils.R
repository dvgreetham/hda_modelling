library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(msm)
library(stringr)
library(caret)
library(tibble)


read_data <- function(file, drop_na = TRUE) {
    if (tools::file_ext(file) == "feather") {
        # TODO: fallback to CSV if arrow not installed or feather file doesn't exist?
        requireNamespace("arrow", quietly = TRUE)
        df <- arrow::read_feather(file)
    } else {
        df <- read.csv(file, sep = "\t", header = TRUE) %>%
            prepare_factors()
    }
    
    # Select only subjects that can be used for survival analysis
    if ("fluct" %in% colnames(df)) {
        df <- df %>%
            dplyr::filter(fluct == 0) %>% # subjects without fluctuations
            select(-c("fluct")) # drop constant column
    }
    if ("hdcat_0" %in% colnames(df)) {
        df <- df %>%
            dplyr::filter(hdcat_0 == 3) %>% # subjects with manifest HD at enrollment
            select(-c("hdcat_0")) # drop constant column
    }
    if (drop_na) {
        df <- df %>%
            tidyr::drop_na(any_of(c("duration_min", "duration_max")))
    }
    
    df
}


create_first_order_forward_chain_intensity_matrix <- function(n_states) {
    Q <- matrix(0, nrow = n_states, ncol = n_states)
    diag(Q[, -1]) <- rep.int(0.1, n_states-1)
    diag(Q) <- c(rep.int(-0.1, n_states-1), 0)
    return(Q)
}


generate_pca_plots <- function(pca_results) {
    # Function for producing plots of the variance of the PCs
    # and the proportion of variance explained.
    
    pc_plot_data <- tibble(pc = seq_along(pca_results$sdev),
                           variances = pca_results$sdev^2,
                           cumulative_proportion = cumsum(variances / sum(variances)))
    
    var_plot <- ggplot(pc_plot_data, aes(x = pc, y = variances)) +
        geom_col() +
        ggtitle("Variance of principal components") +
        scale_x_continuous(breaks = seq_along(pca_results$sdev))
    
    var_explained_plot <- ggplot(pc_plot_data, aes(x = pc, y = cumulative_proportion)) +
        geom_point() +
        geom_line() +
        ggtitle("Proportion of variance explained") +
        scale_x_continuous(breaks = seq_along(pca_results$sdev))
    
    return(list(var_plot = var_plot, var_explained_plot = var_explained_plot))
}


plot_principal_component_variability <- function(observations_data,
                                                 grouping_var,
                                                 n_pcs_to_plot = NULL) {
    # Produce a plot of the variability of the PCA loadings
    # of the data in observations_data between groups defined by grouping_var.
    if (is.null(n_pcs_to_plot)) {
        n_pcs_to_plot <- ncol(observations_data$PCs)
    }
    
    if (n_pcs_to_plot > 9) {
        pc_names <- paste0("pc_", c(paste0("0", 1:9), as.character(10:n_pcs_to_plot)))
    } else {
        pc_names <- paste0("pc_", seq_len(n_pcs_to_plot))
    }
    
    for (i_pc in seq_len(n_pcs_to_plot)) {
        observations_data %<>%
            mutate(!!pc_names[i_pc] := PCs[, i_pc])
    }
    
    observations_data %<>%
        select(!!grouping_var, emission, all_of(pc_names)) %>%
        distinct() %>%
        pivot_longer(all_of(pc_names), names_to = "PC", values_to = "loading")
    
    sign_flip_data <- find_PCA_sign_flips(observations_data, grouping_var)
    observations_data %<>%
        apply_sign_flips_to_PC_data(sign_flip_data, grouping_var)
    
    loadings_plot <- ggplot(observations_data, aes(x = as.numeric(.data[[grouping_var]]),
                                                   y = loading,
                                                   colour = emission)) +
        geom_point() +
        geom_line() +
        facet_wrap(vars(PC)) +
        scale_x_continuous(breaks = seq_len(nlevels(observations_data[[grouping_var]])),
                           labels = levels(observations_data[[grouping_var]])) +
        labs(x = grouping_var) +
        theme(text = element_text(size = 15),
              panel.spacing = unit(2, "lines"))
    
    return(loadings_plot)
}


apply_sign_flips_to_PC_data <- function(pc_data, sign_flip_data, grouping_var) {
    # Using the sign flip data (see find_PCA_sign_flips),
    # apply the appropriate sign flips to the loadings in pc_data
    # (the output of do_grouped_PCA),
    # where we assume the PCs have been calculated after grouping by grouping_var.
    
    pc_data %<>%
        left_join(sign_flip_data, by = c(grouping_var, "PC")) %>%
        mutate(loading = ifelse(apply_sign_flip, -loading, loading)) %>%
        select(-apply_sign_flip)
    
    return(pc_data)
}


find_PCA_sign_flips <- function(observations_data, grouping_var) {
    # The signs of the PCA loadings may not be consistent between groups
    # even if the PCs are very similar -- principal components are not unique,
    # and small numerical differences can lead to a flip in the direction.
    # The PCs may flip sign between groups, which leads to confusing plots.
    # This function attempts to clarify comparisons of loadings between groups
    # by making the signs of the loadings more consistent.
    
    # The idea is for each PC, take the emission whose loading for that PC
    # is the largest in absolute value on average across groups,
    # and make the sign of that loading consistent across groups.
    # This function returns a data frame indicating which PCs should have
    # their signs flipped.
    biggest_loadings <- observations_data %>%
        group_by(PC, emission) %>%
        summarise(mean_abs = mean(abs(loading)), .groups = "drop_last") %>%
        summarise(biggest_loading = emission[which.max(mean_abs)],
                  .groups = "drop_last")
    
    sign_data <- observations_data %>%
        group_by(PC) %>%
        filter(.data[[grouping_var]] == first(.data[[grouping_var]])) %>%
        summarise(sign_of_biggest_loading =
                      sign(loading[emission ==
                                       biggest_loadings$biggest_loading[biggest_loadings$PC ==
                                                                            first(PC)]]),
                  .groups = "drop_last")
    
    biggest_loadings %<>% left_join(sign_data, by = "PC")
    
    sign_flip_data <- observations_data %>%
        group_by(PC) %>%
        filter(emission ==
                   biggest_loadings$biggest_loading[biggest_loadings$PC ==
                                                        first(PC)]) %>%
        group_by(.data[[grouping_var]]) %>%
        mutate(apply_sign_flip =
                   sign(loading) !=
                   biggest_loadings$sign_of_biggest_loading[biggest_loadings$PC == first(PC)]) %>%
        select(!!grouping_var, PC, apply_sign_flip)
    
    return(sign_flip_data)
}


do_grouped_PCA <- function(obs_data, grouping_var, preproc,
                           n_pcs_to_keep = NULL) {
    # PCA is performed in groups, where obs_data is grouped by grouping_var.
    # A data frame containing various summary information,
    # including the loadings ("PCs"), is returned in pc_data,
    # and the transformed data (scores) are returned in transformed_data,
    # The observations in obs_data are preprocessed prior to grouping
    # using the caret::preProcess object preproc,
    # and the PCA is subsequently performed in groups without further
    # centering or scaling.
    # n_pcs_to_keep is the number of principal components to keep.
    
    emission_names <- colnames(obs_data$obs)
    
    obs_data$obs_z <- predict(preproc, obs_data$obs)
    
    obs_data_grouped <- obs_data %>%
        select(subjid, age, dssage, !!grouping_var, obs, obs_z) %>%
        group_by(.data[[grouping_var]])
    
    transformed_data <- obs_data_grouped %>%
        mutate(obs_transformed_grouped =
                   extract_PC_outputs(obs_z, n_pcs_to_keep)$transformed_data) %>%
        ungroup()
    
    pc_data <- obs_data_grouped %>%
        summarise(means = colMeans(obs_z),
                  covariance = cov(obs_z),
                  extract_PC_outputs(obs_z, n_pcs_to_keep)$loadings_data,
                  .groups = "keep") %>%
        mutate(emission = emission_names,
               emission_num = as.factor(seq_along(emission_names))) %>%
        ungroup()
    
    return(list(transformed_data = transformed_data,
                pc_data = pc_data))
}


extract_PC_outputs <- function(obs_z, n_pcs_to_keep = NULL) {
    # Utility function for extracting the outputs from prcomp
    # to use in dplyr::summarise() or dplyr::mutate()
    
    prcomp_out <- prcomp(obs_z, center = FALSE, scale. = FALSE,
                         rank. = n_pcs_to_keep)
    
    loadings_data <- tibble(PCs = prcomp_out$rotation)
    
    if (length(prcomp_out$sdev) == dim(prcomp_out$rotation)[1]) {
        loadings_data$PC_vars <- prcomp_out$sdev^2
    } else {
        loadings_data$PC_vars <- NA
    }
    
    return(list(transformed_data = prcomp_out$x,
                loadings_data = loadings_data))
    
}


construct_survival_plot_data <- function(msm_fit, time_vec, state_vec, label) {
    # Function for constructing a tibble used for plotting survival curves
    # based on an msm fit.
    # Expects first state is called state 0 in state_vec
    sojourn_times <- sojourn.msm(msm_fit)
    
    plot_data <- tibble(time = time_vec,
                        state = state_vec,
                        survival_prob = exp(-time / sojourn_times$estimates[state+1]),
                        curve_type = label)
    
    if (is.null(sojourn_times$L)) {
        plot_data %<>%
            mutate(lower = survival_prob,
                   upper = survival_prob)
    } else {
        plot_data %<>%
            mutate(lower = exp(-time / sojourn_times$L[state+1]),
                   upper = exp(-time / sojourn_times$U[state+1]))
    }
    return(plot_data)
}


plot_medication_costs <- function(pharmacotx_data, grouping_var, n_subjs) {
    
    plot_data <- pharmacotx_data %>%
        group_by(.data[[grouping_var]]) %>%
        summarise(cost_per_patient_avg = sum(total_cost) / n_subjs,
                  ingredient = first(cmtrt__ing),
                  .groups = "drop_last") %>%
        arrange(cost_per_patient_avg) %>%
        tail(n = 15)
    
    grouping_name <- switch(grouping_var,
                            cmtrt__ing = "active ingredient",
                            cmtrt__modify = "product",
                            cmindc__modify = "indication")
    
    labels <- switch(grouping_var,
                     cmtrt__modify = str_trunc(paste0(plot_data[[grouping_var]], " (",
                                                      plot_data$ingredient, ")"), 30),
                     cmtrt__ingredient = str_trunc(plot_data$ingredient),
                     str_trunc(plot_data[[grouping_var]], 30))
    
    cost_plot <- ggplot(plot_data, aes(x = .data[[grouping_var]],
                                       y = cost_per_patient_avg)) +
        geom_col() +
        coord_flip() +
        scale_x_discrete(limits = plot_data[[grouping_var]], labels = labels) +
        theme(text = element_text(size = 15)) +
        labs(y = "Cost incurred so far by \n the average patient in Enroll (GBP)",
             x = grouping_name,
             title = paste0("Cost of medication by ", grouping_name))
    
    return(cost_plot)
}


extract_quantities <- function(string_vec) {
    # Given an input vector of strings with drug names and a previously
    # extracted mass, extract the measurement units
    
    quantity <- str_extract_all(string_vec, "\\d*,?\\d+\\.?\\d*([:alpha:]|/)+")
    
    amount <- lapply(quantity,
                     function(s) unlist(str_extract_all(s, "\\d*,?\\d+\\.?\\d*")))
    
    unit <- lapply(seq_along(quantity), function(idx) {
        if (is.null(amount[[idx]])) {
            return(NULL)
        }
        
        return(str_remove_all(quantity[[idx]], amount[[idx]]))
    })
    
    unit <- lapply(unit, function(s) {
        if (is.null(s)) {
            return(NULL)
        }
        
        out <- str_replace_all(s, c("micrograms" = "microgram",
                                    "microg(?!ram)" = "microgram",
                                    "mcg" = "microgram",
                                    "unit" = "IU"))
        out <- str_replace_all(out, c("microgram/dose" = "microgram"))
        out <- str_remove_all(out, "/$")
        return(out)
    })
    
    amount <- lapply(amount, function(x) {
        if (is.null(x)) {
            return(NULL)
        }
        
        out <- as.numeric(str_remove_all(x, ","))
        
        return(out)
    })
    
    # Handle ethinylestradiol
    idx_ethinylestradiol <- which(str_detect(str_to_lower(string_vec),
                                             "ethinylestradiol"))
    for (idx in idx_ethinylestradiol) {
        amount[[idx]] <- amount[[idx]][1]
        unit[[idx]] <- unit[[idx]][1]
    }
    
    
    return(tibble(amount = amount, unit = unit))
}


use_last_visit_as_end_date <- function(therapy_data, last_visits) {
    # Update the therapy data to impute missing therapy end dates
    # with the date of the patient's last visit
    therapy_data %<>%
        group_by(subjid) %>%
        mutate(cmendy = ifelse(
            is.na(cmendy),
            last_visits$last_visit[last_visits$subjid == first(subjid)],
            cmendy
        )) %>%
        ungroup()
    
    return(therapy_data)
}


add_observations_of_state_1_at_age_0 <- function(visit_data, obs_var = "obs") {
    # Set up data to have extra observation of state 1 at age 0
    n_obs <- ncol(visit_data[[obs_var]])
    
    if ("obstrue" %in% colnames(visit_data)) {
        visit_data_extra_obs <- visit_data %>%
            select(subjid, age, obstrue, !!obs_var)
    } else {
        visit_data_extra_obs <- visit_data %>%
            select(subjid, age, !!obs_var) %>%
            mutate(obstrue = NA)
    }
    
    if ("obstype" %in% colnames(visit_data)){
        visit_data_extra_obs$obstype <- visit_data$obstype
    } else {
        visit_data_extra_obs$obstype <- 1
    }
    
    extra_rows <- tibble(subjid = unique(visit_data$subjid), age = 0) %>%
        mutate(obstrue = 1, obstype = 2)
    
    if (is.null(n_obs)) {
        extra_rows[[obs_var]] <- NA
    } else {
        extra_rows[[obs_var]] <- matrix(rep.int(NA, n_obs * nrow(extra_rows)),
                                        nrow = nrow(extra_rows),
                                        ncol = n_obs, byrow = FALSE)
    }
    
    visit_data_extra_obs %<>%
        rbind(extra_rows) %>%
        arrange(subjid, age)
    
    return(visit_data_extra_obs)
}


add_death_observations <- function(visit_data, label_death, idx_deathstate,
                                   age_death_max, obs_var = "obs") {
    # For patients with an observed age of death,
    # add this as an extra row in the data
    n_obs <- ncol(visit_data[[obs_var]])
    
    if ("obstrue" %in% colnames(visit_data)) {
        visit_data_extra_obs <- visit_data %>%
            select(subjid, age, obstrue, !!obs_var)
    } else {
        visit_data_extra_obs <- visit_data %>%
            select(subjid, age, !!obs_var) %>%
            mutate(obstrue = NA)
    }
    
    visit_data_extra_obs %<>% mutate(obstype = 1)
    
    # Add 0.5 years on to age of death to reduce maximum error to 6 months
    extra_rows <- visit_data %>%
        group_by(subjid) %>%
        summarise(has_death_date = any(!is.na(dssage)),
                  age = first(dssage) + 0.5,
                  .groups = "drop_last") %>%
        filter(has_death_date) %>%
        select(-has_death_date) %>%
        mutate(obstrue = NA,
               obstype = 3)
    
    extra_rows[[obs_var]] <- matrix(rep.int(label_death, n_obs * nrow(extra_rows)),
                                    nrow = nrow(extra_rows),
                                    ncol = n_obs, byrow = FALSE)
    
    for (subj in extra_rows$subjid) {
        visit_ages_this_subj <- visit_data %>%
            filter(subjid == subj) %>%
            pull(age)
        
        if (any(visit_ages_this_subj >=
                extra_rows %>% filter(subjid == subj) %>% pull(age))) {
            # If age of last visit is greater than or equal to age of death,
            # assume patient died 3 months after last visit
            extra_rows[extra_rows$subjid == subj, "age"] <- 
                max(visit_ages_this_subj) + 0.25
        }
    }

    visit_data_extra_obs %<>%
        rbind(extra_rows) %>%
        arrange(subjid, age)
    
    # Add artificial interval-censored death observation at age_death_max
    # for patients with no true death observation
    extra_rows <- visit_data %>%
        group_by(subjid) %>%
        summarise(has_death_date = any(!is.na(dssage)),
                  max_visit_age = max(age),
                  .groups = "drop_last") %>%
        filter(!has_death_date) %>%
        select(-has_death_date) %>%
        mutate(age = pmax(age_death_max, max_visit_age + 1),
               obstrue = NA,
               obstype = 1) %>%
        select(-max_visit_age)
    
    extra_rows[[obs_var]] <- matrix(rep.int(label_death, n_obs * nrow(extra_rows)),
                                    nrow = nrow(extra_rows),
                                    ncol = n_obs, byrow = FALSE)
    
    visit_data_extra_obs %<>%
        rbind(extra_rows) %>%
        arrange(subjid, age)
    
    return(visit_data_extra_obs)
    
}


calculate_training_log_likelihood <- function(msm_fits, n_states_vec) {
    # Calculate the log-likelihood of the msm models given in the list msm_fits
    # on the training set. n_states_vec is the vector of the number of states
    # used in the models in msm_fits.
    training_llhd <- numeric(length(n_states_vec))
    for (i_fit in seq_along(n_states_vec)) {
        training_llhd[i_fit] <- logLik(msm_fits[[i_fit]])
    }
    plot_data <- data.frame(n_states = n_states_vec,
                            training_llhd = training_llhd)
    
    return(plot_data)
}


calculate_test_log_likelihood <- function(msm_fits, test_data, n_states_vec,
                                          plot_data = NULL, obs_var = "obs",
                                          label_death) {
    # Calculates the log-likelihood of the msm models given in the list msm_fits
    # on the test data set test_data.
    # n_states_vec is the vector of the number of states in each of the models;
    # plot_data is an optional plot data.frame
    # such as the one produced by calculate_training_log_likelihood,
    # and obs_var is a string giving the variable name to take as the
    # observations (emissions).
    if (is.null(plot_data)) {
        plot_data <- data.frame(n_states = n_states_vec)
    }
    plot_data$test_llhd <- numeric(length(n_states_vec))
    for (i_fit in seq_along(n_states_vec)) {
        Q <- qmatrix.msm(msm_fits[[i_fit]])$estimate
        hmm_model <- set_up_hmodel_from_msm_fit(msm_fits[[i_fit]], label_death)
        .GlobalEnv$test_data_fit <- test_data[[i_fit]]
        msm_fit_test <- msm(as.formula(paste0(obs_var, " ~ age")),
                            subject = subjid,
                            data = test_data_fit,
                            qmatrix = Q, hmodel = hmm_model,
                            obstrue = test_data_fit$obstrue,
                            deathexact = n_states_vec[i_fit],
                            est.initprobs = FALSE,
                            method = "BFGS",
                            fixedpars = TRUE,
                            control = list(maxit = 1, fnscale = 5e4))
        
        plot_data$test_llhd[i_fit] <- logLik(msm_fit_test)
    }
    
    return(plot_data)
}


make_PCA_signs_more_consistent <- function(pc_data, obs_data, n_pcs_to_keep,
                                           grouping_var = "age_bracket") {
    # Given data.frames pc_data containing information on the PC loadings and
    # obs_data (these two can be produced by do_grouped_PCA),
    # this function attempts to make the directions of the principal components
    # more consistent between groups according to the grouping variable grouping_var.
    # n_pcs_to_keep is the number of leading principal components to keep.
    for (i_pc in seq_len(n_pcs_to_keep)) {
        pc_data %<>%
            mutate(!!glue("pc_{i_pc}") := PCs[, i_pc])
        
        obs_data %<>%
            mutate(!!glue("pc_{i_pc}") := obs_transformed_grouped[, i_pc])
    }
    
    pc_data %<>%
        select(!!grouping_var, emission, glue("pc_{seq_len(n_pcs_to_keep)}")) %>%
        distinct() %>%
        pivot_longer(glue("pc_{seq_len(n_pcs_to_keep)}"), names_to = "PC",
                     values_to = "loading")
    
    # Signs of PC scores may be inconsistent between age brackets,
    # so make them more consistent for easier comparison
    sign_flip_data <- find_PCA_sign_flips(pc_data, grouping_var)
    
    obs_data %<>%
        ungroup() %>%
        pivot_longer(glue("pc_{seq_len(n_pcs_to_keep)}"), names_to = "PC",
                     values_to = "score") %>%
        left_join(sign_flip_data, by = c(grouping_var, "PC")) %>%
        mutate(score = ifelse(apply_sign_flip, -score, score)) %>%
        select(-apply_sign_flip) %>%
        pivot_wider(names_from = PC, values_from = score)
    
    obs_data %<>%
        mutate(obs_transformed_grouped =
                   as.matrix(obs_data[, glue("pc_{seq_len(n_pcs_to_keep)}")]))
    
    pc_data %<>% apply_sign_flips_to_PC_data(sign_flip_data, grouping_var)
    
    return(list(pc_data = pc_data, obs_data = obs_data))
}


calculate_AIC_BIC <- function(plot_data, msm_fits, n_obs) {
    # Calculates the Akaike Information Criterion (AIC) and Bayesian Information
    # Criterion (BIC) for the models given in msm_fits.
    # plot_data is a data.frame with at least a variable n_states giving the
    # number of states for each model.
    # n_obs is the number of observations (i.e. records, not emissions) in the
    # training data.
    plot_data_AIC <- plot_data %>%
        select(n_states) %>%
        mutate(criterion = sapply(msm_fits, AIC),
               type = "AIC")
    
    # N.b. the R function BIC doesn't work,
    # as there is no nobs method implemented for the msm class.
    plot_data_BIC <- plot_data_AIC %>%
        mutate(
            criterion = sapply(msm_fits, function(msm_fit) {
                -2*logLik(msm_fit) +
                    msm_fit$paramdata$npars * log(n_obs)
            }),
            type = "BIC"
        )
    
    plot_data_AIC_BIC <- bind_rows(plot_data_AIC, plot_data_BIC)
    
    return(plot_data_AIC_BIC)
}


do_grouped_PCA_list <- function(obs_data, preproc, n_states_vec, n_pcs_to_keep) {
    # Given a list of numbers of states n_states_vec and caret::preProcess
    # object preproc, this function applies grouped PCA on the data given in
    # obs_data, where the grouping is into the number of age brackets equal to
    # the number of states being fitted.
    # 
    # We then correct sign differences
    # between groups as much as possible, and return the results in two lists:
    # pc_data, which is the corresponding loadings data, and
    # transformed_data, which contains the scores after keeping the first
    # n_pcs_to_keep principal components.
    transformed_data <- list()
    pc_data <- list()
    
    for (i_model in seq_along(n_states_vec)) {
        grouping_var <- glue("age_bracket_{n_states_vec[i_model]}")
        this_train <- obs_data %>%
            mutate(!!grouping_var := cut(age, n_states_vec[i_model] - 1))
        
        pc_results <- this_train %>%
            do_grouped_PCA(grouping_var, preproc,
                           n_pcs_to_keep = n_pcs_to_keep)
        
        pc_data[[i_model]] <- pc_results$pc_data
        
        transformed_data[[i_model]] <- pc_results$transformed_data
        
        out <- make_PCA_signs_more_consistent(pc_data[[i_model]],
                                              transformed_data[[i_model]],
                                              n_pcs_to_keep,
                                              grouping_var)
        
        pc_data[[i_model]] <- out$pc_data
        transformed_data[[i_model]] <- out$obs_data
    }
    
    return(list(pc_data = pc_data, transformed_data = transformed_data))
}


transform_raw_observations_to_grouped_PCs <- function(obs_data, pc_data, preproc,
                                                      n_states_vec, n_pcs) {
    # This function is for transforming a set of data obs_data into the
    # principal component basis given by the loadings in pc_data.
    # N.b. the obs_data does not need to be the data from which the PCs were
    # calculated, unlike in do_grouped_PCA.
    # This is used for applying a consistent PC tranformation to a validation
    # or test set.
    # preproc is the caret::preProcess object used to center and scale the data,
    # n_states_vec is the vector of numbers of states corresponding to the
    # sequence of models we wish to fit,
    # and n_pcs is the number of principal components we want to keep.
    transformed_data <- list()
    
    for (i_model in seq_along(n_states_vec)) {
        pc_data[[i_model]] %<>%
            pivot_wider(names_from = PC, values_from = loading)
        
        grouping_var <- glue("age_bracket_{n_states_vec[i_model]}")
        
        transformed_data[[i_model]] <- obs_data %>%
            mutate(obs_z = predict(preproc, obs),
                   !!grouping_var := cut(age, n_states_vec[i_model] - 1)) %>%
            group_by(.data[[grouping_var]]) %>%
            mutate(
                obs_transformed_grouped =
                    obs_z %*%
                    as.matrix(
                        # as.numeric here because age boundaries of test data set
                        # will not exactly match those of training set
                        pc_data[[i_model]][as.numeric(pc_data[[i_model]][[grouping_var]]) ==
                                               as.numeric(first(.data[[grouping_var]])),
                                           glue("pc_{seq_len(n_pcs)}")]
                    )
            ) %>%
            ungroup()
        
        colnames(transformed_data[[i_model]]$obs_transformed_grouped) <-
            paste0("pc_", seq_len(n_pcs))
    }
    
    return(transformed_data)
}


plot_state_transition_matrix <- function(msm_fits, k_folds_xval = 1, 
                                         i_model = NULL, n_years = 1) {
    plot_data <- data.frame()
    
    for (i_fold in seq_len(k_folds_xval)) {
        if (k_folds_xval > 1) {
            P_est <- pmatrix.msm(msm_fits[[i_fold]][[i_model]], t = n_years)
        } else {
            P_est <- pmatrix.msm(msm_fits, t = n_years)
        }
        class(P_est) <- "matrix"
        probs_data <- as_tibble(P_est)
        colnames(probs_data) <- as.character(seq_len(ncol(probs_data)))
        
        probs_data %<>%
            rownames_to_column("from_state") %>%
            pivot_longer(!from_state, names_to = "to_state", values_to = "prob") %>%
            mutate(fold = i_fold)
        
        plot_data %<>%
            bind_rows(probs_data)
    }
    
    prob_plot <- ggplot(plot_data, aes(x = to_state, y = from_state)) +
        geom_tile(aes(fill = prob)) +
        geom_text(aes(label = round(prob, 2))) +
        scale_y_discrete(limits = rev) +
        scale_fill_gradient(low = "white", high = "blue")
    
    if (k_folds_xval > 1) {
        prob_plot <- prob_plot + facet_wrap(vars(fold))
    }
    
    return(prob_plot)
}


get_hd_life_expectancy <- function() {
    # Using data from NPD/10724/TD/TL/009,
    # estimate the mean age of death of HD patients
    
    # Use the NCDR figure, as this had more patients.
    # Not using an average of the two data sources,
    # because there will be significant overlap in included patients.
    mean_death_age <- 63.9
    
    return(mean_death_age)
}


get_death_ages_and_last_visits <- function(enroll_train) {
    # T_i is the vector of observed death ages for patients that have them,
    # T_j is the vector of ages at last visit for patients without observed
    # death ages.
    
    subjs_with_death_obs <- enroll_train %>%
        group_by(subjid) %>%
        summarise(has_death_date = any(!is.na(dssage)), .groups = "drop_last") %>%
        filter(has_death_date) %>%
        pull(subjid)
    
    T_i <- enroll_train %>%
        filter(subjid %in% subjs_with_death_obs) %>%
        group_by(subjid) %>%
        summarise(dssage = first(dssage), .groups = "drop_last") %>%
        pull(dssage)
    
    T_j <- enroll_train %>%
        filter(!(subjid %in% subjs_with_death_obs)) %>%
        group_by(subjid) %>%
        summarise(age_last_visit = max(age), .groups = "drop_last") %>%
        pull(age_last_visit)
    
    return(list(T_i = T_i, T_j = T_j))
}


calculate_artificial_max_death_age <- function(target_life_expectancy, T_i, T_j) {
    # See analysis/5_model_fitting/3_impact_of_artificial_death_observations
    # for explanation and derivation
    
    m <- length(T_i)
    n <- length(T_j)
    
    q <- 1 / target_life_expectancy
    
    A <- 1/12 * q^2 * n
    B <- -1/2 * q * (1/3 * q * sum(T_j) + n)
    C <- 1/12 * q^2 * sum(T_j^2) - q * (sum(T_i) + 1/2 * sum(T_j)) + m + n
    
    T_max <- (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    
    return(T_max)
}