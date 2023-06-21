# Add parallel optimisation 

msm.optim.optimParallel <- function(p, gr, hessian, msmdata, qmodel, qcmodel, cmodel, hmodel, ...) {
    optim.args <- list(...)
    
    if (is.null(optim.args$control)) optim.args$control <- list()
    
    optim.args <- c(optim.args, list(par=p$inits, fn=msm:::lik.msm, hessian=hessian, gr=gr,
                                     msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                     cmodel=cmodel, hmodel=hmodel, paramdata=p))
    opt <- do.call("optimParallel", optim.args)
    if (opt$convergence==1)
        warning("Iteration limit in optim() reached without convergence. Reported estimates are not the maximum likelihood. Increase \"maxit\" or change optimisation method - see help(optim) and help(msm).")
    else if (opt$convergence==10)
        warning("Not converged: Nelder-Mead simplex is degenerate. Reported estimates are not the maximum likelihood.")
    else if (opt$convergence %in% 51:52)
        warning("Not converged: error in L-BFGS-B, see help(optim). Reported estimates are not the maximum likelihood.")
    ctrace <- !is.null(optim.args$control$trace) && optim.args$control$trace > 0
    if (ctrace){
        cat("Used", opt$counts[1], "function and", opt$counts[2], "gradient evaluations\n")
    }
    if (!is.null(opt$message)) warning("optimParallel() returned a message: ",opt$message)
    p$lik <- opt$value
    p$params[p$optpars] <- opt$par
    p$opt <- opt
    p
}


calculate_initial_parameter_values_for_msm <- function(df_extra_obs,
                                                       n_states,
                                                       label_death,
                                                       obs_var = "obs") {
    # function for calculating some rough initial parameter values for the
    # msm fit based on the data.
    # The input data frame df_extra_obs should have had an extra observation
    # added to it at age 0 for each patient.
    # n_states is the number of states in the model, including death
    # n_obs is the number of observed variables.
    
    n_hidden_states = n_states - 1  # Exclude death
    
    n_obs <- ncol(df_extra_obs[[obs_var]])
    
    age_summary <- df_extra_obs %>%
        filter(age > 0) %>%
        summarise(min_age = min(age), max_age = max(age))
    
    age_step <- (age_summary$max_age - age_summary$min_age) / n_hidden_states
    age_lower <- seq(age_summary$min_age,
                     age_summary$min_age + (n_hidden_states-1) * age_step,
                     length.out = n_hidden_states)
    age_upper <- age_lower + age_step
    
    init_means <- matrix(0, nrow = n_hidden_states, ncol = n_obs)
    init_sds <- matrix(0, nrow = n_hidden_states, ncol = n_obs)
    
    for (i_state in seq_len(n_hidden_states)) {
        data_in_this_age_bracket <- df_extra_obs %>%
            filter(df_extra_obs[[obs_var]][, 1] != label_death,
                   age > age_lower[i_state],
                   age <= age_upper[i_state])
        for (j_obs in seq_len(n_obs)) {
            init_means[i_state, j_obs] <- mean(data_in_this_age_bracket[[obs_var]][, j_obs],
                                               na.rm = TRUE)
            init_sds[i_state, j_obs] <- sd(data_in_this_age_bracket[[obs_var]][, j_obs],
                                           na.rm = TRUE)
        }
        
    }

    colnames(init_means) <- colnames(df_extra_obs[[obs_var]])
    colnames(init_sds) <- colnames(df_extra_obs[[obs_var]])
    
    sojourn_times <- c(age_summary$min_age + age_step,
                      rep.int(age_step, n_states - 2))
    
    # Scale down initial sojourn times based on expected lifespan
    expected_lifespan <- 60  # years
    sojourn_times <- sojourn_times * expected_lifespan / sum(sojourn_times)
    
    return(list(init_means = init_means, init_sds = init_sds,
                sojourn_times = sojourn_times))
}


set_up_hmodel_from_msm_fit <- function(msm_fit, label_death) {
    # Function for setting up a hmodel that can be input to msm
    # from the output of msm.
    
    n_states <- msm_fit$hmodel$nstates
    n_obs <- msm_fit$hmodel$nout[1]
    init_emis_params <- list(init_means = matrix(0, nrow = n_states - 1, ncol = n_obs),
                             init_sds = matrix(0, nrow = n_states - 1, ncol = n_obs))
    
    for (i_state in seq_len(n_states-1)) {
        for (i_obs in seq_len(n_obs)) {
            init_emis_params$init_means[i_state, i_obs] <-
                msm_fit$hmodel$pars[2*(i_state-1)*n_obs + 2*(i_obs-1) + 1]
            init_emis_params$init_sds[i_state, i_obs] <-
                msm_fit$hmodel$pars[2*(i_state-1)*n_obs + 2*(i_obs-1) + 2]
        }
    }
    
    hmodel <- set_up_hmm_model(init_emis_params, seq_len(n_obs),
                               label_death)
    
    return(hmodel)
}


set_up_hmm_model <- function(init_emis_params, obs_to_select, label_death) {
    # Set initial hmodel parameters -- the hidden Markov model fit
    # requires initialisation of the "hmodel" argument.
    #
    # Arguments:
    # init_emis_params -- initial emission parameters from
    #                       calculate_initial_parameter_values_for_msm
    # obs_to_select -- indices of observations to use
    # label_death -- a number to use as a death label
    
    n_states <- nrow(init_emis_params$init_means) + 1
    hmodel <- list()
    for (i_state in seq_len(n_states-1)) {
        hmmMV_args <- list()
        for (i_obs in seq_along(obs_to_select)) {
            .GlobalEnv$this_mean <- init_emis_params$init_means[i_state, obs_to_select[i_obs]]
            .GlobalEnv$this_sd <- init_emis_params$init_sds[i_state, obs_to_select[i_obs]]
            hmmMV_args[[i_obs]] <- hmmNorm(this_mean, this_sd)
        }
        hmodel[[i_state]] <- do.call(hmmMV, hmmMV_args)
    }
    
    hmodel[[n_states]] <- hmmIdent(label_death)
    
    return(hmodel)
}


create_FC_intensity_matrix_from_sojourn_times <- function(sojourn_times, order = 1) {
    # Given a set of sojourn times (expected times in each state),
    # create a forward-chain intensity matrix (i.e. transitions to earlier
    # states are not allowed) for a Markov process.
    # The order indicates how many states forward one can instantaneously go,
    # i.e. order = 1 means you can only transition to the next state,
    # while order = 2 means you can jump two states, e.g. go from state 1 to 3.
    # This is used for intialising the continuous-time hidden Markov model fit.
    
    # Last state is absorbing, so assume sojourn time for that is not provided
    n_states <- length(sojourn_times) + 1
    Q <- matrix(0, nrow = n_states, ncol = n_states)
    if (n_states > 2) {
        diag(Q[, -1]) <- 1 / sojourn_times
    } else if (n_states == 2) {
        Q[1, 2] <- 1 / sojourn_times
    } else {
        stop("Fewer than two states is not supported")
    }
    
    if (order == 1) {
        diag(Q) <- c(-1 / sojourn_times, 0)
    } else if (order == 2) {
        # Assume half intensity for skipping a state
        if (n_states > 3) {
            diag(Q[, -c(1, 2)]) <- 1 / sojourn_times[seq_len(n_states - 2)] / 2
            diag(Q) <- c(-diag(Q[, -1])[seq_len(n_states - 2)] - diag(Q[, -c(1, 2)]),
                         -diag(Q[, -1])[n_states - 1],
                         0)
        } else if (n_states == 3) {
            Q[1, 3] <- 1 / sojourn_times[1] / 2
        }
        
    } else {
        stop("Only orders 1 and 2 are supported")
    }
    
    return(Q)
}


create_2nd_order_FC_intensity_matrix_with_likely_skips <- function(sojourn_times) {
    # This is an experimental function for creating a 2nd-order forward-chain
    # intensity matrix in which there is a strong chance of skipping some states.
    # This is to explore whether there are local optimum solutions for the
    # intensity matrix.
    
    # Start with usual 2nd-order FC intensity matrix and swap things around
    Q <- create_FC_intensity_matrix_from_sojourn_times(sojourn_times, order = 2)
    
    # Make probability of 2 -> 4 higher
    Q[2, 4] <- Q[2, 3]
    Q[2, 3] <- Q[2, 4] / 2   
    
    # Make probability of 3 -> 4 very low
    Q[3, 5] <- Q[3, 4]
    Q[3, 4] <- Q[3, 5] / 1e5
    Q[3, 3] <- -Q[3, 4] - Q[3, 5]
    
    return(Q)
}


create_2nd_order_FC_intensity_matrix_with_3_4_forbidden <- function(sojourn_times) {
    # This is an experimental function for creating a 2nd-order forward-chain
    # intensity matrix in which the transition between states 3 and 4 is
    # forbidden. This is to explore the possibility of parallel states,
    # i.e. a pair of states that are mutually exclusive for each patient.
    
    # Start with usual 2nd-order FC intensity matrix and swap things around
    Q <- create_FC_intensity_matrix_from_sojourn_times(sojourn_times, order = 2)
    
    # Make probability of 3 -> 4 zero
    Q[3, 5] <- Q[3, 4] + Q[3, 5]
    Q[3, 4] <- 0
    
    return(Q)
}