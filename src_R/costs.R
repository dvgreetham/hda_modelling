library("dplyr")
library("magrittr")
library("foreach")
library("doParallel")


state_pension_age <- 68


# function for calculating HSC costs
estimate_HSC_costs_incurred <-
    function(visit_data,
             benefits_thresholds,
             quantile_data) {
        # quantile_data columns contain
        # random samples from U(0, 1) distributions,
        # representing what quantile this patient is at
        # in terms of speed of benefit uptake
        # and baseline economic circumstances (for means-testing purposes).
        # These indicate the inter-patient variability
        # in economics and benefit uptake, but remain constant for
        # an individual patient.
        
        # Set up comparator functions for determining how the thresholds are applied
        # on age, TFC etc.
        comparators <- tribble(
            ~ variable,
            ~ threshold_name,
            ~ comparator,
            "tfc",
            "tfc_threshold",
            `<=`,
            "occupatn",
            "occupatn_upper",
            `<=`,
            "carehome",
            "carehome_value",
            `==`,
            "carelevl",
            "carelevl_value",
            `==`,
            "age",
            "age_max",
            `<`,
            "is_living_with_partner",
            "is_living_with_partner",
            `==`,
            "partner_will_stay_at_work",
            "partner_will_stay_at_work",
            `==`
        )
        
        # Make sure data is in expected order
        visit_data %<>% arrange(subjid, visdy)
        
        # Randomly assign partners if maristat is missing
        p_partner <- visit_data %>%
            drop_na(maristat) %>%
            group_by(subjid) %>%
            summarise(is_with_partner =
                          sum(maristat %in% c("partnership", "married")) / n() > 0.5,
                      .groups = "drop_last") %>%
            pull(is_with_partner) %>%
            mean()
        
        set.seed(20220506)
        visit_data %<>%
            group_by(subjid) %>%
            mutate(maristat = ifelse(
                is.na(maristat),
                ifelse(rbinom(1, 1, p_partner) == 1, "partnership", "single"),
                maristat
            )) %>%
            ungroup()
        
        benefits_thresholds %<>%
            mutate(age_max = ifelse(under_state_pension_age,
                                    state_pension_age,
                                    NA))
        
        care_worker_hours_max <- 48  # per week
        
        # We identify two scenarios in which a care worker comes in:
        # Care comes in 1: Patient is not living with a partner but needs care
        # Care comes in 2: Patient is living with a partner but the partner stays at work.
        # Here we set up the parameters for a linear increase in the number of
        # hours of a care worker needed as TFC decreases,
        # up to a maximum number of hours at a TFC of 1.
        # (N.b. patients with a TFC of 0 are definitely in a nursing home,
        # so don't use care workers.)
        # We use the upper limit of TFC for "Care comes in 1" --
        # could equally have used "Care comes in 2" as the TFC limit is the same.
        
        tfc_care_upper <- benefits_thresholds %>%
            filter(name == "Care comes in 1") %>%
            pull(tfc_upper)
        tfc_max_care <- 1
        
        n_cores <- detectCores()
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)
        
        visit_data <- foreach(
            subj = unique(visit_data$subjid),
            .combine = bind_rows,
            .multicombine = TRUE,
            .packages = "dplyr",
            .export = "is_patient_receiving_benefit_at_this_visit"
        ) %dopar% {
            quantiles_this_patient <- quantile_data %>%
                filter(subjid == subj)
            
            q_patient_benefit_uptake <-
                quantiles_this_patient$q_patient_benefit_uptake
            
            # Calculate patient-specific TFC thresholds
            for (i_benefit in seq_along(benefits_thresholds$name)) {
                benefits_thresholds[i_benefit, 'tfc_threshold'] <-
                    q_patient_benefit_uptake * benefits_thresholds[i_benefit, 'tfc_upper'] +
                    (1 - q_patient_benefit_uptake) * benefits_thresholds[i_benefit, 'tfc_lower']
            }
            
            visit_data_this_subj <-
                visit_data %>% filter(subjid == subj)
            
            for (i_visit in seq_len(nrow(visit_data_this_subj))) {
                comparators %<>%
                    select(-any_of("value")) %>%
                    left_join(
                        tribble(
                            ~ variable,
                            ~ value,
                            "age",
                            visit_data_this_subj[i_visit, "age"],
                            "tfc",
                            visit_data_this_subj[i_visit, "tfcscore"],
                            "occupatn",
                            visit_data_this_subj[i_visit, "occupatn"],
                            "carehome",
                            visit_data_this_subj[i_visit, "carehome"],
                            "carelevl",
                            visit_data_this_subj[i_visit, "carelevl"],
                            "is_living_with_partner",
                            visit_data_this_subj[["maristat"]][i_visit] %in%
                                c("partnership", "married"),
                            "partner_will_stay_at_work",
                            quantiles_this_patient$partner_will_stay_at_work
                        ),
                        by = "variable"
                    )
                
                visit_data_this_subj[i_visit, 'hours_of_care_needed_per_week'] <-
                    pmin(pmax(
                        0,
                        care_worker_hours_max *
                            (tfc_care_upper - visit_data_this_subj[i_visit, "tfcscore"]) /
                            (tfc_care_upper - tfc_max_care)
                    ),
                    care_worker_hours_max)
                
                for (i_benefit in seq_along(benefits_thresholds$name)) {
                    benefit <- benefits_thresholds$name[i_benefit]
                    visit_data_this_subj[i_visit, benefit] <-
                        is_patient_receiving_benefit_at_this_visit(
                            benefits_thresholds,
                            i_benefit,
                            comparators,
                            visit_data_this_subj,
                            i_visit
                        )
                }
            }
            
            visit_data_this_subj
        }
        
        stopCluster(cl)
        
        # If getting enhanced PIP, can't also get standard PIP
        visit_data %<>%
            mutate(
                `PIP daily living standard` = `PIP daily living standard` &
                    !`PIP daily living enhanced`,
                `PIP mobility standard` = `PIP mobility standard` &
                    !`PIP mobility enhanced`
            )
        
        # If getting High ESA, can't also get low ESA
        visit_data %<>%
            mutate(`Employment support allowance - low` =
                       `Employment support allowance - low` &
                       !`Employment support allowance - high`)
        
        return(visit_data)
    }


plot_TFC_cost_assumptions <- function(benefit_thresholds) {
    plot_data <- benefit_thresholds %>%
        pivot_longer(
            cols = c(tfc_upper, tfc_lower),
            names_to = "threshold",
            values_to = "tfc"
        ) %>%
        mutate(prop_people = ifelse(threshold == "tfc_upper", 0, 1)) %>%
        select(-threshold) %>%
        bind_rows(
            tibble(
                name = benefit_thresholds$name,
                tfc = 0,
                prop_people = 1
            ),
            tibble(
                name = benefit_thresholds$name,
                tfc = 13,
                prop_people = 0
            )
        ) %>%
        arrange(name, tfc)
    
    threshold_plot <-
        ggplot(plot_data, aes(x = tfc, y = prop_people, colour = name)) +
        geom_line(size = 1.5) +
        labs(x = "TFC score", y = "Proportion of people claiming") +
        scale_x_continuous(breaks = 0:13) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        theme(text = element_text(size = 20))
    
    return(threshold_plot)
}


assign_costs_to_cost_items <-
    function(visit_data, cost_data, quantile_data) {
        # Given the logical indicators showing which costs apply in visit_data,
        # determine the value of those costs, and the total cost per month
        # at each patient visit.
        days_per_year_avg <- 365.25
        days_per_week <- 7
        months_per_year <- 12
        weeks_per_year_avg <- days_per_year_avg / days_per_week
        weeks_per_month_avg <- weeks_per_year_avg / months_per_year
        
        
        unit_scale_factor <- c(
            "week" = weeks_per_month_avg,
            "month" = 1,
            "year" = 1 / months_per_year
        )
        
        # Calculate socio-economic-dependent contributions
        n_subjs <- n_distinct(visit_data$subjid)
        salary <- numeric(n_subjs)
        quantile_data %<>% arrange(subjid)
        salary_percentile <- floor(quantile_data$q_socioeconomic * 10) * 10
        age_summary <- visit_data %>%
            group_by(subjid) %>%
            summarise(age_at_baseline = min(age), .groups = "drop_last") %>%
            mutate(age_group = cut(age_at_baseline,
                                   breaks = c(18, 22, 30, 40, 50, 60, Inf),
                                   right = FALSE,
                                   labels = c("18-21", "22-29", "30-39",
                                              "40-49", "50-59", "60+"))) %>%
            arrange(subjid)
        
        quantile_data %<>%
            left_join(age_summary, by = "subjid")
        
        for (i_subj in seq_len(n_subjs)) {
            if (salary_percentile[i_subj] == 0) {
                salary[i_subj] <- 0
            } else {
                cost_item <- paste0("gross salary ",
                                    quantile_data$age_group[i_subj],
                                    " percentile ",
                                    salary_percentile[i_subj])
                salary[i_subj] <-
                    cost_data[['cost_typical_gbp']][cost_data$item == cost_item] *
                    unit_scale_factor[cost_data[['cost_per']][cost_data$item == cost_item]]
            }
        }
        
        GDP_out <- estimate_contributions_to_GDP(cost_data, salary,
                                                 quantile_data)
        contributions_to_GDP <- GDP_out$individual_GDP_contribution
        k <- GDP_out$k
        
        # Calculate carer's allowance and UC for caring for disabled
        carers_allowance <-
            cost_data[['cost_typical_gbp']][cost_data$item ==
                                                "carer's allowance"] *
            unit_scale_factor[cost_data[['cost_per']][cost_data$item ==
                                                          "carer's allowance"]]
        universal_credit_care_for_severely_disabled <-
            cost_data[['cost_typical_gbp']][cost_data$item ==
                                                "universal credit - care for severely disabled person"] *
            unit_scale_factor[cost_data[['cost_per']][cost_data$item ==
                                                          "universal credit - care for severely disabled person"]]
        
        
        
        # Fixed costs - costs that do not vary either intra-patient or inter-patient
        fixed_costs_mapping <- get_fixed_costs_mappings(cost_data)
        
        # Inter-patient variable costs --
        # for these we assume a fixed amount for each patient,
        # but it can vary between patients.
        # N.b. whether the cost is incurred can still vary between a patient's
        # visits.
        living_with_partner <- visit_data %>%
            group_by(subjid) %>%
            summarise(
                is_living_with_partner =
                    sum(maristat %in% c("partnership", "married")) / n() > 0.5,
                .groups = "drop_last"
            )
        
        quantile_data %<>% left_join(living_with_partner, by = "subjid")
        
        inter_patient_variable_costs_mapping <- tibble(
            subjid = quantile_data$subjid,
            `Partner stops working` = contributions_to_GDP + carers_allowance +
                universal_credit_care_for_severely_disabled,
            `Contribution to GDP halved` = 0.5 * contributions_to_GDP,
            `No contribution to GDP` = 0.5 * contributions_to_GDP,  # N.b. cumulative with contribution to GDP halved
            `Universal credit - children` =
                cost_data[['cost_low_gbp']][cost_data$item == "universal credit - children"] *
                unit_scale_factor[cost_data[['cost_per']][cost_data$item == "universal credit - children"]] *
                quantile_data$n_children,
            # Assume almost everyone is over 25, so use cost_high_gbp
            `Universal credit - standard` =
                ifelse(
                    quantile_data$is_living_with_partner,
                    cost_data[['cost_high_gbp']][cost_data$item == "universal credit - standard, couple"] *
                        unit_scale_factor[cost_data[['cost_per']][cost_data$item == "universal credit - standard, couple"]],
                    cost_data[['cost_high_gbp']][cost_data$item == "universal credit - standard, single"] *
                        unit_scale_factor[cost_data[['cost_per']][cost_data$item == "universal credit - standard, single"]]
                )
        )
        
        inter_patient_variable_costs <-
            setdiff(names(inter_patient_variable_costs_mapping), "subjid")
        
        # Intra-patient variable costs
        low <-
            cost_data %>% filter(item == "home care") %>% pull(cost_low_gbp)
        mid <-
            cost_data %>% filter(item == "home care") %>% pull(cost_typical_gbp)
        high <-
            cost_data %>% filter(item == "home care") %>% pull(cost_high_gbp)
        cost_care_hour_expected <-
            (low + 4 * mid + high) / 6  # 3-point estimate
        visit_data %<>%
            mutate(
                cost_of_care_worker_per_month =
                    ifelse(
                        `Care comes in 1` | `Care comes in 2`,
                        hours_of_care_needed_per_week * cost_care_hour_expected *
                            unit_scale_factor['week'],
                        0
                    ),
                cost_of_care_home_per_month =
                    ifelse(
                        `Local authority funding for care home`,
                        cost_data[["cost_typical_gbp"]][cost_data$item ==
                                                            "local authority own-provision care home"] *
                            unit_scale_factor[cost_data[["cost_per"]][cost_data$item ==
                                                                          "local authority own-provision care home"]],
                        0
                    )
            )
        
        # Change all units to per month
        fixed_costs_mapping$scale_factor <-
            unit_scale_factor[fixed_costs_mapping$cost_per]
        fixed_costs_mapping %<>%
            mutate(cost_per_month = round(cost * scale_factor, 2))
        
        # Attach fixed costs to visit data
        for (cost_item in fixed_costs_mapping$cost_item) {
            visit_data[[paste0('cost ', cost_item)]] <-
                ifelse(visit_data[[cost_item]],
                       fixed_costs_mapping$cost_per_month[fixed_costs_mapping$cost_item == cost_item],
                       0)
        }
        
        # Attach inter-patient variable costs to visit data
        for (i_visit in seq_len(nrow(visit_data))) {
            subj <- visit_data$subjid[i_visit]
            for (cost_name in inter_patient_variable_costs) {
                visit_data[i_visit, paste0("cost ", cost_name)] <-
                    ifelse(
                        visit_data[i_visit, cost_name],
                        (
                            inter_patient_variable_costs_mapping %>%
                                filter(subjid == subj) %>%
                                pull(cost_name)
                        ),
                        0
                    )
            }
        }
        
        # PIP is intended to help cover costs of care workers/care homes,
        # so local authority funding will be less as a result
        visit_data$funding_for_care_worker_or_home_per_month <-
            NA_real_
        visit_data %<>%
            group_by(subjid) %>%
            mutate(
                funding_for_care_worker_or_home_per_month =
                    pmax(
                        0,
                        (1 - quantile_data[quantile_data$subjid == first(subjid),
                                           'q_socioeconomic']) *
                            pmax(
                                cost_of_care_worker_per_month,
                                cost_of_care_home_per_month
                            )
                        -
                            `cost PIP daily living standard` -
                            `cost PIP mobility standard` -
                            `cost PIP daily living enhanced` -
                            `cost PIP mobility enhanced`
                    )
            ) %>%
            ungroup()
        
        cost_names <- c(paste0(
            "cost ",
            c(
                fixed_costs_mapping$cost_item,
                inter_patient_variable_costs
            )
        ),
        "funding_for_care_worker_or_home_per_month",
        "total_expected_cost_per_month_pharma",
        "total_expected_cost_per_month_nonpharma")
        
        visit_data[['total_costs_per_month']] <-
            rowSums(visit_data[, paste0("cost ",
                                        c(
                                            fixed_costs_mapping$cost_item,
                                            inter_patient_variable_costs
                                        ))]) +
            visit_data$funding_for_care_worker_or_home_per_month +
            visit_data$total_expected_cost_per_month_pharma +
            visit_data$total_expected_cost_per_month_nonpharma
        
        
        
        return(list(visit_data = visit_data,
                    cost_names = cost_names,
                    k = k))
    }


get_fixed_costs_mappings <- function(cost_data) {
    fixed_costs_mapping <- tribble(
        ~ cost_item,
        ~ cost,
        ~ cost_per,
        "PIP daily living standard",
        cost_data[['cost_low_gbp']][cost_data$item ==
                                        "personal independence payments - daily living"],
        cost_data[['cost_per']][cost_data$item ==
                                    "personal independence payments - daily living"],
        "PIP daily living enhanced",
        cost_data[['cost_high_gbp']][cost_data$item ==
                                         "personal independence payments - daily living"],
        cost_data[['cost_per']][cost_data$item ==
                                    "personal independence payments - daily living"],
        "PIP mobility standard",
        cost_data[['cost_low_gbp']][cost_data$item ==
                                        'personal independence payments - mobility'],
        cost_data[['cost_per']][cost_data$item ==
                                    'personal independence payments - mobility'],
        "PIP mobility enhanced",
        cost_data[['cost_high_gbp']][cost_data$item ==
                                         'personal independence payments - mobility'],
        cost_data[['cost_per']][cost_data$item ==
                                    'personal independence payments - mobility'],
        "Care home costs - NHS continuing care",
        cost_data[['cost_typical_gbp']][cost_data$item ==
                                            'local authority own-provision care home'],
        cost_data[['cost_per']][cost_data$item ==
                                    'local authority own-provision care home'],
        "Nursing in care home",
        cost_data[['cost_typical_gbp']][cost_data$item ==
                                            "nursing in care home"],
        cost_data[['cost_per']][cost_data$item ==
                                    "nursing in care home"],
        "Employment support allowance - low",
        cost_data[['cost_low_gbp']][cost_data$item == "employment and support allowance"],
        cost_data[['cost_per']][cost_data$item == "employment and support allowance"],
        "Employment support allowance - high",
        cost_data[['cost_high_gbp']][cost_data$item == "employment and support allowance"],
        cost_data[['cost_per']][cost_data$item == "employment and support allowance"],
    )
    
    return(fixed_costs_mapping)
}


is_patient_receiving_benefit_at_this_visit <-
    function(benefits_thresholds,
             i_benefit,
             comparators,
             visit_data_this_subj,
             i_visit) {
        # Based on the conditions required for each benefit given in benefits thresholds,
        # determine whether this subject is receiving benefit i_benefit at visit i_visit.
        
        benefit <- benefits_thresholds$name[i_benefit]
        
        meets_conditions <- logical(nrow(comparators))
        for (i_condition in seq_len(nrow(comparators))) {
            meets_conditions[i_condition] <-
                is.na(benefits_thresholds[i_benefit, comparators$threshold_name[i_condition]]) ||
                comparators$comparator[[i_condition]](comparators$value[i_condition],
                                                      benefits_thresholds[i_benefit,
                                                                          comparators$threshold_name[i_condition]])
        }
        
        # Determine which costs apply based on visit data
        if (i_visit == 1) {
            return(all(meets_conditions))
        } else {
            # Assume that once cost has been incurred once,
            # it will continue to be incurred,
            # unless it depends on the patient not being in a care home
            # (e.g. partner may go back to work if patient moves to care home)
            i_carehome <- which(comparators$variable == "carehome")
            return((visit_data_this_subj[[benefit]][i_visit - 1] ||
                        all(meets_conditions[-i_carehome])) &&
                       meets_conditions[i_carehome])
        }
    }


estimate_contributions_to_GDP <- function(cost_data, salary, quantile_data) {
    # Assume contribution to GDP is (1 + k) * salary,
    # where k * salary is the assumed contribution to savings generated by businesses.
    
    total_UK_GDP <- 2317054e6  # From https://www.ons.gov.uk/economy/grossdomesticproductgdp/datasets/uksecondestimateofgdpdatatables
                               # Sheet P, 2021 figure.
    n_subjs <- nrow(quantile_data)
    UK_working_population <- 23865e3  # From source 14 in data/raw/costs_of_care.xlsx

    expected_salary_enroll <- mean(salary)
    expected_salary_genpop <- cost_data[["cost_typical_gbp"]][cost_data$cost_item ==
                                                                  "gross salary all mean"]
    age_demographic_scale_factor <- expected_salary_enroll / expected_salary_genpop
    
    total_GDP_contribution_enroll <-
        total_UK_GDP *
        age_demographic_scale_factor *
        n_subjs / UK_working_population
    
    total_salary_enroll <- sum(salary)
    
    # Calculate constant of proportionality k
    k <- max(0, total_GDP_contribution_enroll / total_salary_enroll - 1)
    
    individual_GDP_contribution <- (1 + k) * salary
    
    return(list(individual_GDP_contribution = individual_GDP_contribution,
                k = k))
}