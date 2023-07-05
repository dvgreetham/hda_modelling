from hem import simulation
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pyprojroot import here
from scipy.stats import gamma

state_pension_age = 68
gdp_name = "Lower contribution to GDP"
MONTHS_PER_YEAR = 12


def generate_patients_with_costs(Q, target_means):

    # Set up progression simulation
    sim = simulation.Progression.from_Q(Q)

    # Create list of subjects
    n_patients = 10_000
    subjids = [f"examp{i:07d}" for i in range(n_patients)]

    params_all = simulation.prepare_simulation_parameters(subjids)

    state_true = simulation.create_synthetic_disease_progression_data(sim, params_all)

    state_true = state_true.rename(columns={"date": "age"})

    (
        age_vec,
        costs_by_age,
        cost_names,
        costs_over_lifetime,
    ) = assign_costs_to_patients(subjids, state_true, target_means)

    return age_vec, costs_by_age, cost_names, costs_over_lifetime


def assign_costs_to_patients(subjids, state_true, target_means):
    n_patients = len(subjids)
    np.random.seed(20220520)
    patient_quantile = np.random.uniform(size=n_patients)
    quantile_data = pd.DataFrame({"subjid": subjids, "quantile": patient_quantile})

    cost_data = state_true.merge(quantile_data, on="subjid", how="left").rename(
        columns={"stage": "state"}
    )

    cost_params = pd.read_csv(here() / "models" / "total_cost_MLE_fit.csv")

    # Get total monthly costs for each patient for hidden states
    cost_data["total_cost_per_month"] = 0
    n_hidden_states = cost_params.shape[0]

    for idx in cost_data.index:
        if cost_data.loc[idx, "state"] < n_hidden_states:
            shape = cost_params.loc[cost_data.loc[idx, "state"], "alpha_MLE"]
            scale = 1 / cost_params.loc[cost_data.loc[idx, "state"], "lambda_MLE"]
            cost_data.loc[idx, "total_cost_per_month"] = gamma.ppf(
                cost_data.loc[idx, "quantile"], a=shape, scale=scale
            )
            # else leave at zero

    cost_groups = {
        "PIP": [
            "cost PIP daily living standard",
            "cost PIP daily living enhanced",
            "cost PIP mobility standard",
            "cost PIP mobility enhanced",
        ],
        "Universal Credit + ESA": [
            "cost Employment support allowance - low",
            "cost Employment support allowance - high",
            "cost Universal credit - children",
            "cost Universal credit - standard",
        ],
        "Partner stops working/care worker/care home": [
            "cost Care home costs - NHS continuing care",
            "cost Nursing in care home",
            "cost Partner stops working",
            "funding_for_care_worker_or_home_per_month",
        ],
        "Lower contribution to GDP": ["Lower contribution to GDP"],
        "Pharmacologic and non-pharmacologic therapies": [
            "total_expected_cost_per_month_pharma",
            "total_expected_cost_per_month_nonpharma",
        ],
    }

    target_means_grouped = pd.concat(
        [
            pd.DataFrame(target_means.loc[cols].sum(axis=0).to_dict(), index=[grp])
            for (grp, cols) in cost_groups.items()
        ]
    )

    idx_gdp = (target_means_grouped.index == gdp_name).nonzero()

    indiv_costs, total_costs = get_individual_costs_from_total_costs(
        cost_data["total_cost_per_month"].values,
        target_means_grouped,
        cost_params["lambda_MLE"].values,
        cost_data["state"].values,
        idx_gdp,
        cost_data["age"].values,
    )

    # Include loss of GDP due to death in total costs
    cost_data["total_cost_per_month"] = total_costs

    cost_names = list(indiv_costs.columns)

    cost_data = pd.merge(cost_data, indiv_costs, left_index=True, right_index=True)

    age_vec = np.linspace(18, 100, 83)
    max_age = age_vec.max()

    costs_by_age = Parallel(n_jobs=-1, verbose=1)(
        delayed(get_costs_per_month_at_age)(cost_data, age_vec[i_age])
        for i_age in range(len(age_vec))
    )

    costs_by_age = pd.concat(costs_by_age)

    costs_by_age["total_cost_per_year"] = (
        costs_by_age["total_cost_per_month"] * MONTHS_PER_YEAR
    )

    costs_over_lifetime = costs_by_age.groupby("subjid").apply(
        lambda s: pd.Series(
            {
                "mean_cost_per_year": np.trapz(s["total_cost_per_year"], s["age"])
                / max_age,
                "mean_sq_cost_per_year": np.trapz(
                    s["total_cost_per_year"] ** 2, s["age"]
                )
                / max_age,
            }
        )
    )

    costs_over_lifetime["var_cost_per_year"] = (
        costs_over_lifetime["mean_sq_cost_per_year"]
        - costs_over_lifetime["mean_cost_per_year"] ** 2
    )

    return age_vec, costs_by_age, cost_names, costs_over_lifetime


def calculate_mean_and_variance_of_cost_saving(costs_by_age, costs_by_age_delayed):
    max_age = costs_by_age["age"].max()

    subjids = costs_by_age["subjid"]
    ages = costs_by_age["age"]
    cols_to_drop = ["subjid", "age", "state", "quantile"]
    cost_saving = costs_by_age.drop(columns=cols_to_drop) - costs_by_age_delayed.drop(
        columns=cols_to_drop
    )
    cost_saving["subjid"] = subjids
    cost_saving["age"] = ages

    grouped_summary_stats = cost_saving.groupby("subjid").apply(
        lambda s: pd.Series(
            {
                "mean_cost_saving_per_year": np.trapz(
                    s["total_cost_per_year"], s["age"]
                )
                / max_age,
                "mean_sq_cost_saving_per_year": np.trapz(
                    s["total_cost_per_year"] ** 2, s["age"]
                )
                / max_age,
            }
        )
    )

    grouped_summary_stats["var_cost_saving_per_year"] = (
        grouped_summary_stats["mean_sq_cost_saving_per_year"]
        - grouped_summary_stats["mean_cost_saving_per_year"] ** 2
    )

    mean_cost_saving_per_patient = grouped_summary_stats[
        "mean_cost_saving_per_year"
    ].mean()

    var_cost_saving_per_patient = (
        grouped_summary_stats["mean_cost_saving_per_year"].var()
        + grouped_summary_stats["var_cost_saving_per_year"].mean()
    )

    return mean_cost_saving_per_patient, var_cost_saving_per_patient


def get_costs_per_month_at_age(cost_data, age):

    costs_per_month_at_age = (
        cost_data.loc[cost_data["age"] <= age]
        .groupby("subjid")
        .tail(1)
        .reset_index(drop=True)
    )

    costs_per_month_at_age["age"] = age

    # Can't get universal credit or ESA if over state pension age,
    # and wouldn't be contributing to gdp
    idx_over_pension = costs_per_month_at_age["age"] >= state_pension_age
    UC_ESA_GDP_cost_names = ["Universal Credit + ESA", gdp_name]
    for cost_name in UC_ESA_GDP_cost_names:
        costs_per_month_at_age.loc[
            idx_over_pension, "total_cost_per_month"
        ] += -costs_per_month_at_age.loc[idx_over_pension, cost_name]

    costs_per_month_at_age.loc[idx_over_pension, UC_ESA_GDP_cost_names] = 0

    return costs_per_month_at_age


def get_individual_costs_from_total_costs(
    total_costs, target_means, beta, i_state, idx_gdp, age
):
    n_costs = target_means.shape[0]
    n_hidden_states = target_means.shape[1]
    (idx_death,) = (i_state == n_hidden_states).squeeze().nonzero()
    (idx_alive,) = (i_state < n_hidden_states).squeeze().nonzero()

    n_samples = len(total_costs)

    np.random.seed(20220526)

    alpha_indiv = np.tile(beta, (n_costs, 1)) * target_means.to_numpy()

    gammas = np.zeros((n_costs, n_samples))

    for i_cost in range(n_costs):
        gammas[i_cost, idx_alive] = gamma.rvs(
            a=alpha_indiv[i_cost, i_state[idx_alive]],
            scale=1 / beta[i_state[idx_alive]],
        )

    dirichlets = gammas / np.tile(np.sum(gammas, axis=0), (n_costs, 1))

    indiv_costs = (dirichlets * np.tile(total_costs, (n_costs, 1))).T

    # Set cost due to lost contribution to GDP in death
    # equal to the cost for the last living stage,
    # if the death age was below state pension age,
    # 0 otherwise
    for i in idx_death:
        indiv_costs[i, :] = 0
        if age[i] < state_pension_age:
            indiv_costs[i, idx_gdp] = indiv_costs[i - 1, idx_gdp]
            total_costs[i] = indiv_costs[i, idx_gdp]
        # else leave at zero

    indiv_costs = pd.DataFrame(
        indiv_costs, index=range(indiv_costs.shape[0]), columns=target_means.index
    )

    return indiv_costs, total_costs
