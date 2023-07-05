from __future__ import annotations

import numbers
import typing

import attr
import numpy
import pandas
import scipy.stats
from joblib import Parallel, delayed

from .model import distribution as dist
from .model import process


@attr.s(auto_attribs=True)
class Progression:
    stages: process.JumpProcess

    def simulate(
        self, random_state=None, subjid: str = None, mean_age: float = 0.0,
        sd_age: float = 1.0, **kwargs
    ) -> pandas.DataFrame:
        random_state = numpy.random.default_rng(random_state)

        # Draw a random number for later
        age_offset = mean_age + sd_age * random_state.standard_normal()

        # Create series
        series = self.stages.simulate_series(time_start=0, random_state=random_state)

        # Offset time vector to randomize position of t=0.
        # Entry into first state corresponds to birth,
        # and t=0 corresponds to enrolment in a study,
        # so assume a normal distribution of ages at baseline.
        series.index = series.index - age_offset

        # Convert to data frame
        out = series.reset_index(drop=False)
        if subjid:
            out["subjid"] = subjid

        return out

    @classmethod
    def from_Q(cls, Q: numpy.ndarray, **kwargs):
        stages = process.ContinuousTimeMarkovChain.from_partial_Q(
            Q, name="stage", time_name="date"
        )
        return cls(stages=stages, **kwargs)


@attr.s(auto_attribs=True)
class Observations:
    visit_interval: dist.DistributionType
    time_max: float
    emission_means: numpy.ndarray = numpy.array(
        []
    )  # Rows correspond to variables, columns to states
    emission_covs: numpy.ndarray = numpy.array([])
    emission_names: list = tuple()

    def simulate(
        self, subj_data: pandas.DataFrame, random_state=None, **kwargs
    ) -> pandas.DataFrame:
        random_state = numpy.random.default_rng(random_state)

        # TODO: Randomise time_max
        time_max = self.time_max

        visit_process = self.visit_process(random_state=random_state)
        visits = visit_process.simulate_series(
            time_start=0,
            stopping=lambda time, stage: time >= time_max,
            random_state=random_state,
        )

        # TODO: Remove random fraction of visits to model missed visits
        # Never remove the first one

        # Calculate subject's age at baseline
        age_at_baseline = -subj_data["date"].values[0]

        # Sample progression data at visit dates
        out = (
            subj_data.set_index("date")
            .reindex(visits.index, method="ffill")
            .join(visits)
            .reset_index(drop=False)
        )

        out = out.assign(age=out["date"] + age_at_baseline)

        if len(self.emission_names) > 0:
            # Simulate emissions for each visit
            emissions = numpy.zeros((out.shape[0], len(self.emission_names)))
            for i_row in range(out.shape[0]):
                emissions[i_row, :] = scipy.stats.multivariate_normal.rvs(
                    mean=self.emission_means[:, out.stage.values[i_row]],
                    cov=self.emission_covs[out.stage.values[i_row], :, :],
                    random_state=random_state,
                )

        if len(self.emission_names) > 0:
            for i_emission in range(len(self.emission_names)):
                out[self.emission_names[i_emission]] = emissions[:, i_emission]

        return out

    def visit_process(self, **kwargs):
        return process.JumpProcess(
            holding_time=process.FixedDistribution(self.visit_interval),
            next_state=process.FixedJumps(),
            name="seq",
            time_name="date",
        )


@attr.s(auto_attribs=True)
class IntervalImputer:
    start_col: str = "time_min"
    end_col: str = "time_max"
    time_col: str = "time"
    event_col: str = "event"
    sample: typing.Optional[numbers.Real] = None
    how: typing.Union[str, numbers.Real] = "uniform"

    def transform(
        self, data: pandas.DataFrame, *, random_state=None
    ) -> pandas.DataFrame:
        """
        Convert interval censored data to right-censored via imputation.

        Left-censored data are not supported.
        """

        # Resample data (if requested)
        if self.sample:
            random_state = numpy.random.default_rng(random_state)
            # Note that you can't pass a Generator directly to DataFrame.sample
            # https://github.com/pandas-dev/pandas/issues/38100
            data = data.sample(
                frac=self.sample, replace=True, random_state=random_state.bit_generator
            )

        # Extract start and end times
        time1 = data.loc[:, self.start_col]
        time2 = data.loc[:, self.end_col]

        # Check for left-censored data
        if not numpy.isfinite(time1).all():
            raise ValueError("Left-censored data not supported")

        # Type of event
        #   False = right-censored
        #   True = event at time
        event = numpy.isfinite(time2).rename(self.event_col)

        # Time of event
        if isinstance(self.how, numbers.Real):
            u = self.how
        elif self.how == "uniform":
            # Note that `random_state.uniform(low=time1, high=time2)` doesn't work
            random_state = numpy.random.default_rng(random_state)
            u = random_state.uniform(size=time1.shape)
        else:
            methods = {"minimum": 0, "median": 0.5, "maximum": 1}
            u = methods[self.how]

        time = (time1 + u * (time2 - time1)).where(event, time1).rename(self.time_col)

        return pandas.concat([time, event], axis=1)


def prepare_simulation_parameters(subjids):
    ss_entropy = 280_550_565_495_846_231_969_507_303_885_667_512_834
    ss = numpy.random.SeedSequence(ss_entropy)

    n_subjs = len(subjids)
    params_all = [
        {"subjid": s, "random_state": seed, "mean_age": 0, "sd_age": 0}
        for s, seed in zip(subjids, ss.spawn(n_subjs))
    ]

    return params_all


def create_synthetic_disease_progression_data(sim, params_all):
    state_true = Parallel(n_jobs=-1, verbose=1)(
        delayed(sim.simulate)(**params)
        for params in params_all
    )
    state_true = pandas.concat(state_true, ignore_index=True)

    return(state_true)


def get_intensity_matrix_for_delayed_state_entry(Q, state, t_delay, t_LE_max=81.0):
    # Q is the original intensity matrix. We assume it's a first-order forward-chain progression.
    # state is the state that we are delaying entry into (n.b. 0-indexed)
    # t_delay is the amount we are delaying it by (years)
    # t_LE_max is the maximum life expectancy that the model may have.
    # The default is approximately the UK's life expectancy
    
    Q_delayed = Q.copy()
    Q_delayed[state - 1, state] = 1 / (1 / Q_delayed[state - 1, state] + t_delay)

    Q_delayed[state - 1, state - 1] = -Q_delayed[state - 1, state]

    if sum(1 / numpy.diag(Q, 1)) > t_LE_max:
        Q_delayed[-2, -1] = 1 / (t_LE_max - sum(1 / numpy.diag(Q_delayed[:-1, :-1], 1)))
        Q_delayed[-2, -2] = -Q_delayed[-2, -1]
        if Q_delayed[-2, -1] < 0:
            raise ValueError('Cannot delay state transition this much with current model')

    return Q_delayed 
