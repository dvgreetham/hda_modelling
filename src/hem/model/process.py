"""
Definitions here are focussed on simulation, so, internally,
they use the jump chain/holding time definition of a Markov Chain.
This also allows modelling of other kinds of processes with holding
times that have general, rather than just memoryless, distributions,
such as renewal-reward processes and certain types of queue.

[1]: J.R.Norris. Markov Chains. Cambridge University Press, 1997.
"""

# TODO: Consider improving with labelled stages, rather than numeric ones

from __future__ import annotations

import typing

import attr
import numpy
import pandas
import scipy

from . import distribution as dist

StateType = int
DistributionFunction = typing.Callable[[StateType], dist.DistributionType]


### Abstract processes
@attr.s(auto_attribs=True)
class JumpProcess:
    holding_time: DistributionFunction
    next_state: DistributionFunction
    initial_state: StateType = 0
    name: str = "state"
    time_name: str = "time"
    time_format: str = ".1f"

    def get_params(self) -> dict:
        return {
            "initial_state": self.initial_state,
            "name": self.name,
            "time_name": self.time_name,
            "time_format": self.time_format,
        }

    def simulate(
        self,
        initial_state: StateType = None,
        *,
        time_start: float = 0.0,
        stopping: typing.Callable[[float, StateType], bool] = lambda time, state: False,
        random_state=None,
        verbose: bool = False,
    ) -> typing.Iterator[typing.Tuple[float, StateType]]:
        random_state = numpy.random.default_rng(random_state)
        if initial_state is None:
            initial_state = self.initial_state

        time: float = time_start
        state: StateType = None
        delta: float = 0.0

        while numpy.isfinite(time):
            # Determine next state
            if state is None:
                state = initial_state
                if verbose:
                    print(f"{time:{self.time_format}}: start in {state}")
            else:
                state = self.next_state(state).rvs(random_state=random_state)
                if verbose:
                    print(
                        f"{time:{self.time_format}}: jump to {state} after {delta:{self.time_format}}"
                    )

            # Yield jump time and next state
            yield time, state

            # Determine next jump time
            delta = self.holding_time(state).rvs(random_state=random_state)
            time += delta

            if stopping(time, state):
                # Time of next jump is reached or target state is reached
                # or some combination
                break

        if verbose:
            if numpy.isfinite(time):
                print(f"{time:{self.time_format}}: stop in {state}")
            else:
                print(f"terminate in {state}")

    def simulate_series(self, *args, **kwargs) -> pandas.Series:
        index, values = zip(*self.simulate(*args, **kwargs))
        return pandas.Series(
            values, index=index, name=self.name, dtype=StateType
        ).rename_axis(self.time_name)


### Generic transition processes
@attr.s(auto_attribs=True, frozen=True)
class JumpMatrix:
    jump_matrix: numpy.ndarray = numpy.array([[]])

    def __call__(self, state: StateType) -> dist.DistributionType:
        return dist.categorical_distribution(tuple(self.jump_matrix[state, :].flat))


@attr.s(auto_attribs=True, frozen=True)
class FixedJumps:
    step: int = 1

    def __call__(self, state: StateType) -> dist.DistributionType:
        next_state = state + self.step
        return dist.delta_distribution(next_state)


### Generic holding time distributions
@attr.s(auto_attribs=True, frozen=True)
class FixedTime:
    timestep: float = 1.0

    def __call__(self, state: StateType) -> dist.DistributionType:
        return dist.delta_distribution(self.timestep)


@attr.s(auto_attribs=True, frozen=True)
class FixedDistribution:
    rv: dist.DistributionType

    def __call__(self, state: StateType) -> dist.DistributionType:
        return self.rv


@attr.s(auto_attribs=True, frozen=True)
class GeometricTime:
    jump_prob: float = 1.0

    def __call__(self, state: StateType) -> dist.DistributionType:
        return dist.geometric_distribution(p=self.jump_prob)


@attr.s(auto_attribs=True, frozen=True)
class ExponentialTime:
    jump_rate: float = 1.0

    def __call__(self, state: StateType) -> dist.DistributionType:
        return dist.exponential_distribution(rate=self.jump_rate)


@attr.s(auto_attribs=True, frozen=True)
class AbsorbingStates:
    """
    Override an arbitrary holding-time distribution with fixed absorbing states.
    """

    holding_time: DistributionFunction
    absorbing_states: typing.Set[StateType] = attr.Factory(set)

    def __call__(self, state: StateType) -> dist.DistributionType:
        if state in self.absorbing_states:
            return dist.delta_distribution(numpy.inf)
        else:
            return self.holding_time(state)


@attr.s(auto_attribs=True, frozen=True)
class DiscreteMarkov:
    jump_prob: numpy.ndarray = numpy.array([])

    def __call__(self, state: StateType) -> dist.DistributionType:
        return dist.geometric_distribution(p=self.jump_prob[state])


@attr.s(auto_attribs=True, frozen=True)
class ContinuousMarkov:
    jump_rate: numpy.ndarray = numpy.array([])

    def __call__(self, state: StateType) -> dist.DistributionType:
        return dist.exponential_distribution(rate=self.jump_rate[state])


### Specific processes
@attr.s(auto_attribs=True)
class PoissonProcess(JumpProcess):
    @classmethod
    def from_rate(cls, rate: float, **kwargs):
        return cls(
            next_state=FixedJumps(1), holding_time=ExponentialTime(rate), **kwargs
        )


@attr.s(auto_attribs=True)
class FixedTimeMarkovChain(JumpProcess):
    @classmethod
    def from_P(cls, P: numpy.ndarray, timestep: float = 1, **kwargs):
        return cls(next_state=JumpMatrix(P), holding_time=FixedTime(timestep), **kwargs)


@attr.s(auto_attribs=True)
class DiscreteTimeMarkovChain(JumpProcess):
    @classmethod
    def from_P(cls, P: numpy.ndarray, **kwargs):
        if not P.ndim == 2:
            raise ValueError("Input array must be 2-d")
        if not numpy.alltrue(numpy.diff(P.shape) == 0):
            raise ValueError("All dimensions of input must be of equal length")
        if not numpy.alltrue(P >= 0):
            raise ValueError("Elements of input array must all be non-negative")
        if not numpy.allclose(P.sum(axis=1), 1):
            raise ValueError("Row sums of input array must all be 1")

        # Shortcut for square matrix
        diag_slice = slice(None, None, P.shape[1] + 1)

        # Calculate jump probabilities
        # = 1 - p_ii
        jump_prob = 1 - P.flat[diag_slice].copy()
        jump_prob[jump_prob <= 0] = 0  # Fix negative zero values (and -eps)

        # Create jump matrix
        jump_matrix = P.copy()
        jump_matrix.flat[diag_slice] = 0
        row_sums = jump_matrix.sum(axis=1)
        jump_matrix /= numpy.where(row_sums <= 0, 1, row_sums)[:, None]
        jump_matrix.flat[diag_slice] = row_sums == 0

        return cls(
            next_state=JumpMatrix(jump_matrix),
            holding_time=DiscreteMarkov(jump_prob),
            **kwargs,
        )

    @property
    def jump_matrix(self) -> numpy.ndarray:
        return self.next_state.jump_matrix

    @property
    def jump_prob(self) -> numpy.ndarray:
        return self.holding_time.jump_prob


@attr.s(auto_attribs=True)
class ContinuousTimeMarkovChain(JumpProcess):
    @classmethod
    def from_Q(cls, Q: numpy.ndarray, **kwargs):
        # Based on section 2.6 of [1]
        if not Q.ndim == 2:
            raise ValueError("Input array must be 2-d")
        if not numpy.alltrue(numpy.diff(Q.shape) == 0):
            raise ValueError("All dimensions of input must be of equal length")
        if not numpy.allclose(Q.sum(axis=1), 0):
            raise ValueError("Row sums of input array are not all 0")

        # Shortcut for square matrix
        diag_slice = slice(None, None, Q.shape[1] + 1)

        # Extract jump rates
        # = -q_ii
        jump_rate = -Q.flat[diag_slice].copy()
        jump_rate[jump_rate == 0] = 0  # Fix negative zero values
        if not numpy.alltrue(jump_rate >= 0):
            raise ValueError(
                "Diagonal elements of input array must all be non-positive"
            )

        # Create jump matrix
        #   See p87 of [1]
        jump_matrix = Q.copy()
        jump_matrix /= numpy.where(jump_rate == 0, 1, jump_rate)[:, None]
        jump_matrix.flat[diag_slice] = jump_rate == 0
        if not numpy.alltrue(jump_matrix >= 0):
            raise ValueError(
                "Off-diagonal elements of input array must all be non-negative"
            )

        # Done
        return cls(
            next_state=JumpMatrix(jump_matrix),
            holding_time=ContinuousMarkov(jump_rate),
            **kwargs,
        )

    @classmethod
    def from_partial_Q(cls, Q: numpy.ndarray, **kwargs):
        """
        Like `from_Q` but ignores diagonal elements of input array.
        """

        if not Q.ndim == 2:
            raise ValueError("Input array must be 2-d")
        if not numpy.alltrue(numpy.diff(Q.shape) == 0):
            raise ValueError("All dimensions of input must be of equal length")

        # Shortcut for square matrix
        diag_slice = slice(None, None, Q.shape[1] + 1)

        # Copy Q matrix and reset diagonal elements to 0
        jump_matrix = Q.copy()
        jump_matrix.flat[diag_slice] = 0
        if not numpy.alltrue(jump_matrix >= 0):
            raise ValueError(
                "Off-diagonal elements of input array must all be non-negative"
            )

        # Calculate jump_rates from off-diagonal elements of Q
        jump_rate = jump_matrix.sum(axis=-1)
        if not numpy.alltrue(jump_rate >= 0):
            raise ValueError(
                "Diagonal elements of input array must all be non-positive"
            )

        # Create jump matrix
        jump_matrix /= numpy.where(jump_rate == 0, 1, jump_rate)[:, None]
        jump_matrix.flat[diag_slice] = jump_rate == 0

        # Done
        return cls(
            next_state=JumpMatrix(jump_matrix),
            holding_time=ContinuousMarkov(jump_rate),
            **kwargs,
        )

    @property
    def jump_matrix(self) -> numpy.ndarray:
        return self.next_state.jump_matrix

    @property
    def jump_rate(self) -> numpy.ndarray:
        return self.holding_time.jump_rate

    def to_Q(self) -> numpy.ndarray:
        Q = self.jump_matrix.copy()
        numpy.fill_diagonal(Q, -1)
        Q *= self.jump_rate[:, None]
        Q[Q == 0] = 0  # Fix negative zero values

        return Q

    def to_P(self, timestep: float = 1) -> numpy.ndarray:
        return scipy.linalg.expm(self.to_Q() * timestep)

    def to_discrete(self, timestep: float = 1) -> DiscreteTimeMarkovChain:
        # Note that the "time" vector for the DiscreteTimeMarkovChain
        # is "steps" rather than "time"
        P = self.to_P(timestep)
        return DiscreteTimeMarkovChain.from_P(P, **self.get_params())
