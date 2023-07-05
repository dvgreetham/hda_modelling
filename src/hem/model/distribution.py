import functools
import typing

import numpy
import scipy.stats

DistributionType = scipy.stats.distributions.rv_frozen


def prepare_seed(seed=None):
    if seed is None:
        seed = numpy.random.SeedSequence()
        print(f"Entropy is {seed.entropy:_d}")
    elif isinstance(seed, int):
        seed = numpy.random.SeedSequence(seed)
    return seed


### Factories for generic distributions with caching
@functools.lru_cache(maxsize=100)
def delta_distribution(loc: typing.Union[float, int] = 0.0) -> DistributionType:
    # CDF: F(x) = I(x >= loc)
    if isinstance(loc, int):
        # return integer
        return scipy.stats.rv_discrete(values=((loc,), (1.0,)))
    else:
        # return float
        return scipy.stats.norm(loc=loc, scale=0)


@functools.lru_cache(maxsize=100)
def exponential_distribution(rate=1.0) -> DistributionType:
    if rate <= 0:
        return delta_distribution(numpy.inf)
    else:
        return scipy.stats.expon(scale=1 / rate)


@functools.lru_cache(maxsize=100)
def uniform_distribution(low=0.0, high=1.0) -> DistributionType:
    if low > high:
        raise ValueError
    elif low == high:
        return delta_distribution(loc=low)
    else:
        return scipy.stats.uniform(loc=low, scale=high - low)


@functools.lru_cache(maxsize=100)
def geometric_distribution(p=1.0) -> DistributionType:
    if p <= 0:
        return delta_distribution(numpy.inf)
    else:
        return scipy.stats.geom(p=p)


@functools.lru_cache(maxsize=100)
def categorical_distribution(pvals: typing.Tuple[float, ...]) -> DistributionType:
    # This assumes that `pvals` is already normalized
    # Note that `pvals` must be hashable in order for the cache to work
    return scipy.stats.rv_discrete(
        values=zip(*((v, p) for (v, p) in enumerate(pvals) if p > 0))
    )
