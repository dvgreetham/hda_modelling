import typing

import attr
import numpy as np
import pandas as pd


def event_index(
    series: pd.Series, index: int = 0, default: typing.Any = np.nan
) -> typing.Any:
    try:
        return series[series].index[index]
    except IndexError:
        return default


@attr.s(auto_attribs=True)
class Demographics:
    cols: typing.Sequence[str] = attr.ib()
    subjid: str = "subjid"
    visdy: str = "visdy"

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract values for demographic features per subject.

        Since some demographic features can vary over time (e.g. employment
        status), the following filters are applied in order, assuming they match
        at least 1 row:

        * study = "ENR"
        * visit = "Baseline"
        * minimum number of missing values
        * earliest visdy

        Missing features are filled with NaN.

        Output is one row per subject.
        """
        # Check which columns are available
        present, missing = self.check_cols(df)

        # Select one row per subject
        out = df.groupby(self.subjid).apply(
            self.extract, present=present, missing=missing
        )

        # Revert datatypes
        for c in out:
            out[c] = out[c].astype(df[c].dtype)

        return out

    def check_cols(
        self, df: pd.DataFrame
    ) -> typing.Tuple[typing.Sequence[str], typing.Sequence[str]]:
        present = [c for c in self.cols if c in df]
        missing = [c for c in self.cols if c not in df]

        return present, missing

    # Given the augmented visits data, want to return the one row per subject with the demographic data
    def extract(self, df: pd.DataFrame, *, present=None, missing=None) -> pd.Series:
        if present is None or missing is None:
            present, missing = self.check_cols(df)

        # Reference to input
        out = df.copy()

        # If there's data for the ENR study, then limit attention to those
        if out.shape[0] > 1 and "studyid" in out and "ENR" in out["studyid"].values:
            out = out[out["studyid"] == "ENR"]

        # If there's data for the Baseline visits, then limit attention to those
        if out.shape[0] > 1 and "visit" in out and "Baseline" in out["visit"].values:
            out = out[out["visit"] == "Baseline"]

        # Choose rows with the minimum number of missing values
        if out.shape[0] > 1:
            missing_vals = out[present].isnull().sum(axis=1)
            out = out[missing_vals == missing_vals.min()]

        # Copy the first (by date)
        if out.shape[0] > 1:
            out = out.sort_values(self.visdy)
        out = out.iloc[0, :].copy()

        # Fill missing columns with NA
        if missing:
            out[missing] = np.nan

        # Return requested columns in order
        return out.loc[self.cols]

    def summary(self, demo: pd.DataFrame, max_cats: int = 10) -> None:
        for c in self.cols:
            if c in demo:
                counts = demo[c].value_counts(dropna=False)
                if counts.size <= max_cats:
                    print(f"{c.title()} frequencies:")
                    print(counts)
                else:
                    print(f"{c.title()} has {counts.size} unique values")
                print("")


@attr.s(auto_attribs=True)
class ProgressEvents:
    target: str = attr.ib()
    subjid: str = "subjid"
    visdy: str = "visdy"

    @property
    def increasing(self) -> bool:
        return self.target in {"stage"}

    def is_monotonic(self, x: pd.Series) -> bool:
        # Better to be explicit than rely on pandas.Series.is_monotonic, etc
        # So can be sure about strict vs non-strict monotonicity
        deltas = np.diff(x.values)
        if not self.increasing:
            deltas = -deltas
        return np.all(deltas >= 0)

    def transform(self, df: pd.DataFrame, stage_0, stage_1=None) -> pd.DataFrame:
        """
        Extract the intervals for the start and finish events of a given stage
        (or range of stages) per subject and derive interval for duration of stage.

        Stages of interest are between `stage_0` and `stage_1` inclusive.

        Output is one row per subject.
        """
        if stage_1 is None:
            stage_1 = stage_0

        out = df.groupby(self.subjid).apply(self.extract, stage_0, stage_1)

        return out

    def extract(
        self,
        data: typing.Union[pd.DataFrame, pd.Series],
        stage_0: int,
        stage_1: int = None,
    ) -> pd.Series:
        if stage_1 is None:
            stage_1 = stage_0

        out = pd.Series(
            {
                "fluct": None,
                "previous_visits": None,
                "required_visits": None,
                "following_visits": None,
                "start_min": None,
                "start_max": None,
                "finish_min": None,
                "finish_max": None,
                "duration_min": None,
                "duration_max": None,
            },
            dtype="object",
        )

        # TODO: If subject is deceased due to disease, then date_max should be the date of death.
        # Technically date_min can be inferred from age at day 0, but it's not helpful.
        date_min = np.NINF
        date_max = np.inf

        # Extract relevant time series
        if isinstance(data, pd.DataFrame):
            # This checks that visdy are all unique
            stages = data.set_index(self.visdy, verify_integrity=True)[self.target]
        else:
            stages = data
        stages = stages.dropna().sort_index()

        # Check for fluctuations
        if stages.size < 2:
            # Less than 2 visits
            out["fluct"] = 2
        else:
            out["fluct"] = 1 - self.is_monotonic(stages)

        # Take only subjects with monotonic stages
        if out["fluct"] > 0:
            return out

        # Find visits matching the required stages
        if self.increasing:
            previous = stages < stage_0
            required = (stage_0 <= stages) & (stages <= stage_1)
            following = stages > stage_1
        else:
            previous = stages > stage_0
            required = (stage_1 <= stages) & (stages <= stage_0)
            following = stages < stage_1

        # Count visits (for debugging)
        out["previous_visits"] = previous.sum()
        out["required_visits"] = required.sum()
        out["following_visits"] = following.sum()

        # Last day in earlier range
        out["start_min"] = event_index(previous, index=-1, default=date_min)

        # First day not in earlier range
        out["start_max"] = event_index(~previous, index=0, default=date_max)

        # Last day not in later range
        out["finish_min"] = event_index(~following, index=-1, default=date_min)

        # First day in later range
        out["finish_max"] = event_index(following, index=0, default=date_max)

        # Calculate durations
        out["duration_min"] = max(out["finish_min"] - out["start_max"], 0)
        out["duration_max"] = out["finish_max"] - out["start_min"]

        return out


@attr.s(auto_attribs=True)
class Fluctuations:
    target: str = attr.ib()
    subjid: str = "subjid"
    visdy: str = "visdy"

    @property
    def increasing(self) -> bool:
        return self.target in {"stage"}

    def transform(self, df: pd.DataFrame) -> pd.Series:
        return df.groupby(self.subjid).apply(self.extract).rename("fluct")

    def extract(self, data: typing.Union[pd.DataFrame, pd.Series]) -> int:
        if isinstance(data, pd.DataFrame):
            # Extract relevant time series
            #   This checks that visdy are all unique
            data = data.set_index(self.visdy, verify_integrity=True)[self.target]

        # Drop null values and sort by index
        data = data.dropna().sort_index()

        # Check for fluctuations
        if data.size < 2:
            # Not enough data to evaluate
            return 2
        else:
            # Better to be explicit than rely on pandas.Series.is_monotonic, etc
            # So can be sure about strict vs non-strict monotonicity
            deltas = np.diff(data.values)
            if not self.increasing:
                deltas = -deltas

            if np.all(deltas >= 0):
                # No fluctuations
                return 0
            else:
                # Fluctuations present
                return 1
