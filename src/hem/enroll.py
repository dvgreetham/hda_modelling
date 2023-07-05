"""
Read CSV tables from the Enroll-HD dataset.

Three kinds of tables are supported: profiles, participation and study visits.
Filename for study visits can be "enroll.csv" or "registry.csv".
"""

import pathlib
import typing

import numpy as np
import pandas as pd

from .data import fix_numerical_column, restore_or_create_backup

FilenameType = typing.Union[pathlib.PurePath, str, None]

region_codes = {
    "Europe": 0,
    "Northern America": 1,
    "Australasia": 2,
    "Latin America": 2,
}


class Enroll:
    def __init__(self, *, raw: bool = False, backup: bool = False) -> None:
        self.raw = raw
        self.backup = backup

        self.missing = {
            "Unknown": {"numeric": 9999, "string": "UNKNOWN"},
            "Missing": {"numeric": 9998, "string": "MISSING"},
            "N/A": {"numeric": 9997, "string": "NOTAPPL"},
            "Wrong": {"numeric": 9996, "string": "WRONG"},
        }
        self.codes = {
            "region": ["Europe", "Northern America", "Australasia", "Latin America"],
            "sex": {"m": "male", "f": "female"},
            "jobclas": {
                1: "full-time",
                2: "part-time",
                3: "self-employed",
                4: "not working",
            },
            "res": {1: "rural", 2: "village", 3: "town", 4: "city"},
            "maristat": {
                1: "single",
                2: "married",
                3: "partnership",
                4: "divorced",
                5: "widowed",
                6: "separated",
            },
            "race": {
                1: "Causasian",
                2: "American Black",
                3: "Hispanic",
                8: "Native American",
                16: "Asian",
                15: "Mixed",
                6: "Other",
            },
            "handed": {1: "right", 2: "left", 3: "mixed"},
        }
        self.replace = {
            "age": {"<18": "18"},
            "age_0": {"<18": "18"},
            "caghigh": {">70": "75"},
            "caglow": {">28": "32"},
        }

        self.subjects_ = None
        self.studies_ = None
        self.visits_ = None

    def read_data(
        self,
        data_path: typing.Union[pathlib.Path, str],
        visit_filename: typing.Union[FilenameType, typing.Collection[FilenameType]],
        profile_filename: FilenameType = "profile.csv",
        participation_filename: FilenameType = "participation.csv",
        *,
        backup: bool = False,
        raw: bool = False,
    ) -> pd.DataFrame:
        """
        Read and process data files for analysis.

        Visits file name can be "enroll.csv" or "registry.csv".
        Assumes that "profile.csv" and "participation.csv" are in the same
        directory as the visits data file.
        """
        if not isinstance(data_path, pathlib.Path):
            data_path = pathlib.Path(data_path)
        data_path = data_path.resolve()

        # Read visits data
        if isinstance(visit_filename, (pathlib.PurePath, str)):
            visit_filename = [visit_filename]

        for filename in visit_filename:
            if filename:
                self.load_visits(data_path / filename)

        # Read profile data
        if profile_filename:
            self.load_subjects(data_path / profile_filename)

        # Read participation data
        if participation_filename:
            self.load_studies(data_path / participation_filename)

        # Return visits data
        return self.get_visits()

    def read_csv(self, filename: pathlib.Path) -> pd.DataFrame:
        df = pd.read_csv(filename, sep="\t", low_memory=False)
        print(f"Read {filename}: {df.shape}")
        return df

    def prepare_data(self, df: pd.DataFrame, index_cols) -> pd.DataFrame:
        # Fix numerical columns
        if not self.raw:
            for name, replace in self.replace.items():
                if name in df:
                    if self.backup:
                        restore_or_create_backup(df, name)
                    df[name] = fix_numerical_column(df[name], replace)

        # Replace missing data
        if not self.raw:
            df = self.mask_missing_values(df)

        # Convert to categoricals
        if not self.raw:
            for name, codes in self.codes.items():
                if name in df and not hasattr(df[name], "cat"):
                    if self.backup:
                        restore_or_create_backup(df, name)
                    if isinstance(codes, typing.Mapping):
                        codes = pd.Series(codes, dtype="category")
                        codes.index = codes.index.astype(df[name].dtype)
                        df[name] = df[name].map(codes, na_action="ignore")
                    else:
                        df[name] = pd.Categorical(
                            df[name], categories=list(codes), ordered=False
                        )

        # Check for duplicated rows and missing values
        index = df[index_cols]
        if index.duplicated().any():
            raise ValueError("Duplicated rows found in index")
        if index.isnull().any(axis=None):
            raise ValueError("Index columns contain missing values")

        # Done
        return df

    def load_visits(self, filename: pathlib.Path) -> None:
        """
        Read per-visit data file.

        Visits data file can be "enroll.csv" or "registry.csv".

        Contains one row per subject per study per visit.
        """
        df = self.read_csv(filename)

        # Prepare data
        df = self.prepare_data(df, ["subjid", "studyid", "seq"])

        # Merge with previous data
        if self.visits_ is None:
            self.visits_ = df
        else:
            self.visits_ = self.visits_.merge(df, how="outer")

    def load_subjects(self, filename: pathlib.Path) -> None:
        """
        Read per-subject data file.

        Profile data file is "profile.csv".

        Contains one row per subject.
        """
        df = self.read_csv(filename)

        # Prepare data
        df = self.prepare_data(df, ["subjid"])

        # Merge with previous data
        if self.subjects_ is None:
            self.subjects_ = df
        else:
            self.subjects_ = self.subjects_.merge(df, how="outer")

    def load_studies(self, filename: pathlib.Path) -> None:
        """
        Read per-study data file.

        Profile data file is "participation.csv".

        Contains one row per subject per study.
        """
        df = self.read_csv(filename)

        # Prepare data
        df = self.prepare_data(df, ["subjid", "studyid"])

        # Merge with previous data
        if self.studies_ is None:
            self.studies_ = df
        else:
            self.studies_ = self.studies_.merge(df, how="outer")

    def get_subjects(
        self,
        subject_cols: typing.Sequence[str] = None,
    ) -> pd.DataFrame:
        df = self._extract_cols(
            self.subjects_, index_cols=["subjid"], value_cols=subject_cols
        )
        if df.empty:
            raise ValueError("No subjects data")

        # Prepare data
        df = self.prepare_data(df, ["subjid"]).sort_values(["subjid"])

        return df

    def get_studies(
        self,
        study_cols: typing.Sequence[str] = None,
        subject_cols: typing.Sequence[str] = None,
    ) -> pd.DataFrame:
        df = self._extract_cols(
            self.studies_, index_cols=["subjid", "studyid"], value_cols=study_cols
        )
        if df.empty:
            raise ValueError("No studies data")

        df = self._merge_cols(
            df, self.subjects_, index_cols=["subjid"], value_cols=subject_cols
        )

        # Prepare data
        df = self.prepare_data(df, ["subjid", "studyid"]).sort_values(
            ["subjid", "studyid"]
        )

        return df

    def get_visits(
        self,
        visit_cols: typing.Sequence[str] = None,
        subject_cols: typing.Sequence[str] = None,
        study_cols: typing.Sequence[str] = None,
    ) -> pd.DataFrame:
        df = self._extract_cols(
            self.visits_, index_cols=["subjid", "studyid", "seq"], value_cols=visit_cols
        )
        if df.empty:
            raise ValueError("No visits data")

        df = self._merge_cols(
            df, self.subjects_, index_cols=["subjid"], value_cols=subject_cols
        )
        df = self._merge_cols(
            df, self.studies_, index_cols=["subjid", "studyid"], value_cols=study_cols
        )

        # Prepare data
        df = self.prepare_data(df, ["subjid", "studyid", "seq"]).sort_values(
            ["subjid", "studyid", "seq"]
        )

        return df

    def _extract_cols(
        self,
        df: pd.DataFrame,
        index_cols: typing.Sequence[str] = None,
        value_cols: typing.Sequence[str] = None,
    ) -> pd.DataFrame:
        if df is None:
            # Nothing to select
            return pd.DataFrame()
        elif value_cols is not None:
            # Only select specified index and value columns
            cols = []
            if index_cols is not None:
                cols.extend(index_cols)
            cols.extend(value_cols)
            return df[cols]
        else:
            # Select all columns
            return df

    def _merge_cols(
        self,
        lhs,
        rhs,
        index_cols: typing.Sequence[str] = None,
        value_cols: typing.Sequence[str] = None,
        *,
        how: str = "left",
        **kwargs,
    ) -> pd.DataFrame:
        data = self._extract_cols(rhs, index_cols=index_cols, value_cols=value_cols)
        if lhs is None:
            return data
        elif data.empty:
            return lhs
        else:
            return lhs.merge(data, on=index_cols, how=how, **kwargs)

    def missing_value_summary(
        self, df: pd.DataFrame, axis=0, present=True
    ) -> pd.DataFrame:
        """
        Data frame containing summary of missing values per column.
        """
        status = df.apply(self.missing_status)
        labels = list(self.missing.keys())
        labels.append("Null")
        if present:
            labels.append("Present")

        return pd.DataFrame(
            {label: (status == label).sum(axis=axis) for label in labels}
        )

    def missing_status(
        self, series: pd.Series, present_value: str = "Present", na_value: str = "Null"
    ) -> pd.Series:
        """
        Data frame containing summary of missing values per column.
        """

        data_type = "numeric" if pd.api.types.is_numeric_dtype(series) else "string"
        mapping = {codes[data_type]: label for label, codes in self.missing.items()}
        return series.map(
            lambda x: na_value if pd.isnull(x) else mapping.get(x, present_value)
        )

    def mask_missing_values(
        self, df: pd.DataFrame, na_value: typing.Any = np.nan, *, inplace: bool = False
    ) -> typing.Optional[pd.DataFrame]:
        """
        Replace missing data codes with numpy.nan. Default is in-place.
        """
        if not inplace:
            df = df.copy()

        status = df.apply(self.missing_status)
        df[status != "Present"] = na_value

        return None if inplace else df


def convert_region(df: pd.DataFrame, col="region"):
    """
    Convert "region" to numerical codes.

    Uses `region_codes` dictionary to convert values to codes.
    """
    return df[col].map(region_codes, na_action="ignore").astype("Int64")
