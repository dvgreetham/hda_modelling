# hda_modelling
Code for the paper "Data-Driven Huntington’s Disease Progression Modelling and Estimation of Societal Cost in the UK"

## Project Organization

```raw
├── README.md          <- The top-level README for developers using this project.
├── analysis           <- Analysis scripts and notebooks. Naming convention is a number (for
│   │                     ordering) and a short description with `_` separators, e.g.
│   │                     `1_0_initial_data_exploration`. Add directories to organise and
|   │                     separate purposes.
│   ├── 0_prepare_data.ipynb <- Create enroll-all.feather which is needed to run the R analysis
|   │                     and to be able to reproduce work from the raw Enroll data
│   ├── ...
|
├── data
│   ├── external       <- Data from third party sources.
│   ├── interim        <- Intermediate data that has been transformed.
│   ├── processed      <- The final, canonical data sets for modelling.
│   └── raw            <- The original, immutable data dump.
│
├── src                <- Python source code for use in this project.
│   ├── hem
│   │   ├── __init__.py
│   │   ├── _version.py
│   │   ├── ...
│
├── src_R              <- R source code for use in this project.
|
```

## Prerequisites
 - Miniconda for Python 3.6 (or newer):
    Project uses conda for package management and virtual environments.

  - Conda version must be >= 4.6.
  - Use `conda init` to ensure that conda is correctly installed and configured.
  - Ensure that `conda activate` works even with no conda environment activated.

- R version 4.0 or later

- Rstudio

## Set-up

1. Obtain the code: either a git clone or an export
2. Download the [data]
3. In the virtual environment, install the additional R dependencies:

```bash
Rscript install.R
```

## Data
**Do not commit the original data files to source control.**


## Model fitting
In order to completely reproduce the fitted model,
the following analyses must be run in order:

1. [0_prepare_data.ipynb](analysis/0_prepare_data.ipynb)
2. [1_data_preparation/1_prepare_enroll_data](analysis/1_data_preparation/1_prepare_enroll_data)
3. [1_data_preparation/3_link_cost_of_care_to_enroll](analysis/1_data_preparation/3_link_cost_of_care_to_enroll)
4. [2_model_fitting/2_fit_CTHMM_to_enroll_training_set](analysis/2_model_fitting/2_fit_CTHMM_to_enroll_training_set)

N.b. more analyses would need to be run if the data changes
in order to reassess decisions,
such as the number of states.


## Style conventions

RMarkdown documents in this project
are written using [semantic linefeeds](https://rhodesmill.org/brandon/2012/one-sentence-per-line/).

--------

