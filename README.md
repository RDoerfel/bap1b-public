# Brain Age Prediction in ADNI

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains code for validating multiple publicly available brain-age estimation packages on the Alzheimer's Disease Neuroimaging Initiative (ADNI) dataset. The analysis explores relationships between predicted brain age deviation (PAD) and cognitive decline, brain atrophy, and diagnosis status in the context of Alzheimer's disease. This repository accompanies the paper ["Prediction of Brain Age Using Structural Magnetic Resonance Imaging: A Comparison of Clinical Validity of Publicly Available Software Packages"](https://www.medrxiv.org/content/10.1101/2025.03.13.25323902v1).

## Project Structure

- `bap1b/`: Python package containing core functionality
  - `figures.py`: Plotting utilities for consistent visualization
- `scripts/`: Analysis and visualization scripts
  - `*.R`: Statistical analysis scripts in R
  - `*.py`: Data processing and visualization scripts in Python
  - `tab_*.py`: Scripts for generating formatted tables
  - `vis_*.py`: Scripts for generating figures
definitions
- `data/`: Data storage (not included in repository)
  - `megamastersheet_simulated.xlsx`: Main data file
  - `stats/`: Derived statistical data
- `results/`: Output figures and tables (not included in repository)

## Analysis Overview

The project conducts several analyses:

1. **Cross-Sectional Analysis**:
   - `cs_a.r`: Tests whether PAD at baseline can differentiate between clinical groups (CN, MCI, AD)
   - `cs_b.r`: Association between PAD and memory performance at baseline
   - `cs_summary.R`: Summary statistics for the cross-sectional analysis

2. **Longitudinal Analysis**:
   - `long_c.R`: Association between baseline PAD and conversion from CN to MCI/AD within four years
   - `long_de.R`: Association between baseline PAD and future decline in memory performance or atrophy
   - `long_fg.R`: Association between changes in PAD and changes in memory performance/atrophy

3. **Data Summary and Visualization**:
   - `data_summary.r`: Summary statistics across all analyses
   - Various `vis_*.py` scripts: Create publication-ready figures coresponding to the respective `*.R` analysis
   - Various `tab_*.py` scripts: Create formatted tables coresponding to the respective `*.R` analysis

## Installation

### Dependencies

This project uses both Python and R. The Python dependencies are managed through Poetry and include:
- pandas
- numpy
- matplotlib
- seaborn
- nibabel
- scikit-learn
- bids-validator
- clinica

R dependencies include:
- data.table
- dplyr
- LMMstar
- mets
- riskRegression
- boot
- mmrm
- writexl

### Setup

```bash
# Clone the repository
git clone https://github.com/RDoerfel/bap1b-public.git
cd bap1b-public

# For Python dependencies
pip install -e .  # Or use Poetry
poetry install

# For R packages
# Open R and run:
# install.packages(c("data.table", "dplyr", "writexl", "mmrm", "boot", "riskRegression", "mets", "LMMstar"))
```

### Running the Analysis
All files in scripts atart with an index, i.e., `01_` or `02_`. Those indicate the order in which one should execute the scripts. This will reproduce all tables and figures in the manuscript.  

## Data

This project uses the ADNI dataset, which is not included in this repository due to data usage agreements. To run the code, you need to:

1. Request access to the ADNI dataset: [ADNI Data Access](http://adni.loni.usc.edu/)
2. Process the data according to the methods described in our paper
3. Create a `megamastersheet_simulated.xlsx` file in the `data/` directory with the required structure. To run all scripts, we provide an example sheet that provides all necessary files with the required structure. Please note, the values are random samples drawn from distirbutions that match the original data, but this is not the data used for the publications. We were not allowed to publish the original data. The script `00_randomize_data.py` provides and example on how the randomization took place.

## Results

Results are stored in the `results/` directory but not tracked by git. Figures and tables include:
- Brain age deviation by diagnostic group
- Correlation between PAD and cognitive measures
- Longitudinal analysis of PAD changes
- ROC curves for prediction of conversion to MCI/AD

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{Doerfel2025,
  title = {Prediction of Brain Age Using Structural Magnetic Resonance Imaging: {{A}} Comparison of Clinical Validity of Publicly Available Software Packages},
  shorttitle = {Prediction of Brain Age Using Structural Magnetic Resonance Imaging},
  author = {D{\"o}rfel, Ruben P. and Ozenne, Brice and Ganz, Melanie and Svensson, Jonas and {Plav{\'e}n-Sigray}, Pontus},
  year = {2025},
  month = mar,
  publisher = {medRxiv},
  doi = {10.1101/2025.03.13.25323902},
  urldate = {2025-03-18}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- ADNI for providing the dataset
- Contributors to the brain age estimation packages evaluated in this study
