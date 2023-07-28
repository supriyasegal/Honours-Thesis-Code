# Honours-Thesis-Code

## Folder: synthetic-study
This folder contains all of the code used in the synthetic study.

- `1-simulated-data-max-stable.R`: Run this code first, to generate the synthetic data and evaluate the max-stable approach.
- `2-kriging-synthetic-study.py`: Run this code second, to evaluate the kriging approach on the synthetic data, and create the box-and-whisker plots of the RMSE score.

## Folder: cross-validation
This folder contains all of the code used in leave-one-out cross-validation (LOOCV) on the Australian data.

- `kriging-cross-validation.py`: Run this code to perform LOOCV for the kriging approach.
- `max-stable-cross-validation.R`: Run this code to perform LOOCV for the max-stable approach.
- `data`
  - This folder contains extreme yearly rainfall data for the weather stations used in the study, for each climate region. It was extracted from NOAA's GHCN-D data set and transformed as outlined in the thesis.
  - NOAA's GHCN-D data set can be found here: https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/ 
