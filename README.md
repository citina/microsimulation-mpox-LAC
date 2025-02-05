# microsimulation-mpox-LAC

## Description
This repository hosts the MATLAB code for the microsimulation model used in the study titled **_"Viral introductions and return to baseline sexual behaviors maintain low-level mpox incidence in Los Angeles County, USA, 2023-2024."_** The model was employed to generate the results and figures for the associated journal article. The code is compatible with MATLAB R2024a.

## Installation

### Prerequisites
Ensure you have MATLAB R2024a installed on your machine to run the simulation code. 

### Setup
1. Clone this repository to your local machine: `git clone https://github.com/citina/microsimulation-mpox-LAC.git`
2. Navigate to the cloned directory: `cd microsimulation-mpox-LAC`
   
## Structure and Usage
The repository is structured as follows:
- **mpox_simulation_scripts/**: Contains the main simulation scripts and helper functions.
- **input/**: Includes the Excel spreadsheet `Inputs_mpox2024_set2.xlsx`, which contains model parameters and transition probabilities used in the model.

### Running Simulations
To run a base case scenario:
1. Open MATLAB.
2. Navigate to the `mpox_simulation_scripts/` directory.
3. Execute the desired scenario by running one of the following shell scripts:
   * `Mpox2024_ShellScript_sq1.m`: For scenarios with a probability of isolation set at 0.2 and varying Force of Infection (FoI) values of 0.7, 1.45, 2.2.
   * `Mpox2024_ShellScript_sq2.m`: For scenarios with varying probabilities of isolation by adjusting `piso_input` and a fixed FoI of 2.2.
   * `Mpox2024_ShellScript_sq3.m`: For scenarios with a fixed FoI of 0.7 and adjustments to `X_input` to alter the constant number of disease importations.
   * `Mpox2024_ShellScript_sq4.m`: For scenarios targeting specific periods with lower transmission rates.

## Contact
For queries or collaboration requests, please contact Citina Liang at [citinal@usc.edu](mailto:citinal@usc.edu). Please reach out before using the model to ensure it is applied appropriately.


