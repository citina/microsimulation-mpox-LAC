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
- **Code/**: Contains the main simulation scripts and helper functions.
- **Data/**: Includes an Excel spreadsheet `Inputs_mpox2024_set2.xlsx` which contains model parameters.
- **Scenarios/**: Contains shell scripts for running different simulation scenarios.

### Running Simulations
To run a base case scenario:
1. Open MATLAB.
2. Navigate to the `Scenarios/` directory.
3. Run `Mpox2024_ShellScript_sq1.m` script: `run('mpoxshell.m')`

## Contact
For queries or collaboration requests, please contact Citina Liang at [citinal@usc.edu](mailto:citinal@usc.edu). Please reach out before using the model to ensure it is applied appropriately.


