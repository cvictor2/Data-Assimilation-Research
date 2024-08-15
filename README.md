# 2D NSE Solver with Data Assimilation (Matlab)

This repository contains a 2D Navier-Stokes Equation (NSE) solver implemented in Matlab, along with various data assimilation methods. 
The solver is designed to simulate fluid flow in two dimensions and assimilate observational data to improve the accuracy and stability of the simulations.

## Features

- Efficient 2D NSE solver using pseudo-spectral methods
- Inbuilt data assimilation algorithms:
  - Azouani-Olson-Titi
  - Synchronization
- Customizable simulation settings

## Getting Started

### Prerequisites

To run the solver, you need to have Matlab installed on your system. You can download Matlab [here](https://www.mathworks.com/products/matlab.html).

### Installation

1. Clone the repository

2. Navigate to the project folder in Matlab.

3. Run scripts in the `code` folder to run simulations.



## Usage
The `mutual_nudging.m` is set up to run nudging intertwinements. The `mu_1` and `mu_2` variables affect the coupling between the two solutions.


The `self_synchronous.m` is set up to run synchronization intertwinements. The relevant variable to change for the coupling is `theta_1`, the `theta_2` value is set to be the Holder conjugate of `theta_1`.


Note that the scripts use an instance of the `DA_obs` class to store errors and relevant parameters. The `mu` variable in the `DA_obs` class has been replaced with the value of `mu_2` or `theta_2`.


You can configure the simulation parameters by modifying the parameters at the start of 'simulation_runner.m`

Every run requires at least 1 Settings and Parameters object, which dictate the general settings of the simulation, such as how often to plot and save variables, as well as the important parameters used in simulations, such as the kinematic viscosity coefficient.

The variables in the declaration of settings and parameters objects are instantiated in name-value pairs. That is, one needs to specify the parameter or setting name for the first argument, followed by the value for that variable, and this must be repeated for each value that is specified. Any variable not specified will be set to a default value.

Note that each data assimilation method will simulate it's own solution to NSE, so the more methods you test in a single run the longer the run will take.

### Notes
This code is custom made and features 4 different custom classes and I use an inputParser in each of their constructors. For simplicity all the tools and custom classes are stored in the helpers folder, which is loaded into the path at the beginning of the simulation_runner script.
The DA_obs class inherits from the handler class in order to make it pass by reference. This means that any function that any changes made to a DA_obs object are applied globally and not just inside the calling function. This will be used in a future update when the plotting functionality is adjusted.
The helpers folder contains useful functions that are utilized by the classes and the main program. This is necessary to run and is automatically added to the path at the beginning of the program.

