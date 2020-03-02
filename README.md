# Inferring transmission trees to quantify spatiotemporal spread of visceral leishmaniasis

This repository contains MCMC code and simulation code for estimating the parameters of, and simulating, the individual-level spatial kernel visceral leishmaniasis transmission model described in ['Inferring transmission trees to guide targeting of interventions against visceral leishmaniasis and post-kala-azar dermal leishmaniasis'](https://doi.org/10.1101/2020.02.24.20023325) [1].

## Prerequisites

* MATLAB R2017b+, available from <https://uk.mathworks.com/downloads/>, requires a users licence. Installation and activation instructions: <https://www.mathworks.com/help/install/index.html>.

* The following MATLAB toolboxes are required to run the MCMC code:
  * Optimization Toolbox v8.0+
  * Statistics and Machine Learning Toolbox v11.2+

* Julia v1.0.5+, available from <https://julialang.org/downloads/>. Installation instructions: <https://julialang.org/downloads/platform/>.

* The following Julia packages are required to run the simulation code:
  * CSV v0.5.20+
  * DataFrames v0.19.4+
  * DelimitedFiles
  * Distributions v0.21.11+
  * FileIO v1.2.0+
  * JLD v0.9.2+
  * JLD2 v0.1.11+
  * LinearAlgebra
  * Plots v0.20.3+
  * SparseArrays
  * StatsBase v0.32.0+


## Data

The MCMC code requires data (files 'data_final.mat' and 'data_final2.mat') which cannot be made publicly available as they contain personally identifiable information. If you would like to obtain a copy of the data please contact <lloyd.chapman@lshtm.ac.uk> in the first instance. Simulated data and code for testing the MCMC  algorithm on it will be added to this repository in due course. All data necessary to run the Julia simulation code is contained within this repository.

## Installing

Clone/download this project into a folder on your machine using the green button at the top right of this page.

### MCMC code

Once MATLAB has been installed and you have data in the required format, the MCMC code can be run by opening MATLAB, changing the working directory to the folder in which the files are saved, setting the parameters for the model(s) you would like to test in the SpecifyModel.m script, and then running each model by entering

```matlab
>> RunMCMC(i)
```

at the command prompt, where i is the number of the model in SpecifyModel.m you wish to run. Note that the code takes ~2.5 weeks to run 100,000 iterations on a 3.2GHz 64GB RAM 32 CPU Dell Precision Tower 7910.

The MCMC output can then be processed to compare different models tested, derive estimates for the model parameters, calculate the contribution of different infection states, reconstruct transmission trees, and estimate numbers of secondary infections/cases by typing

```matlab
>> RunPostProcessing
```

at the command prompt.

### Julia code

Once Julia has been installed, the simulation code can be run by opening a Command Prompt (Windows)/Terminal window (Mac/Linux) and typing `julia` at the prompt, changing the working directory to where the code is saved, then installing the required packages by entering the following commands:

```julia
julia> using Pkg
julia> Pkg.add(["CSV","DataFrames","DelimitedFiles","Distributions","FileIO","JLD","JLD2","LinearAlgebra","Plots","SparseArrays","StatsBase"])
```

then running 

```julia
julia> include("RunSims.jl") 
```

with the values of the simulation parameters, and number of posterior samples and simulations per sample set in RunSims.jl (currently 100 samples maximum). Note that the code takes ~16 hours to run 10,000 simulations on a desktop with the above specs.

The output of the simulations can then be processed and plotted by running

```julia
julia> include("PlotSims.jl")
```

## Built With

* MCMC code: [MATLAB 9.3.0.713579 (R2017b)](https://uk.mathworks.com/downloads/)
* Simulation code: [Julia 1.0.5](https://julialang.org/downloads/)

## Authors

* **Lloyd Chapman:** <lloyd.chapman@lshtm.ac.uk>

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments

* The techniques used to accelerate the adaptation and convergence of the MCMC algorithm are due to Simon Spencer. Further details can be found in the  supporting information for the paper and have been submitted for publication in a separate paper [2].
* The code for plotting the shaded arrows in the transmission trees uses the MATLAB functions "quiver_thick.m" and "arrow_thick.m" from [Open Earth Tools](https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/plot_fun/).

## References
1. Chapman LAC, Spencer SEF, Pollington TM, Jewell CP, Mondal D, Alvar J, Hollingsworth TD, Cameron MM, Bern C, Medley GF. Inferring transmission trees to guide targeting of interventions against visceral leishmaniasis and post-kala-azar dermal leishmaniasis. medRxiv 2020; doi: [10.1101/2020.02.24.20023325](https://doi.org/10.1101/2020.02.24.20023325)
 

2. SEF Spencer. Accelerating adaptation in the adaptive Metropolis Hastings random walk algorithm. (Submitted), 1–19 (2020).
