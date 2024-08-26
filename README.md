# Fences

Agent-based model and analysis code to simulate the causes and effects of different choices of land ownership, cattle production, and fencing decisions on biodiversity and trophic networks, meant to represent threats to biodiversity in Kenya. This is currently just a prototype. Model code is located in [src/simplefences.jl](src/simplefences.jl) and analysis code is in [scripts/plot.R](scripts/plot.R).

## Get the code and install dependencies

To get this code, clone the repository at the command line
```sh
git clone https://github.com/mt-digital/fences.git
```
or otherwise (e.g. GUI in RStudio, Visual Studio, GitHub Desktop, etc.).

You must have [Julia](https://julialang.org/downloads/) and [R](https://cloud.r-project.org/) installed to run the model and analysis, respectively.

To install Julia dependencies, open a Julia console and run the following to first install [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) globally, then activate the DrWatson project and install dependencies via `Pkg.instantiate()`:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

Currently R dependencies must be installed (quasi-)manually. Personally I would open [scripts/plot.R](scripts/plot.R) in RStudio, expecting that RStudio will ask if I want to install any dependencies that are `require`d but not currently installed. I would click yes, and let RStudio run the installation. 

## The model

The model implements a simple trophic network in a landscape where vegetation grows and is consumed by cattle and herbivores, which are in turn consumed by carnivores. Landowners own cattle and build fences to protect their cattle. At this prototype stage, trophic "food web" interactions are _ad hoc_, with predation probabilities and birth/death rates set primarily to produce minimally-interesting dynamics. Similarly, landowner . Furthermore, in Kenya there are land conservancies who incentivize landowners to allow wildlife in their landholdings, which exposes landowners not in conservancies to wildlife that threaten either cattle or crops. When conservancies encroach and non-conservancy neighbors build fences, a non-conservancy landowner may be more likely to build a fence that may disrupt trophic networks. These are important factors that would be included were this project to move beyond prototype stage.

The model is implemented in [src/simplefences.jl](src/simplefences.jl) using the [Agents.jl](https://juliadynamics.github.io/Agents.jl/stable/) library in the [Julia programming language](https://julialang.org/learning/).

### Quickstart to run the model

First load model code
```julia
# Use DrWatson to activate dependencies of reproducible environment.
using DrWatson
@quickactivate "fences"

# DataFrames contains innerjoin function used below.
using DataFrames

# Load model code.
include("src/simplefences.jl")

# Create a new simplefences model with default initialization.
m = initialize_simplefences()

# Run model for fifty time steps. See src/simplefences.jl for definitions of
# individual- and global-level stepping functions `agent_step!` and `model_step!`,
# and data collection specifications, `adata` and `mdata`.
nsteps = 50
adf, mdf = run!(m, agent_step!, model_step!, 50; adata, mdata)
resdf = innerjoin(adf, mdf, on = :step)

# Print joined results dataframe for inspection.
println(resdf)
```
The printed dataframe should look something like this:
<img width="1880" alt="image" src="https://github.com/mt-digital/fences/assets/2425472/6524ae10-2d43-4ae0-9a9d-bf95719f6f5a">

Continue reading for how to analyze the time series data contained in this dataframe, and how to make a tile plot (heatmap) of vegetation coverage.


## Analysis

We use R and the [`tidyverse`](https://www.tidyverse.org/) ecosystem for analyzing model output. 
Note a full-length book by Hadley Wickham, Mine Çetinkaya-Rundel, and Garrett
Grolemund on this analytical approach is available online, 
[https://r4ds.hadley.nz](R for Data Science).


### .mp4 to .gif

Use the following script to convert an mp4 to a .gif:


```bash
src="movie/default.mp4"
dest="output.gif"
palette="movie/tmp/palette.png"

ffmpeg -i $src -vf palettegen -y $palette
ffmpeg -i $src -i $palette -lavfi paletteuse -y $dest
```
