##
# scripts/run_simplefences.jl to run the simplefences prototype model of 
# biodiversity change and fencing practices with the target environment of
# Kenya. A product of the Santa Fe Institute workshop "Biodiversity Protection
# – Forecasting Success and Reversing Challenges in a Complex
# Socio-economic-ecological World"
# 
# For more info visit https://www.santafe.edu/events/biodiversity-protection-forecasting-success-and-reversing-challenges-complex-socio-economic-ecological-world
#
# Author: Matthew A. Turner <maturner01@gmail.com>
# Date: 2024-04-30
#

# Use DrWatson to activate dependencies of reproducible environment.
using DrWatson
@quickactivate "fences"

# DataFrames contains innerjoin function used below.
using DataFrames

# Load model code.
include("../src/simplefences.jl")

# Run 
function run_default_simplefences(nsteps = 50; kwargs...)

  # Create a new simplefences model with default initialization.
  m = initialize_simplefences(; kwargs...)

  # Run model for fifty time steps. See src/simplefences.jl for definitions of
  # individual- and global-level stepping functions 
  # `agent_step!` and `model_step!`, and data collection 
  # specifications, `adata` and `mdata`.
  nsteps = 50
  adf, mdf = run!(m, agent_step!, model_step!, 50; adata, mdata)
  resdf = innerjoin(adf, mdf, on = :step)

  # Split off grass dataframe and transform into df with x and y coordinates.
  col_flat(x) = collect(Iterators.flatten(x))

  # This grid_size is the default in initialize_simplefences.
  grid_size = 100

  # x-coord indices repeat 1 to grid_size sequence.
  x = repeat(1:grid_size, grid_size*nsteps)

  # y-coord indices repeat next index grid_size times for each step, 
  # repeated over nsteps.
  y = col_flat(repeat([repeat(ii, grid_size) for ii in 1:grid_size], nsteps))
  
  # The step needs to be repeated over all coordinate combinations.
  new_step = col_flat([repeat([ii], grid_size*grid_size) for ii in 0:nsteps])

  grass_val = col_flat([col_flat(gmat) for gmat in resdf([:, :grass])])

  rename!(colname -> replace(colname, r"#41_s" => ""), resdf)

  return grass_val

#   return resdf
end
