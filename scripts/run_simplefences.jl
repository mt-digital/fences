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

    landscape_df = make_landscape_df(resdf, nsteps, 100)

    # Prevalence df is just renamed and without the grass column.
    prevalence_df = rename(colname -> 
                           replace(
                           replace(colname, r"#\d*_s=" => ""), 
                               r"_species" => ""
                           ), 
                           select(resdf, Not([:grass])))

    return prevalence_df, landscape_df

end


function make_landscape_df(result_df, maxstep = 50, grid_size = 100)

    # Split off grass dataframe and transform into df with x and y coordinates.
    col_flat(x) = collect(Iterators.flatten(x))

    # Number of steps is one more than max step, index starts at step 0.
    nsteps = maxstep + 1

    # x-coord indices repeat 1 to grid_size sequence.
    x = repeat(1:grid_size, grid_size*nsteps)

    # y-coord indices repeat next index grid_size times for each step, 
    # repeated over nsteps.
    y = col_flat(repeat([repeat([ii], grid_size) for ii in 1:grid_size], nsteps))

    # The step needs to be repeated over all coordinate combinations.
    new_step = col_flat([repeat([ii], grid_size*grid_size) for ii in 0:maxstep])

    # Now flatten each grass matrix from each step in the result_df.
    grass_flat = col_flat([col_flat(gmat) for gmat in result_df[:, :grass]])

    # TODO create flattened fence location matrix just like grass;
    # TODO the matrix will contain 0's for no fence, 1 where there's fence.

    return DataFrame(:step => new_step, :x => x, :y => y, 
                     :grass_layer => grass_flat)
    # TODO
    # return DataFrame(:step => new_step, :x => x, :y => y, 
    #                  :grass_layer => grass_vec, :fence_layer = fence_vec)
end


# Create fence perimeter data from output fenced area data. 

