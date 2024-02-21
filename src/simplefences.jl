##
# Simple fencing model to support policy paper
# 
# 
using Agents


@enum Species insect bird herbivore carnivore livestock

# Do not list cattle 
default_birth_rates = 
    Dict(insect => 0.2, bird => 0.05, herbivore => 0.02, carnivore => 0.01)

    prey_lookup => Dict(insect => [], bird => [insect], herbivore => [], 
                        carnivore => [herbivore])

@agent Critter GridAgent{2} begin

    energy::Float64
    Δenergy::Float64

    # Birth rate can be a function of grass availability for, 
    # e.g., livestock, insects, birds.
    birth_rate::Union{Float64,Function}

    species::Species
end




@with_kw mutable struct Landholder begin
    
    id::Int
    holding::Array{Tuple{Int,Int}}
    livestock::Array{Int} = []
    fence_probability::Function = 0.5

end


function livestock_factory!(model, holding; 
                           default_energy = 10, default_Δenergy = 0.1)

    # Select random point in holding to place cattle.
    pos = rand(holding)

    # Create and return new livestock Critter.
    return add_agent!(pos, model, default_energy, default_Δenergy, 
                      grass_frac -> 0.2 * grass_frac, livestock)

end


function add_livestock!(model, landholder::Landholder)

    # Extract holding for convenience.
    holding = landholder.holding

    # Integer number of livestock based on parameter and holding size.
    n_livestock = 
        Int(round(model.init_livestock_per_area * sqrt(length(holding))))

    # Initialize livestock and assign to landholder. 
    landholder.livestock = 
        [livestock_factory!(model, holding) for _ in n_livestock]
end


function make_possible_holding(model)  #, available_cells)

    plot_size = 1 + rand(model.plot_size_distro)
    xmin, ymin = draw(model.available_cells)

    xmax = plot_location[1] + plot_size
    ymax = plot_location[2] + plot_size 

    possible_plot = [(x, y) for x in xmin:xmax for y in ymin:ymax]
end

# Helper function to make a valid holding.
function make_valid_holding(available_cells)

    space_width = model.space_width

    holding_empty = true
    holding = make_possible_holding(available_cells)
    while holding_empty

        holding = filter(coord -> 
                         (coord[1] ≤ space_width) && 
                         (coord[2] ≤ space_width) &&
                         (coord ∈ available_cells),
                        possible_holding)

        holding_empty = isempty(holding)

    end            

    delete!(available_cells, 
            findall(coord -> coord ∈ holding, available_cells))
end


function initialize_simplefences(;
        init_pops = Dict(insect => 20, bird => 10, herbivore => 10, 
                         carnivore => 5, livestock => 50, landholder => 5),
        space_width = 100,
        grass_regrowth_rate = 0.05,
        grass_consumption_rate = 0.1,  # same for livestock and herbivores 
        plot_size_distro_mean_min1 = 1,
        seed = 42
    )

    rng = MersenneTwister(seed)
    space = GridSpace(dims, periodic = true)

    # Model properties contain the grass as two arrays: whether it is fully grown
    # and the time to regrow. Also have static parameter `regrowth_time`.
    # Notice how the properties are a `NamedTuple` to ensure type stability.
    properties = (
        space_width = space_width,
        grass = ones(Float64, dims),
        grass_regrowth_rate = 0.1,
        landowners::Array{Landowner} = [],
        plot_size_distro = Poisson(plot_size_distro_mean_min1),
        available_cells = [(x, y) for x in 1:space_width for y in 1:space_width]
    )

    model = ABM(Critter, space;
        properties, rng, scheduler = Schedulers.randomly, warn = false
    )
        
    critters = initialize_critters!(model, init_pops)
    landholders = initialize_landholders!(model, init_pops.landholder)

    for p in positions(model)
        model.grass[p...] = 1.0
    end

end


function initialize_critters!(model, init_props)

    for species in instances(Species)
        for ii in 1:init_props[species]
            add_agent!(Critter, model, model.init_energy[species], 
                       model.init_Δenergy[species], 
                       model.init_birth_rate[species],
                       species)
        end
    end
end


function initialize_landholders!(model, n_landholders; n_livestock_coeff)

    for l_id in 1:n_landholders
        holding = make_valid_holding!(available_cells)
        new_landholder = Landholder(id = l_id, holding = holding)
        add_livestock!(new_landholder)
        push!(model.landholders, new_landholder)
    end

end
