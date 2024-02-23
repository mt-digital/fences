##
# Simple fencing model to support policy paper
# 
# 
using Agents
using Memoize
using Iterators: filter


@enum Species insect bird herbivore carnivore livestock

# Do not list cattle 
default_birth_rates = 
    Dict(insect => 0.2, bird => 0.1, herbivore => 0.05, livestock => 0.05, carnivore => 0.02)

default_init_energy = 
    Dict(insect => 0.1, bird => 0.4, herbivore => 1, livestock => 2, carnivore => 3)

default_Δenergy = 
   Dict(insect => .01, bird => 0.04, herbivore => 0.1, livestock => 0.2, carnivore => 0.3)

prey_lookup = Dict(insect => [], bird => [insect], herbivore => [], 
                   livestock => [], carnivore => [herbivore])

@with_kw @agent Critter GridAgent{2} begin

    energy::Float64
    Δenergy::Float64

    # Birth rate can be a function of grass availability for, 
    # e.g., livestock, insects, birds.
    birth_rate::Union{Float64,Function}

    species::Species
    landowner::Union{Landholder,Nothing} = nothing

end


# Handle mortal dynamics (use birth rates for each agent for new births and 
# check for energy <= 0.0 for die-off). Grow vegetation depending on presence
# of pollinators (insects).
function model_step!(model)

    for pos in positions(model)
        # Grow vegetation depending on insects.
        
    end

end

#: Helper function for agent_step to find and eat available prey at agent position.
function eat_available_prey!(agent, model)
    if agent.species == carnivore
        prey_species = [herbivore, cattle]
    elseif agent.species == bird
        prey_species = [insect]
    else
        error("agent species, ${agent.species} has no known prey animals")
    end

    # find available prey
    available_prey = filter(a -> a.species ∈ prey_species,
                            agents_in_position(agent.pos, model))

    # take one of the available prey at random if any are available
    if !isempty(available_prey)
        prey = rand(available_prey)
        agent.energy += prey.energy
        remove_agent!(prey, model)
    end

end

# Handle movement, feeding. 
function agent_step!(agent, model)

    ## Move to a location depending on species and state...
    # ...livestock go to random location in landowner's holding...
    if agent.species == livestock
        new_pos = rand(agent.landowner.holding) 
        move_agent!(agent, new_pos, model)

    # ...birds and insects can random-walk anywhere (1 cell by default)...
    elseif agent.species ∈ [bird, insect]
        new_pos = rand(nearby_positions(agent.pos, model))

    # ...finally, herbivores and carnivores can't go to fenced locations.
    else
        new_pos = filter(pos -> check_pos_available(pos, model), 
                         nearby_positions(agent.pos, model))
    end
    move_agent!(agent, new_pos, model)

    # Feed, depending on species, location, and location's state...
    # ...if not a carnivore or bird, eat grass at location up to 4*Δenergy...
    if agent.species ∈ [insect, herbivore, livestock]
        grass_avail = grass[agent.pos]
        Δenergy = agent.Δenergy
        energy_to_eat = min(4*Δenergy, grass_avail)
        grass_avail -= energy_to_eat

    # ...if it is a carnivore, try to eat herbivore or livestock at location...
    else
        eat_available_prey!(agent, model)
    
    end

    # Lose Δenergy required to live for one time step.
    agent.energy -= agent.Δenergy

    # Die off if energy is negative. 
    if agent.energy ≤ 0
        remove_agent!(agent, model)
    end

end


@memoize function check_pos_available(pos, model)

   return pos ∈ model.available_cells

end


@with_kw mutable struct Landholder begin
    
    id::Int
    holding::Array{Tuple{Int,Int}}
    livestock::Array{Int} = []
    fence_probability::Function = 0.5
    fenced::Bool = false
end


function livestock_factory!(model, landholder; 
                           default_energy = 10, default_Δenergy = 0.1)

    # Select random point in holding to place cattle.
    pos = rand(landholder.holding)

    # Create and return new livestock Critter.
    return add_agent!(pos, model, default_energy, default_Δenergy, 
                      grass_frac -> 0.2 * grass_frac, livestock, landholder)

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
