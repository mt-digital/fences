##
# Simple fencing model to support policy paper
# 
# 
using Agents
using DataStructures
using Distributions
using Memoize
using Parameters
using Random


@enum Species insect bird herbivore carnivore livestock landholder


# Do not list cattle 
prey_lookup = Dict(insect => [], bird => [insect], herbivore => [], 
                   livestock => [], carnivore => [herbivore])


@with_kw mutable struct Landholder
    
    id::Int
    holding::Array{NTuple{2,Int}} = []
    livestock::Array{AbstractAgent} = []
    fence_probability::Float64 = 0.05
    fenced::Bool = false
end


@with_kw mutable struct Critter <: AbstractAgent

    id::Int 
    pos::NTuple{2,Int}
    energy::Float64
    Δenergy::Float64

    # Birth rate can be a function of grass availability for, 
    # e.g., livestock, insects, birds.
    birth_rate::Float64

    species::Species
    landholder::Union{Landholder,Nothing} = nothing

end


# Handle mortal dynamics (use birth rates for each agent for new births and 
# check for energy <= 0.0 for die-off). Grow vegetation depending on presence
# of pollinators (insects).
function model_step!(model)

    # Reset total grass count.
    model.total_grass = 0.0

    # Iterate through all grass, grow more based on insects, calc new total.
    for pos in positions(model)

        # Grow vegetation depending on insects.
        c = model.insect_grass_coeff

        insects_at_pos = count(!isnothing, Iterators.filter(a -> a.species == insect,
                                           agents_in_position(pos, model)))

        curr_grass = model.grass[pos...]

        # Grass grows at least c, but more insects means more pollinators and greater growth.
        model.grass[pos...] = min(curr_grass + (c * (1 + insects_at_pos)), 1)

        # Add current patch position grass to total.
        model.total_grass += model.grass[pos...]
    end

    # Possibly build new fences with certain probability.
    for landholder in filter(l -> !l.fenced, model.landholders)
        if rand() < landholder.fence_probability
            landholder.fenced = true
        end
    end
end

#: Helper function for agent_step to find and eat available prey at agent position.
function eat_available_prey!(agent, model)
    if agent.species == carnivore 
        prey_species = [herbivore, livestock]
    elseif agent.species == bird
        prey_species = [insect]
    else
        error("agent species, $(agent.species) has no known prey animals")
    end

    # find available prey
    # TODO make it prey within a certain larger radius, e.g. within a radius
    # of five cells, or whatever number of cells, using Agents.nearby_agents()
    # see https://juliadynamics.github.io/Agents.jl/stable/api/#Nearby-Agents
    # (in window on other screen).
    available_prey = Iterators.filter(a -> a.species ∈ prey_species,
                                      agents_in_position(agent.pos, model))

    # take one of the available prey at random if any are available
    if !isempty(available_prey)
        prey = rand(collect(available_prey))
        agent.energy += prey.energy
        remove_agent!(prey, model)
    end

end


# Handle movement, feeding. 
function agent_step!(agent, model)

    ## Move to a location depending on species and state...
    # ...livestock go to random location in landholder's holding...
    if agent.species == livestock
        try
            new_pos = rand(agent.landholder.holding) 
            move_agent!(agent, new_pos, model)
        catch
            println(agent)
            error("stop")
        end

    # ...birds and insects can random-walk anywhere (1 cell by default)...
    elseif agent.species ∈ [bird, insect]

        new_pos = rand(collect(nearby_positions(agent.pos, model)))

    # ...finally, herbivores and carnivores can't go to fenced locations.
    else
        available_locations = collect(Iterators.filter(pos -> check_pos_available(pos, model), 
                                                       nearby_positions(agent.pos, model)))
        new_pos = agent.pos
        if !isempty(available_locations)
            new_pos = rand(available_locations)
        end
    end
    move_agent!(agent, new_pos, model)

    # Feed, depending on species, location, and location's state...
    # ...if not a carnivore or bird, eat grass at location up to 4*Δenergy...
    if agent.species ∈ [insect, herbivore, livestock]
        grass_avail = model.grass[agent.pos...]
        Δenergy = agent.Δenergy
        energy_to_eat = min(20*Δenergy, grass_avail)
        model.grass[agent.pos...] -= energy_to_eat

    # ...if it is a carnivore, try to eat herbivore or livestock at location...
    else
        eat_available_prey!(agent, model)
    end

    # Maybe reproduce.
    if rand() < agent.birth_rate

        new_id = model.max_id + 1
        model.max_id = new_id
        species = agent.species
        new_agent = Critter(new_id, agent.pos, model.init_energy[species],
                            model.Δenergy[species], model.birth_rate[species],
                            species, agent.landholder)

        add_agent!(new_agent, model)
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


function livestock_factory!(model, landholder; 
                           default_energy = 10, default_Δenergy = 0.1)

    # Select random point in holding to place cattle.
    a_idx = length(model.agents) + 1
    pos = rand(landholder.holding)
    new_critter = Critter(a_idx, pos, model.init_energy[livestock],
                          model.Δenergy[livestock], model.birth_rate[livestock],
                          livestock, landholder)

    # Create and return new livestock Critter.
    return add_agent!(new_critter, model)

end


function add_livestock!(model, landholder::Landholder)

    # Extract holding for convenience.
    holding = landholder.holding

    # Integer number of livestock based on parameter and holding size.
    n_livestock = 
        Int(round(model.init_livestock_per_area * sqrt(length(holding))))

    # Initialize livestock and assign to landholder. 
    landholder.livestock = 
        [livestock_factory!(model, landholder) for _ in 1:n_livestock]
end


function make_possible_holding(model)  #, available_cells)

    plot_size = 1 + rand(model.plot_size_distro)
    xmin, ymin = rand(model.available_cells)

    space_width = model.space_width
    while (xmin < (0.1 * space_width)) || (ymin < (0.1 * space_width))
        xmin, ymin = rand(model.available_cells)
    end

    xmax = xmin + plot_size
    ymax = ymin + plot_size 

    possible_plot = [(x, y) for x in xmin:xmax for y in ymin:ymax]
end


# Helper function to make a valid holding.
function make_valid_holding!(model)

    space_width = model.space_width

    # Create possible holdings and make valid until a non-empty holding is created.
    holding = []
    while isempty(holding)

        holding = make_possible_holding(model)

        holding = filter(coord -> 
                         (coord[1] ≤ space_width - (space_width * .1)) && 
                         (coord[2] ≤ space_width - (space_width * .1)) &&
                         (coord ∈ model.available_cells),
                        holding)

        deleteat!(model.available_cells, 
                  findall(coord -> coord ∈ holding, model.available_cells))

    end            

    return holding
end


function initialize_simplefences(;
        init_pops = Dict(insect => 1000, bird => 300, herbivore => 50, 
                         carnivore => 25, livestock => 100, landholder => 20),

        birth_rate = 
            Dict(insect => 0.02, bird => 0.005, herbivore => 0.01, 
                 livestock => 0.008, carnivore => 0.008),

        init_energy = 
            Dict(insect => 0.8, bird => 1.0, herbivore => 5, 
                 livestock => 5, carnivore => 5),

        Δenergy = 
           Dict(insect => .002, bird => 0.01, herbivore => 0.01, 
                livestock => 0.02, carnivore => 0.2),

        space_width = 500,
        grass_regrowth_rate = 0.1,
        grass_consumption_rate = 0.2,  # same for livestock and herbivores 
        plot_size_distro_mean_min1 = 15,
        insect_grass_coeff = 0.02,
        init_livestock_per_area = 1,
        seed = 42,
        initial_fenced_prob = 0.5
    )

    dims = (space_width, space_width)
    rng = MersenneTwister(seed)
    space = GridSpace(dims, periodic = false, metric = :manhattan)

    # Model properties contain the grass as two arrays: whether it is fully grown
    # and the time to regrow. Also have static parameter `regrowth_time`.
    # Notice how the properties are a `NamedTuple` to ensure type stability.
    properties = Dict(
        :space_width => space_width,
        :grass => ones(Float64, dims),
        :grass_regrowth_rate => 0.1,
        :landholders => [],
        :plot_size_distro => Poisson(plot_size_distro_mean_min1),
        :available_cells => [(x, y) for x in 1:space_width for y in 1:space_width],
        # Coefficient relating number of insects to vegetation growth (i.e., due
        # to presence of pollinators.
        :insect_grass_coeff => insect_grass_coeff,
        :init_livestock_per_area => init_livestock_per_area,
        :init_pops => init_pops,
        :birth_rate => birth_rate,
        :init_energy => init_energy,
        :Δenergy => Δenergy,
        :total_fenced_area => 0.0,
        :initial_fenced_prob => initial_fenced_prob,
        :total_grass => 0.0,
        :max_id => 0
    )

    model = ABM(Critter, space;
        properties, rng, scheduler = Schedulers.randomly, warn = false
    )
        
    critters = initialize_critters!(model)
    println(length(filter(a -> a.species == herbivore, collect(allagents(model)))))
    landholders = initialize_landholders!(model, init_pops[landholder])
    update_total_fenced_area!(model)

    for p in positions(model)
        model.grass[p...] = 1.0
    end

    model.max_id = maximum(map(a -> a.id, allagents(model)))

    return model
end


function update_total_fenced_area!(model)

    # Holding area equals length of coords; include in total if holding is fenced.
    holding_sizes = [length(l.holding) for l in model.landholders if l.fenced]
    model.total_fenced_area = isempty(holding_sizes) ? 0.0 : sum(holding_sizes)

end


function initialize_critters!(model)

    critter_idx = 1
    println("Going to create $(model.init_pops[herbivore]) herbivores!")
    for species in instances(Species)
        # For each species, initialize the init population size specified in init_pops.
        if species ∉ [landholder, livestock]
            for ii in 1:model.init_pops[species]
                pos = random_position(model)
                new_agent = Critter(critter_idx, pos, model.init_energy[species], 
                                    model.Δenergy[species], model.birth_rate[species],
                                    species, nothing)

                add_agent!(new_agent, model)
                critter_idx += 1
            end
        end
    end
end


function initialize_landholders!(model, n_landholders; n_livestock_coeff = 1)

    for l_id in 1:n_landholders
        holding = make_valid_holding!(model)
        fenced = false

        if rand() < model.initial_fenced_prob
            fenced = true
        end

        new_landholder = Landholder(id = l_id, holding = holding, fenced = fenced)
        add_livestock!(model, new_landholder)
        push!(model.landholders, new_landholder)
    end
    
end


frac_a(v) = sum(v .== a) / length(v)

is_minority(x) = x.group == 1
frac_a_ifdata(v) = isempty(v) ? 0.0 : frac_a(collect(v))

adata = [(:species, v -> sum(v .== s)) for s in instances(Species)]
# adata = [(:species, v -> length(collect(v)), s -> isequal(s,insect))]

mdata = [:total_grass, :total_fenced_area, :grass]
