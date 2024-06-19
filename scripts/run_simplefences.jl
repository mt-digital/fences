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


using CSV
using DataFrames
using JSON
using Glob
using ProgressMeter
using RCall
using VideoIO

# Load model code.
include("../src/simplefences.jl")


# Run model with all defaults, kwargs passed to model init if desired.
function run_default_simplefences(nsteps = 50; 
                                  map_df_csvname = "map.csv",
                                  prevalence_csv_name = "prevalence.csv",
                                  model_to_restart = nothing,
                                  model_kwargs...)

    # If a model is passed to restart, use it as the model run...
    if !isnothing(model_to_restart)
        model = model_to_restart

    # ...otherwise initialize a new model.
    else
        model = initialize_simplefences(; model_kwargs...)
    end

    # Run model nsteps time steps. See src/simplefences.jl for definitions of
    # individual- and global-level stepping functions 
    # `agent_step!` and `model_step!`, and data collection 
    # specifications, `adata` and `mdata`.
    agent_df, model_df = run!(model, agent_step!, model_step!, nsteps; adata, mdata)
    result_df = innerjoin(agent_df, model_df, on = :step)

    # Reshape matrix entries in spatial result_df cells in long format.
    map_df = make_map_df(result_df, nsteps, model.space_width)

    # Prevalence df is just renamed and without the grass column.
    prevalence_df = rename(colname -> 
                           replace(
                               replace(colname, r"#\d*_s=" => ""), 
                               r"_species" => ""
                           ), 
                           select(result_df, Not([:grass])))

    if !isnothing(map_df_csvname)
        CSV.write(map_df_csvname, map_df)
    end

    # if !isnothing(prevalence_csv_name)
    #     CSV.write(prevalence_csv_name, prevalence_df)
    # end

    return prevalence_df, map_df, model
end


function run_make_frame(steps_per_frame, movie_tmp_dir = "movie/tmp", basename = ""; 
                        tstep = 0, prevalence_csv_name = "prevalence.csv", 
                        model_to_restart = nothing, prevalence_df = nothing, model_kwargs...) 

    # Run model for desired number of steps for the frame after save path set.
    csvname = joinpath(movie_tmp_dir, "csv", basename * string(tstep) * ".csv")

    new_prevalence_df, _, model_to_restart = 
        run_default_simplefences(
            steps_per_frame; 
            map_df_csvname = csvname,
            prevalence_csv_name,
            model_to_restart,
            model_kwargs...
        )

    if !isnothing(prevalence_df)
        new_prevalence_df.step .+= tstep
        prevalence_df = vcat(prevalence_df, new_prevalence_df)
    else
        prevalence_df = new_prevalence_df
    end

    CSV.write(prevalence_csv_name, prevalence_df)

    # Make CSVs of other data: animal and fence locations.
    fences_csv_name = joinpath(movie_tmp_dir, "csv", "fences-$basename$tstep.csv")
    animals_csv_name = joinpath(movie_tmp_dir, "csv", "animals-$basename$tstep.csv")

    fences_df, animals_df = make_overlay_data(model_to_restart)

    CSV.write(fences_csv_name, fences_df)
    CSV.write(animals_csv_name, animals_df)
    
    # Create the frame after R code sourced and frame path set.
    reval("source('scripts/plot.R');")
    frame_fname = joinpath(movie_tmp_dir, "png", basename * string(tstep) * ".png")

    # Create the map using R.
    reval("plot_map('$csvname', '$fences_csv_name', 
                    '$animals_csv_name', '$frame_fname');")

    return prevalence_df, model_to_restart

end


function make_overlay_data(model; species_to_plot = [herbivore, carnivore, livestock])

    fences_df = vcat([DataFrame(:landholder_id => landholder.id,
                                :x => [coord[1] for coord in landholder.holding],
                                :y => [coord[2] for coord in landholder.holding])
                      for landholder in filter(l -> l.fenced, model.landholders)]...)

    critters_to_plot = filter(a -> a.species in species_to_plot, 
                              collect(allagents(model)))

    animals_df = vcat([DataFrame(:a_id => critter.id, 
                                    :x => critter.pos[1], 
                                    :y => critter.pos[2], 
                                    :species => critter.species)
                             for critter in critters_to_plot]...)

    return fences_df, animals_df

end


function make_map_df(result_df, maxstep = 50, grid_size = 100)

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
end


function fences_movie(max_steps = 50, 
                      steps_per_frame = 2; 
                      tmp_dir = "movie/tmp",
                      prevalence_csv_name = "prevalence.csv",
                      prevalence_figure_name = "figures/prevalence_dynamics.pdf",
                      movie_file = "movie/default.mp4",
                      framerate = 4, 
                      model_kwargs...)

    # Clear previous tmp directories.
    for tmpfile in glob(joinpath(tmp_dir, "*", "*"))
        rm(tmpfile)
    end

    # Initialize model and make first frame.
    prevalence_df, model_to_restart = run_make_frame(steps_per_frame; model_kwargs...)

    # Initialize tstep at 2 from previous call then build frames.
    tstep = steps_per_frame
    while tstep < max_steps

        prevalence_df, model_to_restart = 
            run_make_frame(steps_per_frame; 
                           tstep, model_to_restart, 
                           prevalence_df,
                           model_kwargs...)

        tstep += steps_per_frame
    end

    # Plot species prevalence dynamics using R before making movie.
    reval("plot_prevalences('$prevalence_csv_name', '$prevalence_figure_name')")

    ## Create movie from frames, following 
    ## https://juliaio.github.io/VideoIO.jl/stable/writing/#Iterative-Encoding
    frame_dir = joinpath(tmp_dir, "png")
    # First load file names...
    frame_names = readdir(frame_dir)
    # ...then sort by index.
    frame_idxs = map(x -> split(x, ".")[1], frame_names)
    p = sortperm(parse.(Int, frame_idxs))
    frame_names = frame_names[p]

    println()
    println("Made it here!")
    println()

    # Set up movie encoder.
    encoder_options = (crf=23, preset="medium")

    # Write movie.
    first_frame = load(joinpath(frame_dir, frame_names[1]))
    open_video_out(movie_file, first_frame, framerate = framerate, 
                   encoder_options = encoder_options, 
                   target_pix_fmt = VideoIO.AV_PIX_FMT_YUV420P) do writer

        @showprogress "Encoding video frames..." for i in eachindex(frame_names)
            frame = load(joinpath(frame_dir, frame_names[i]))
            write(writer, frame)
        end
    end
end
