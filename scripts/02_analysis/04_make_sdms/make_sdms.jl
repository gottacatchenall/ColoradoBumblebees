using DrWatson
@quickactivate :ColoradoBumblebees

const RADIUS = 400.0
const BUFFER = 20.0
const BIAS = 2.


const bias_one_species = [
    "Asclepias speciosa",           # Bias = 1
    "Eremogone fendleri",           # Bias = 1
    "Ipomopsis aggregata",          # Bias = 1
    "Taraxacum officinale",         # Bias = 1
    "Melilotus officinalis",        # Bias = 1
    "Lupinus argenteus",            # Bias = 1
    "Verbascum thapsus",            # Bias = 1
    "Geranium richardsonii",        # Bias = 1
    "Rosa woodsii",                 # Bias = 1
    "Mertensia lanceolata",         # Bias = 1
    "Frasera speciosa",             # Bias = 1
    "Oxytropis lambertii",          # Bias = 1
    "Thermopsis rhombifolia",       # bias = 1
    "Galium boreale",
   ]            


#    "Achillea millefolium",         # Bias = 0.7
#    "Galium boreale"    # Buffer = 40, Bias = 1

function main()
    occurrence_df = load_occurrence_data()
    species = sort(unique(occurrence_df.species))
    cluster = ColoradoBumblebees.CLUSTER

    job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
    this_species = species[job_id]

    # Saving using ArchGDAL sometimes writes scratch files to local dir during
    # the process of writing. Idk if its due to lack of RAM or something else,
    # but running each species in its own dir stops simulataneous overwriting of
    # scratch files
    mkpath("./tmp/$(this_species)")
    cd("./tmp/$(this_species)")
    @info "Job: $job_id, Species: $(this_species), WD: $(pwd())"

    b = this_species ∈ bias_one_species ? 1. : BIAS
    buffer = BUFFER 
    if this_species == "Achillea millefolium"
        b = 0.7
    end 
    if this_species == "Galium boreale"
        buffer = 40.
    end 

    sdms = make_sdms(
        this_species,
        occurrence_df; 
        cluster=cluster, 
        thickening_distance=RADIUS, 
        buffer=buffer,
        bias=b 
    )


    for (i, sdm) in enumerate(sdms)
        sdm_dir = sdmdir(sdm)
        mkpath(sdm_dir)
        ColoradoBumblebees.save(sdm)
    end
end

main()


