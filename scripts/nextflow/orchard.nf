params.simulation_dir       = "/n/fs/ragr-research/projects/fastppm/data/simulations"
params.orchard_star         = "/n/fs/ragr-research/projects/fastppm/dependencies/orchard-star/bin/orchard"
params.orchard              = "/n/fs/ragr-research/projects/fastppm/dependencies/orchard/bin/orchard"
params.create_orchard_input = "/n/fs/ragr-research/projects/fastppm/scripts/processing/make_orchard_input.py"
params.parse_orchard_output = "/n/fs/ragr-research/projects/fastppm/scripts/processing/parse_orchard_output.py"
params.python               = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"

params.nmutations = [100]//, 500, 1000, 2500]
params.nsamples   = [3]
params.coverage   = [30, 100, 1000]
params.seeds      = 1..20

process orchard_star {
    cpus 16
    memory '8 GB'
    time '59m'

    publishDir "nextflow_results/search/orchard-star-l2/${id}", mode: 'copy', overwrite: true

    input:
        tuple val(id), path(variant_matrix), path(total_matrix)

    output:
        tuple val(id), path("results.npz"), path("tree.txt"), path("timing.txt"), path("orchard_mutations.ssm"), path("orchard_params.json")
        
    """
    ${params.python} ${params.create_orchard_input} ${variant_matrix} ${total_matrix} -o orchard
    /usr/bin/time -v ${params.python} ${params.orchard_star} orchard_mutations.ssm orchard_params.json results.npz -k 5 -f 50 -n 16 -p 2> timing.txt
    ${params.python} ${params.parse_orchard_output} results.npz --output tree.txt
    """
}

process orchard {
    cpus 16
    memory '8 GB'
    time '59m'

    publishDir "nextflow_results/search/orchard/${id}", mode: 'copy', overwrite: true

    input:
        tuple val(id), path(variant_matrix), path(total_matrix)

    output:
        tuple val(id), path("results.npz"), path("tree.txt"), path("timing.txt")
        
    """
    python ${params.create_orchard_input} ${variant_matrix} ${total_matrix} -o orchard
    /usr/bin/time -v python ${params.orchard} orchard_mutations.ssm orchard_params.json results.npz -k 5 -f 50 -n 16 -p 2> timing.txt
    python ${params.parse_orchard_output} results.npz --output tree.txt
    """
}


workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))

    sim_files = parameter_channel | map { nmuts, nsamples, coverage, seed ->
        id = "n${nmuts}_s${nsamples}_c${coverage}_r${seed}"
        prefix = "${params.simulation_dir}/${id}/"
        freq_matrix    = "${prefix}/sim_frequency_matrix.txt"
        total_matrix   = "${prefix}/sim_total_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        tree           = "${prefix}/sim_tree.txt"
        usage_matrix   = "${prefix}/sim_usage_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        [id, variant_matrix, total_matrix]
    }

    sim_files | orchard
    sim_files | orchard_star
}
