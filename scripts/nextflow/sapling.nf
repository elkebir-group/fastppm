params.simulation_dir       = "/n/fs/ragr-research/projects/fastppm-data/simulations"
params.sapling              = "/n/fs/ragr-research/projects/fastppm/dependencies/sapling-star/main.py"
params.create_sapling_input = "/n/fs/ragr-research/projects/fastppm/scripts/processing/make_sapling_input.py"
params.sapling_python       = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/sapling/bin/python"
params.python               = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"

params.nmutations = [500] // [50, 100, 250, 500]
params.nsamples   = [50, 100]
params.coverage   = [30, 100, 1000]
params.seeds      = 1..5
params.algorithm  = ["fastppm"]

process sapling {
    cpus 1
    memory '4 GB'
    time '12h'
    clusterOptions '--account=raphael'

    publishDir "nextflow_results/search/sapling-${algorithm}/${id}", mode: 'copy', overwrite: true

    input:
        tuple val(id), path(variant_matrix), path(total_matrix), val(algorithm)

    output:
        tuple val(id), path("sapling_output.txt"), path("timing.txt")
        
    """
    ${params.python} ${params.create_sapling_input} ${variant_matrix} ${total_matrix} > sapling_input.txt
    /usr/bin/time -v ${params.sapling_python} ${params.sapling} -f sapling_input.txt -o sapling_output.txt -a 1 -L ${algorithm} -t 1 2> timing.txt
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

    sim_files | combine(channel.fromList(params.algorithm)) | sapling
}
