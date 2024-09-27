params.simulation_dir = "/n/fs/ragr-research/projects/fastppm/data/simulations"

params.nmutations = [100] //, 500, 1000, 2500]
params.nsamples   = [3]
params.coverage   = [30, 100, 1000]
params.seeds      = 1..20

process orchard {
    cpus 1
    memory '8 GB'
    time '59m'

    publishDir "nextflow_results/search/orchard/${id}", mode: 'copy', overwrite: true

    input:
        tuple val(id), path(variant_matrix), path(total_matrix)

    output:
        val(id)
        
    """
    echo $id
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
}
