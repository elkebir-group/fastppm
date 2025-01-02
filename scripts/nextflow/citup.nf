params.proj_dir       = "/n/fs/ragr-research/projects/fastppm/"
params.simulation_dir = "/n/fs/ragr-research/projects/fastppm-data/simulations"
params.python         = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python"

params.citup_python = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/citupenv/bin/python"
params.citup_iter   = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/citupenv/bin/run_citup_iter.py"
params.citup_input_script = "${params.proj_dir}/scripts/processing/make_citup_input.py"
params.citup_parse_output = "${params.proj_dir}/scripts/processing/parse_citup_output.py"

params.citup_star = "${params.proj_dir}/dependencies/citup-star/citup/run_citup_iter.py"

params.nmutations = [3, 5, 10]
params.nsamples   = [3, 5]
params.coverage   = [30, 100, 1000]
params.seeds      = 1..10            // 20

process citup_star {
    cpus 16
    memory '8 GB'
    time '24h'

    publishDir "nextflow_results/search/citup-star/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(freq_matrix), val(clones), val(id)

    output:
        tuple path("results.hd5"), path("tree.txt"), path("timing.txt"), val(id)
        
    """
    export PATH=/n/fs/ragr-research/projects/fastppm/dependencies/citup-star/src:$PATH
    /usr/bin/time -v ${params.python} ${params.citup_star} ${freq_matrix} results.hd5 \
                  --loglevel INFO --maxjobs 16 --min_nodes ${clones} --max_nodes ${clones} --submit local 2> timing.txt
    '${params.python}' '${params.citup_parse_output}' results.hd5 > tree.txt
    """
}

process citup {
    cpus 16
    memory '8 GB'
    time '24h'
    errorStrategy 'ignore'

    publishDir "nextflow_results/search/citup/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(freq_matrix), val(clones), val(id)

    output:
        tuple path("results.hd5"), path("tree.txt"), path("timing.txt"), val(id)

    """
    export PATH=/n/fs/ragr-data/users/schmidt/miniconda3/envs/citupenv/bin:$PATH
    /usr/bin/time -v '${params.citup_iter}' ${freq_matrix} results.hd5 --loglevel INFO --maxjobs 16 --min_nodes ${clones} --max_nodes ${clones} --submit local 2> timing.txt
    '${params.citup_python}' '${params.citup_parse_output}' results.hd5 > tree.txt
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
        freq_matrix_T  = "${prefix}/sim_frequency_matrix_transpose.txt"
        total_matrix   = "${prefix}/sim_total_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        tree           = "${prefix}/sim_tree.txt"
        usage_matrix   = "${prefix}/sim_usage_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        [freq_matrix_T, nmuts, id]
    }

    // sim_files | citup
    sim_files | citup_star
}
