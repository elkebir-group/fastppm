params.proj_dir           = "/n/fs/ragr-research/projects/fastppm"
params.simulation_dir     = "/n/fs/ragr-research/projects/fastppm-data/simulations"
params.output_dir         = "${params.proj_dir}/nextflow_results/regress"

params.make_projection_input = "${params.proj_dir}/scripts/processing/make_projection_input.py" 

params.projection_command = "${params.proj_dir}/dependencies/projection/projection"
params.fastppm_command    = "${params.proj_dir}/build/src/fastppm-cli"
params.cvxopt_command     = "${params.proj_dir}/scripts/cvxopt_binom_regression.py"
params.gurobi_command     = "${params.proj_dir}/scripts/gurobi_lp_binom_regression.py"
params.time_command       = "/usr/bin/time -v"

params.nmutations  = [100, 500, 1000, 2500]
params.nsamples    = [500]
params.seeds       = 1..20
params.coverage    = [30, 100, 1000]
params.nsegments   = [50, 100]//, 1000] //, 10000]

process regress_l2_projection {
    cpus 1
    memory '1 GB'
    time '59m'

    scratch true
    publishDir "${params.output_dir}/projection_l2/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(clone_tree), path(freq_matrix), path(weight_matrix), val(id)

    output:
        tuple path("output.txt"), path("timing.txt"), val(id)

    """
    python '${params.make_projection_input}' ${clone_tree} ${freq_matrix} --weight_matrix ${weight_matrix} > input.txt
    /usr/bin/time -v '${params.projection_command}' input.txt output.txt 1 2>> timing.txt
    """
}

process regress_l2_fastppm {
    cpus 1
    memory '1 GB'
    time '59m'

    scratch true
    publishDir "${params.output_dir}/fastppm_l2/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(clone_tree), path(variant_matrix), path(total_matrix), path(weight_matrix), val(id)

    output:
        tuple path("output.json"), path("timing.txt"), val(id)

    """
    ${params.time_command} '${params.fastppm_command}' -t ${clone_tree} -v ${variant_matrix} -d ${total_matrix} -w ${weight_matrix} \
                           --output output.json -l l2 2>> timing.txt
    """
}

process regress_binom_fastppm_binomial {
    cpus 1
    memory '32 GB'
    time '4h'

    scratch true
    publishDir "${params.output_dir}/fastppm_binomial/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(clone_tree), path(variant_matrix), path(total_matrix), val(id)

    output:
        tuple path("output.json"), path("timing.txt"), val(id)

    """
    ${params.time_command} '${params.fastppm_command}' -l binomial -v ${variant_matrix} -d ${total_matrix} -t ${clone_tree} -o output.json 2>> timing.txt
    """
}

process regress_binom_fastppm_binomial_K {
    cpus 1
    memory '32 GB'
    time '4h'

    scratch true
    publishDir "${params.output_dir}/fastppm_binomial_K/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(clone_tree), path(variant_matrix), path(total_matrix), val(segments), val(id)

    output:
        tuple path("output.json"), path("timing.txt"), val(id)

    """
    ${params.time_command} '${params.fastppm_command}' -K ${segments} -l binomial_K -v ${variant_matrix} -d ${total_matrix} -t ${clone_tree} -o output.json 2>> timing.txt
    """
}

process regress_binom_cvxopt {
    cpus 1
    memory '32 GB'
    time '4h'

    scratch true
    publishDir "${params.output_dir}/cvxopt_binomial/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(clone_tree), path(variant_matrix), path(total_matrix), val(id)

    output:
        tuple path("output.json"), path("timing.txt"), val(id)

    """
    ${params.time_command} python '${params.cvxopt_command}' ${variant_matrix} ${total_matrix} ${clone_tree} > output.json 2>> timing.txt
    """
}

process regress_binom_gurobi {
    cpus 1
    memory '32 GB'
    time '4h'

    scratch true
    publishDir "${params.output_dir}/gurobi_binomial_K/${id}/", mode: 'copy', overwrite: true

    input:
        tuple path(clone_tree), path(variant_matrix), path(total_matrix), val(segments), val(id)

    output:
        tuple path("output.json"), path("timing.txt"), val(id)

    """
    module load gurobi
    ${params.time_command} python '${params.gurobi_command}' ${variant_matrix} ${total_matrix} ${clone_tree} ${segments} > output.json 2>> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))
                               .combine(channel.fromList(params.nsegments))

    /* Load simulated data. */
    simulations = parameter_channel | map { nmuts, nsamples, coverage, seed, nsegments ->
        id             = "n${nmuts}_s${nsamples}_c${coverage}_r${seed}"
        prefix         = "${params.simulation_dir}/n${nmuts}_s${nsamples}_c${coverage}_r${seed}"
        freq_matrix    = "${prefix}/sim_frequency_matrix.txt"
        total_matrix   = "${prefix}/sim_total_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        tree           = "${prefix}/sim_tree.txt"
        usage_matrix   = "${prefix}/sim_usage_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        weight_matrix  = "${prefix}/sim_weight_matrix.txt"
        [tree, variant_matrix, total_matrix, freq_matrix, weight_matrix, nsegments, id]
    }

    /* Select required files and run methods. */
    // simulations | map { [it[0], it[1], it[2], it[5], "${it[6]}_k${it[5]}"] } | regress_binom_fastppm_binomial_K
    // simulations | map { [it[0], it[1], it[2], it[5], "${it[6]}_k${it[5]}"] } | regress_binom_gurobi
    // simulations | map { [it[0], it[1], it[2], it[6]] } | unique | regress_binom_cvxopt
    // simulations | map { [it[0], it[1], it[2], it[6]] } | unique | regress_binom_fastppm_binomial
    // simulations | map { [it[0], it[3], it[4], it[6]] } | unique | regress_l2_projection 
    simulations | map { [it[0], it[1], it[2], it[4], it[6]] } | unique | regress_l2_fastppm 
}
