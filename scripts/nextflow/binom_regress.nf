params.proj_dir           = "/Users/melkebir/Projects/fastppm"
params.simulation_dir     = "/Users/melkebir/Projects/fastppm-data/simulations"
params.output_dir         = "nextflow_results/binom_regress"
params.fastppm_command    = "${params.proj_dir}/build/src/fastppm-cli"
params.cvxopt_command     = "${params.proj_dir}/scripts/cvxopt_binom_regression.py"
params.gurobi_command     = "${params.proj_dir}/scripts/gurobi_lp_binom_regression.py"
params.time_command       = "gtime -v"

params.nmutations  = [100, 500, 1000, 2500]
params.nsamples    = [3]
params.seeds       = 1..20
params.coverage    = [30, 100, 1000]
params.nsegments   = [100, 1000, 10000]

process regress_binom_fastppm_binomial {
    cpus 1
    memory '2 GB'
    time '59m'

    input:
        tuple path(variant_matrix), path(total_matrix), path(clone_tree), val(id)

    output:
        tuple file("output.json"), file("timing.txt"), val(id)

    """
    ${params.time_command} '${params.fastppm_command}' -l binomial -v ${variant_matrix} -d ${total_matrix} -t ${clone_tree} -o output.json 2>> timing.txt
    """
}

process regress_binom_fastppm_binomial_K {
    cpus 1
    memory '2 GB'
    time '59m'

    input:
        tuple path(variant_matrix), path(total_matrix), path(clone_tree), val(segments), val(id)

    output:
        tuple file("output.json"), file("timing.txt"), val(id)

    """
    ${params.time_command} '${params.fastppm_command}' -K ${segments} -l binomial -v ${variant_matrix} -d ${total_matrix} -t ${clone_tree} -o output.json 2>> timing.txt
    """
}

process regress_binom_cvxopt {
    cpus 1
    memory '2 GB'
    time '59m'

    input:
        tuple path(variant_matrix), path(total_matrix), path(clone_tree), val(id)

    output:
        tuple file("output.json"), file("timing.txt"), val(id)

    """
    ${params.time_command} python '${params.cvxopt_command}' ${variant_matrix} ${total_matrix} ${clone_tree} > output.json 2>> timing.txt
    """
}

process regress_binom_gurobi {
    cpus 1
    memory '2 GB'
    time '59m'

    input:
        tuple path(variant_matrix), path(total_matrix), path(clone_tree), val(segments), val(id)

    output:
        tuple file("output.json"), file("timing.txt"), val(id)

    """
    ${params.time_command} python '${params.gurobi_command}' ${variant_matrix} ${total_matrix} ${clone_tree} ${segments} > output.json 2>> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))

    parameter_segments_channel = channel.fromList(params.nmutations)
                                        .combine(channel.fromList(params.nsamples))
                                        .combine(channel.fromList(params.coverage))
                                        .combine(channel.fromList(params.seeds))
                                        .combine(channel.fromList(params.nsegments))

    sim_files_segments = parameter_segments_channel | map { nmuts, nsamples, coverage, seed, segments ->
        segments_id             = "n${nmuts}_s${nsamples}_c${coverage}_r${seed}_k${segments}"
        segments_prefix         = "${params.simulation_dir}/n${nmuts}_s${nsamples}_c${coverage}_r${seed}"
        segments_freq_matrix    = "${segments_prefix}/sim_frequency_matrix.txt"
        segments_total_matrix   = "${segments_prefix}/sim_total_matrix.txt"
        segments_variant_matrix = "${segments_prefix}/sim_variant_matrix.txt"
        segments_tree           = "${segments_prefix}/sim_tree.txt"
        segments_usage_matrix   = "${segments_prefix}/sim_usage_matrix.txt"
        segments_variant_matrix = "${segments_prefix}/sim_variant_matrix.txt"
        [segments_variant_matrix, segments_total_matrix, segments_tree, segments, segments_id]
    }

    sim_files = parameter_channel | map { nmuts, nsamples, coverage, seed ->
        id             = "n${nmuts}_s${nsamples}_c${coverage}_r${seed}"
        prefix         = "${params.simulation_dir}/${id}/"
        freq_matrix    = "${prefix}/sim_frequency_matrix.txt"
        total_matrix   = "${prefix}/sim_total_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        tree           = "${prefix}/sim_tree.txt"
        usage_matrix   = "${prefix}/sim_usage_matrix.txt"
        variant_matrix = "${prefix}/sim_variant_matrix.txt"
        [variant_matrix, total_matrix, tree, id]
    }

    sim_files_segments_gurobi = sim_files_segments
    gurobi_res = sim_files_segments_gurobi | regress_binom_gurobi
    gurobi_res | map { result, timing, name ->
      result.moveTo("${params.output_dir}/${name}_gurobi_results.json")
      timing.moveTo("${params.output_dir}/${name}_gurobi_timing.txt")
    }

    sim_files_cvxopt = sim_files
    cvxopt_res = sim_files_cvxopt | regress_binom_cvxopt
    cvxopt_res | map { result, timing, name ->
      result.moveTo("${params.output_dir}/${name}_cvxopt_results.json")
      timing.moveTo("${params.output_dir}/${name}_cvxopt_timing.txt")
    }

    sim_files_fastppm_binomial = sim_files
    fastppm_binomial_res = sim_files_fastppm_binomial | regress_binom_fastppm_binomial
    fastppm_binomial_res | map { result, timing, name ->
      result.moveTo("${params.output_dir}/${name}_fastppm_results.json")
      timing.moveTo("${params.output_dir}/${name}_fastppm_timing.txt")
    }

    sim_files_segments_fastppm_binomial = sim_files_segments
    fastppm_binomial_res = sim_files_segments | regress_binom_fastppm_binomial_K
    fastppm_binomial_res | map { result, timing, name ->
      result.moveTo("${params.output_dir}/${name}_fastppm_results.json")
      timing.moveTo("${params.output_dir}/${name}_fastppm_timing.txt")
    }
}
