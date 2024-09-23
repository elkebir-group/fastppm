params.proj_dir           = "/n/fs/ragr-research/projects/fastBE/"
params.sim_script         = "${params.proj_dir}/scripts/simulation.py"
params.conversion_script  = "${params.proj_dir}/scripts/processing/make_projection_input.py" 
params.python_script      = "${params.proj_dir}/scripts/vafpp_lp.py"
params.fastbe_command     = "${params.proj_dir}/build/src/fastbe"
params.projection_command = "${params.proj_dir}/dependencies/projection/projection"

params.nmutations = [125, 250, 375, 500, 625, 750, 875, 1000]
params.nsamples   = [500]
params.seeds      = [0, 1, 2, 3, 4, 5]
params.coverage   = [50] 

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'

    input:
        tuple val(mutations), val(samples), val(coverage), val(seed)

    output:
        tuple file("sim_clonal_matrix.txt"), file("sim_mutation_to_clone_mapping.txt"), file("sim_obs_frequency_matrix.txt"), 
              file("sim_total_matrix.txt"), file("sim_tree.txt"), file("sim_usage_matrix.txt"), file("sim_variant_matrix.txt"),
              val("m${mutations}_n${mutations}_s${samples}_c${coverage}_r${seed}")

    """
    python '${params.sim_script}' --mutations ${mutations} --samples ${samples} --clones ${mutations} --coverage ${coverage} --seed ${seed} --output sim
    """
}

process regress_projection {
    cpus 1
    memory '2 GB'
    time '59m'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(id)

    output:
        tuple file("projection_output.txt"), file("timing.txt"), val(id)

    """
    python '${params.conversion_script}' ${clone_tree} ${freq_matrix} > input.txt
    /usr/bin/time -v '${params.projection_command}' input.txt projection_output.txt 1 2>> timing.txt
    """
}

process regress_fastbe {
    cpus 1
    memory '2 GB'
    time '59m'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(id)

    output:
        tuple file("cpp_results.json"), file("timing.txt"), val(id)

    """
    /usr/bin/time -v '${params.fastbe_command}' regress ${clone_tree} ${freq_matrix} --output cpp 2>> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))

    simulation =  parameter_channel | create_sim 
    fastbe_res = simulation | regress_fastbe
    projection_res = simulation | regress_projection

    fastbe_res | map { result, timing, name ->
      result.moveTo("nextflow_results/timings/${name}_fastbe_results.json")
      timing.moveTo("nextflow_results/timings/${name}_fastbe_timing.txt")
    }

    projection_res | map { result, timing, name ->
      result.moveTo("nextflow_results/timings/${name}_projection_results.txt")
      timing.moveTo("nextflow_results/timings/${name}_projection_timing.txt")
    }

}