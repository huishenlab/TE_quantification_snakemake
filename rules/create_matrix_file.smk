rule create_matrix_file:
    input:
        samples_fc = expand('{output_dir}/featureCounts/{sample}_fc.txt', sample=SAMPLES, output_dir=output_dir),
        intronic_intergenic_tes = f'{intronic_intergenic_tes}'
    output:
        expand('{output_dir}/{count_matrix}', output_dir=output_dir, count_matrix=count_matrix_files)
    params:
        sample_names = expand('{sample}', sample=SAMPLES),
        output_dir = f'{output_dir}',
        single_cell = config['samples']['single_cell'],
        qc = config['samples']['qc'],
        threads = config['resource_specs']['fread_threads']
    conda:
        '../envs/r.yaml'
    resources:
        mem_gb = config['resource_specs']['intermediate_memory_gb'],
        time = config['runtime']['medium']
    script:
        '../scripts/create_matrix_file.R'