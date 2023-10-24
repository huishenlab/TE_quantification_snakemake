rule log_enrichment:
    input:
        tes_cpm = expand('{output_dir}/{tes_cpm_file}', output_dir=output_dir, tes_cpm_file=tes_cpm_file)
    output:
        output_log_enrich = expand('{output_dir}/{log_enrich_file}', output_dir=output_dir, log_enrich_file=log_enrich_file)
    params:
        output_dir = f'{output_dir}',
        sample_name = config['samples']['sample_name'],
        qc = config['samples']['qc']
    conda:
        '../envs/r.yaml'
    resources:
        mem_gb = config['resource_specs']['intermediate_memory_gb'],
        time = config['runtime']['medium']
    script:
        '../scripts/log_enrichment.R'