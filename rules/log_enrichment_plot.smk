rule log_enrichment_plot:
    input:
        log_enrich_file = expand('{output_dir}/{log_enrich_file}', output_dir=output_dir, log_enrich_file=log_enrich_file)
    output:
        output_log_enrich_pdf = expand('{output_dir}/{log_enrich_plot_pdf}', output_dir=output_dir, log_enrich_plot_pdf=log_enrich_plot_pdf)
    conda:
        '../envs/r.yaml'
    resources:
        mem_gb = config['resource_specs']['intermediate_memory_gb'],
        time = config['runtime']['medium']
    script:
        '../scripts/log_enrichment_heatmap.R'