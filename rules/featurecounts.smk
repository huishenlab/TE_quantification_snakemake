rule featurecounts:
    input:
        bam = expand('{output_dir}/aligned_files/{{sample}}.bam', output_dir=output_dir),
        annotation_saf = f'{annotation_saf}'
    output:
        expand('{output_dir}/featureCounts/{{sample}}_fc.txt', output_dir=output_dir)
    params:
        featurecounts_params = config['featurecounts']['params']
    threads: config['resource_specs']['featureCounts_threads']
    conda:
        '../envs/featureCounts.yaml'
    resources:
        mem_gb = config['resource_specs']['intermediate_memory_gb'],
        time = config['runtime']['medium']
    shell:
        """
        featureCounts {params.featurecounts_params} {threads} -a {input.annotation_saf} -o {output} {input.bam} 
        """