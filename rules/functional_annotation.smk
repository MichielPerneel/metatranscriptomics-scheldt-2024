rule functional_annotation:
    input:
        pep=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.pep'),
    output:
        annotated=os.path.join(config['output_dir'], 'annotation', 'functional_eggnog', 'functional_annotation.emapper.annotations'),
    params:
        eggnog_ref_dir=config['eggnog_ref_dir'],
        eggnog_out_dir=directory(os.path.join(config['output_dir'], 'annotation', 'functional_eggnog'))
    log: os.path.join(config['log_dir'], 'functional_annotation.log')
    threads: 60
    shell:'''
    module load eggnog-mapper

    unset OMP_PROC_BIND

    mkdir -p {params.eggnog_out_dir}

    cd {params.eggnog_out_dir}

    emapper.py -m mmseqs --no_file_comments --cpu {threads} -i {input.pep} \
        -o functional_annotation --output_dir {params.eggnog_out_dir} \
        --data_dir {params.eggnog_ref_dir} --dbmem > {log} 2>&1
    '''