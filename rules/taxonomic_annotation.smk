rule create_eukprot_db:
    """
    This rule turns the EukProt sequences into a mmseqs DB
    """
    input: config['eukprot_ref_fa']
    output: os.path.join(config['eukprot_ref_dir'], 'eukprot')
    shell:'''
    unset OMP_PROC_BIND

    module load MMseqs2

    mmseqs createdb {input} {output}
    '''

rule create_metatranscriptome_db:
    """
    This rule turns the metatranscriptome sequences into a mmseqs DB
    """
    input: os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.pep')
    output: os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome_mmseqsDB')
    shell:'''
    unset OMP_PROC_BIND

    module load MMseqs2

    mmseqs createdb {input} {output}
    '''

rule eukprot_annotation:
    input:
        metatranscriptome_mmseqsDB=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome_mmseqsDB'),
        eukprot_mmseqsDB=os.path.join(config['eukprot_ref_dir'], 'eukprot'),
    output:
        os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukprot', 'eukprot_annotation.m8')
    params:
        eukprot_out_dir=os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukprot'),
        tmp=config['scratch_dir']
    threads: 60
    resources:
        mem_mb=250000
    log:
        os.path.join(config['log_dir'], 'eukprot_annotation.log')
    shell:'''
    unset OMP_PROC_BIND

    module load MMseqs2

    # Query assembly against reference
    mmseqs search {input.metatranscriptome_mmseqsDB} {input.eukprot_mmseqsDB} {params.eukprot_out_dir}/resultDB {params.tmp} -s 6  > {log} 2>&1

    # Extract first hit
    mmseqs filterdb {params.eukprot_out_dir}/resultDB {params.eukprot_out_dir}/resultDB.firsthit --extract-lines 1 >> {log} 2>&1
    mmseqs convertalis {input.metatranscriptome_mmseqsDB} {input.eukprot_mmseqsDB} {params.eukprot_out_dir}/resultDB.firsthit {output} >> {log} 2>&1
    '''

rule eukulele_eukprot:
    input:
        pep = os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.pep')
    output:
        eukulele_done = os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukulele', 'EUKulele_done.txt')
    params:
        sampledir = os.path.join(config['output_dir'], 'assembly', 'protein'),
        eukulele_directory = os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukulele'),
        eukulele_reference_dir = config['eukulele_eukprot_ref_dir'],
        reference_pep_fa_gz = config['eukulele_eukprot_ref_fa'],
    conda: 'EUKulele_2'
    log: os.path.join(config['log_dir'], 'eukulele_eukprot.log')
    shell:'''
    which EUKulele
    EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --reference_dir {params.eukulele_reference_dir} --ref_fasta {params.reference_pep_fa_gz} -o {params.eukulele_directory} > {log} 2>&1
    touch {output.eukulele_done}
    '''
