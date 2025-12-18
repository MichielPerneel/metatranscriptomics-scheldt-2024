rule create_pfam_mmseqs_db:
    """
    This rule converts the Pfam-A file into an mmseqs searchable profile database.
    """
    input: config['pfam']
    output: os.path.join(config['pfam_dir'], 'pfam_mmseqs_db')
    log: os.path.join(config['log_dir'], 'protein', 'create_pfam_mmseqs_db.log')
    params:
        msa_db = os.path.join(config['pfam_dir'], 'pfam_msa_db'),
        tmp = os.path.join(config['scratch_dir'], 'pfam_db_tmp')
    shell:'''
    unset OMP_PROC_BIND
    module load MMseqs2
    mmseqs convertmsa {input} {params.msa_db} > {log} 2>&1
    mmseqs msa2profile {params.msa_db} {output} --match-mode 1 > {log} 2>&1
    mmseqs createindex {output} {params.tmp} -k 5 -s 7 > {log} 2>&1
    '''

rule predict_long_orfs:
    """
    Predict long ORFs from the metatranscriptome using TransDecoder.LongOrfs.
    """
    input:
        metatranscriptome=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.fasta')
    output:
        longORFs_pep=os.path.join(config['output_dir'], 'assembly', 'protein', 'final_metatranscriptome.fasta.transdecoder_dir', 'longest_orfs.pep')
    params:
        protein_dir = os.path.join(config['output_dir'], 'assembly', 'protein')
    log:
        os.path.join(config['log_dir'], 'protein', 'predict_long_orfs.log')
    shell:'''
        unset OMP_PROC_BIND
        module load TransDecoder
        module load Perl-bundle-CPAN/5.36.1-GCCcore-12.3.0

        cd {params.protein_dir}

        TransDecoder.LongOrfs -t {input.metatranscriptome} > {log} 2>&1

        touch {output.longORFs_pep}
    '''

rule create_longorfs_mmseqs_db:
    input:
        longORFs_pep=os.path.join(config['output_dir'], 'assembly', 'protein', 'final_metatranscriptome.fasta.transdecoder_dir', 'longest_orfs.pep')
    output:
        longorfs_mmseqs_db=os.path.join(config['scratch_dir'], 'longest_orfs_mmseqs_db')
    log:
        os.path.join(config['log_dir'], 'protein', 'create_longorfs_mmseqs_db.log')
    shell:'''
    unset OMP_PROC_BIND
    module load MMseqs2
    mmseqs createdb {input.longORFs_pep} {output.longorfs_mmseqs_db} > {log} 2>&1
    '''

rule mmseqs_pfam_search:
    """
    Search for protein domains using mmseqs.
    """
    input:
        longORFs_pep=os.path.join(config['scratch_dir'], 'longest_orfs_mmseqs_db'),
        pfam_mmseqs_db=os.path.join(config['pfam_dir'], 'pfam_mmseqs_db')
    output:
        pfam_hits=os.path.join(config['output_dir'], 'assembly', 'protein', 'pfam_mmseqs_hits.m8')
    params:
        longORFs_dir=os.path.join(config['output_dir'], 'assembly', 'protein'),
        tmp=os.path.join(config['scratch_dir'], 'longorfs_pfam_tmp'),
        hits=os.path.join(config['scratch_dir'], 'pfam_mmseqs')
    threads: 60
    resources:
        mem_mb=200000
    log:
        os.path.join(config['log_dir'], 'protein', 'mmseqs_pfam_search.log')
    shell:'''
        unset OMP_PROC_BIND
        module load MMseqs2

        cd {params.longORFs_dir}

        mmseqs search {input.longORFs_pep} {input.pfam_mmseqs_db} {params.hits} {params.tmp} -s 7 -e 1e-2 --cov-mode 0 > {log} 2>&1

        mmseqs filterdb {params.hits} {params.hits}.filtered --extract-lines 10 > {log} 2>&1

        mmseqs convertalis {input.longORFs_pep} {input.pfam_mmseqs_db} {params.hits}.filtered {output.pfam_hits} > {log} 2>&1
        '''

rule predict_proteins:
    """
    Predict proteins using TransDecoder.Predict and retain mmseqs hits.
    """
    input:
        metatranscriptome=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.fasta'),
        pfam_hits=os.path.join(config['output_dir'], 'assembly', 'protein', 'pfam_mmseqs_hits.m8'),
        longORFs_pep=os.path.join(config['output_dir'], 'assembly', 'protein', 'final_metatranscriptome.fasta.transdecoder_dir', 'longest_orfs.pep')
    output:
        bed=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.bed'),
        cds=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.cds'),
        gff3=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.gff3'),
        pep=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.pep')
    params:
        longORFs_dir=os.path.join(config['output_dir'], 'assembly', 'protein'),
    log:
        os.path.join(config['log_dir'], 'protein', 'predict_proteins.log')
    shell:'''
        unset OMP_PROC_BIND
        module load TransDecoder
        module load Perl-bundle-CPAN/5.36.1-GCCcore-12.3.0 

        cd {params.longORFs_dir}

        TransDecoder.Predict -t {input.metatranscriptome} \
                --single_best_only \
                --retain_pfam_hits {input.pfam_hits} > {log} 2>&1

        mv {params.longORFs_dir}/final_metatranscriptome.fasta.transdecoder.bed {output.bed}
        mv {params.longORFs_dir}/final_metatranscriptome.fasta.transdecoder.cds {output.cds}
        mv {params.longORFs_dir}/final_metatranscriptome.fasta.transdecoder.gff3 {output.gff3}
        mv {params.longORFs_dir}/final_metatranscriptome.fasta.transdecoder.pep {output.pep}
    '''