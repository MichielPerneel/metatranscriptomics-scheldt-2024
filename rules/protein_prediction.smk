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
    input:
        metatranscriptome=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.fasta')
    output:
        longORFs_pep=os.path.join(config['output_dir'], 'assembly', 'protein',
                                  'final_metatranscriptome.fasta.transdecoder_dir', 'longest_orfs.pep')
    params:
        protein_dir=os.path.join(config['output_dir'], 'assembly', 'protein'),
        min_aa=150
    log:
        os.path.join(config['log_dir'], 'protein', 'predict_long_orfs.log')
    shell: r"""
        unset OMP_PROC_BIND
        module load TransDecoder
        module load Perl-bundle-CPAN/5.36.1-GCCcore-12.3.0
        set -euo pipefail

        cd {params.protein_dir}
        TransDecoder.LongOrfs -t {input.metatranscriptome} -m {params.min_aa} > {log} 2>&1
    """

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
    input:
        query_db=os.path.join(config['scratch_dir'], 'longest_orfs_mmseqs_db'),
        pfam_db=os.path.join(config['pfam_dir'], 'pfam_mmseqs_db')
    output:
        m8=os.path.join(config['output_dir'], 'assembly', 'protein', 'pfam_mmseqs_hits.m8')
    params:
        tmp=os.path.join(config['scratch_dir'], 'longorfs_pfam_tmp'),
        hits=os.path.join(config['scratch_dir'], 'pfam_mmseqs_hits_db')
    threads: 60
    resources:
        mem_mb=200000
    log:
        os.path.join(config['log_dir'], 'protein', 'mmseqs_pfam_search.log')
    shell: r"""
        unset OMP_PROC_BIND
        module load MMseqs2
        set -euo pipefail

        mmseqs search {input.query_db} {input.pfam_db} {params.hits} {params.tmp} \
            -s 7 -e 1e-2 --cov-mode 0 --threads {threads} > {log} 2>&1

        # keep only best hit per query to reduce file size
        mmseqs filterdb {params.hits} {params.hits}.filtered --extract-lines 1 >> {log} 2>&1

        mmseqs convertalis {input.query_db} {input.pfam_db} {params.hits}.filtered {output.m8} >> {log} 2>&1
    """

rule pfam_hit_transcript_ids:
    input:
        m8=os.path.join(config['output_dir'], 'assembly', 'protein', 'pfam_mmseqs_hits.m8')
    output:
        ids=os.path.join(config['output_dir'], 'assembly', 'protein', 'pfam_hit_transcript_ids.txt')
    shell: r"""
    cut -f1 {input.m8} \
    sed -E 's/\.p[0-9]+$//' \
    sort -u > {output.ids}
    """

rule filter_transcriptome_by_pfam:
    input:
        fasta=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.fasta'),
        ids=os.path.join(config['output_dir'], 'assembly', 'protein', 'pfam_hit_transcript_ids.txt')
    output:
        fasta_filt=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.pfamhit.fasta')
    log:
        os.path.join(config['log_dir'], 'protein', 'filter_transcriptome_by_pfam.log')
    shell: r"""
    module load SeqKit
    seqkit grep -f {input.ids} {input.fasta} > {output.fasta_filt} 2> {log}
    """

rule predict_proteins:
    input:
        metatranscriptome=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.pfamhit.fasta')
    output:
        bed=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.bed'),
        cds=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.cds'),
        gff3=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.gff3'),
        pep=os.path.join(config['output_dir'], 'assembly', 'protein', 'metatranscriptome.pep')
    params:
        workdir=os.path.join(config['output_dir'], 'assembly', 'protein', 'pfamhit_transdecoder')
    log:
        os.path.join(config['log_dir'], 'protein', 'predict_proteins_pfamhit.log')
    threads: 16
    resources:
        mem_mb=200000
    shell: r"""
    unset OMP_PROC_BIND
    module load TransDecoder
    module load Perl-bundle-CPAN/5.36.1-GCCcore-12.3.0
    set -euo pipefail

    mkdir -p {params.workdir}
    cd {params.workdir}

    TransDecoder.LongOrfs -t {input.metatranscriptome} > {log} 2>&1
    TransDecoder.Predict -t {input.metatranscriptome} --single_best_only >> {log} 2>&1

    mv final_metatranscriptome.pfamhit.fasta.transdecoder.bed  {output.bed}
    mv final_metatranscriptome.pfamhit.fasta.transdecoder.cds  {output.cds}
    mv final_metatranscriptome.pfamhit.fasta.transdecoder.gff3 {output.gff3}
    mv final_metatranscriptome.pfamhit.fasta.transdecoder.pep  {output.pep}
    """