rule cluster_assemblies:
        input: expand(os.path.join(config['output_dir'], 'assembly' ,'rnaSPAdes', '{sample}', 'transcripts.fasta'), sample=samples)
        output:
                metatranscriptome=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'metatranscriptome.fasta'),
        params:
                percent_id = "0.98",
                mmseqs_clu=os.path.join(config['scratch_dir'], 'mmseqs_clustering', 'clusterDB'),
                tmp=os.path.join(config['scratch_dir'], 'mmseqs_clustering')
        log: os.path.join(config['log_dir'], 'cluster_assemblies.log')
        shell:'''
        module load MMseqs2
        
        unset OMP_PROC_BIND        
        
        mmseqs easy-linclust {input} {params.mmseqs_clu} {params.tmp} --min-seq-id {params.percent_id} > {log} 2>&1

        mv {params.mmseqs_clu}_rep_seq.fasta {output.metatranscriptome}
        '''

rule rename_transcripts:
        input:  os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'metatranscriptome.fasta')
        output:  os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.fasta')
        params: length = 300
        log: os.path.join(config['log_dir'], 'rename_transcripts.log')
        conda: "anvio-8"
        shell:'''
        anvi-script-reformat-fasta {input} -o {output} -l {params.length} --simplify-names > {log} 2>&1
        '''