rule trimmomatic:
        '''
        This step removes sequencing adapters using Trimmomatic. 
        '''
        input:
                r1 = os.path.join(config['raw_reads'], "{sample}_R1_001.fastq.gz"),
                r2 = os.path.join(config['raw_reads'], "{sample}_R2_001.fastq.gz"),
        output:
                r1 = os.path.join(config['scratch_dir'], "trimmed",
                                "{sample}_1.trimmed.fastq.gz"),
                r2 = os.path.join(config['scratch_dir'], "trimmed",
                                "{sample}_2.trimmed.fastq.gz"),
                # reads where trimming entirely removed the mate
                r1_unpaired = os.path.join(config['scratch_dir'], "trimmed",
                                "{sample}_1.unpaired.fastq.gz"),
                r2_unpaired = os.path.join(config['scratch_dir'], "trimmed",
                                "{sample}_2.unpaired.fastq.gz")
        params:
                trimmer=["ILLUMINACLIP:{}:2:30:7".format(config['adapters']), 
                        "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:50"]
        shell:'''
        module load  Trimmomatic/0.39-Java-11
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
                {input.r1} {input.r2} \
                {output.r1} {output.r1_unpaired} \
                {output.r2} {output.r2_unpaired} \
                {params.trimmer}
        '''

rule sample_rRNA_cleanup:
        '''
        This step removes ribosomal RNA reads using ribodetector.
        '''
        input:
            r1 = os.path.join(config['scratch_dir'], 'trimmed',
                            '{sample}_1.trimmed.fastq.gz'),
            r2 = os.path.join(config['scratch_dir'], 'trimmed',
                            '{sample}_2.trimmed.fastq.gz')
        output:
            r1 = os.path.join(config['scratch_dir'], 'cleanup',
                            '{sample}_1.rRNA_removed.fastq.gz'),
            r2 = os.path.join(config['scratch_dir'], 'cleanup',
                            '{sample}_2.rRNA_removed.fastq.gz'),
            rna1 = os.path.join(config['scratch_dir'], 'rRNA',
                            '{sample}_1.rRNA.fastq.gz'),
            rna2 = os.path.join(config['scratch_dir'], 'rRNA',
                            '{sample}_2.rRNA.fastq.gz')
        conda: os.path.join(config['conda_environments'], 'ribodetector.yaml')
        log:
            os.path.join(config['log_dir'], 'cleanup', '{sample}_ribodetector.log')
        shell:'''
            ribodetector_cpu -t 20 \
                -i {input.r1} {input.r2} \
                -o {output.r1} {output.r2} \
                -r {output.rna1} {output.rna2} \
                -e rrna -l 140 --chunk_size 256 > {log} 2>&1
            '''