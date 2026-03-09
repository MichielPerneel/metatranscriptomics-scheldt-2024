rule kallisto_index:
    input:
        fasta = os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.fasta'),
    output:
        index = os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.idx'),
    log: os.path.join(config['log_dir'], 'quantification', 'index.log'),
    threads: 8
    resources:
        mem_mb = 200000
    wrapper: 'v1.15.0/bio/kallisto/index'

rule kallisto_quant:
    input:
        fastq = [
            os.path.join(config['scratch_dir'], 'cleanup', '{sample}_1.rRNA_removed.fastq.gz'),
            os.path.join(config['scratch_dir'], 'cleanup', '{sample}_2.rRNA_removed.fastq.gz')
        ],
        index = os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.idx'),
    output: directory(os.path.join(config['output_dir'], 'quantification', '{sample}')),
    log: os.path.join(config['log_dir'], 'quantification', '{sample}.log'),
    threads: 24
    resources:
        mem_mb = 400000
    wrapper: 'v1.15.0/bio/kallisto/quant'