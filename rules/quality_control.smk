rule fastqc:
        input:
                os.path.join(config['raw_reads'], '{sample}_R{num}.fastq.gz')
        output:
                html = os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.html'),
                zip = os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.zip')
        log:
                os.path.join(config['output_dir'],'logs', 'fastqc', '{sample}_R{num}.log')
        wrapper:
                '0.27.1/bio/fastqc'

rule multiqc:
        input:
            qc_files = expand(os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.zip'),
                              sample=samples, num=[1, 2])
        output:
            report = os.path.join(config['output_dir'], 'quality_control', 'multiqc_report.html')
        params:
            outdir = os.path.join(config['output_dir'], 'quality_control'),
            tmpdir = os.path.join(config['scratch_dir'], 'tmp_multiqc')
        log:
            os.path.join(config['log_dir'], 'multiqc', 'multiqc.log')
        shell:"""
        module load MultiQC
        export TMPDIR={params.tmpdir}
        mkdir -p $TMPDIR

        multiqc {input.qc_files} \
            --force \
            -o {params.outdir} \
            -n multiqc_report.html \
            > {log} 2>&1
        """