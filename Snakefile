configfile: "config.yaml"

import os
import sys
# Add 'scripts' and its subdirectories to Python's module search path to allow importing custom modules. 
sys.path.insert(1, 'scripts')

from utils import load_config, extract_samplenames, check_config

# Load configuration.
config = load_config('config.yaml')

# Load samples.
samples = extract_samplenames(config["samples"])
print("Samples: ", samples)

# Check setup
check_config(config)

# Include rules from separate modules.
include: "rules/quality_control.smk"
include: "rules/preprocessing.smk"
include: "rules/assembly.smk"
include: "rules/cluster_assemblies.smk"
include: "rules/protein_prediction.smk"
include: "rules/taxonomic_annotation.smk"
include: "rules/functional_annotation.smk"
include: "rules/quantification.smk"

# Run rules.
rule all:
    input:
       expand(os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.html'),
              sample=samples, num=[1, 2]),
       expand(os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.zip'),
              sample=samples, num=[1, 2]),
       #os.path.join(config['output_dir'], 'quality_control', 'multiqc_report.html'),
       expand(os.path.join(config['output_dir'], "assembly", "rnaSPAdes", "{sample}", "transcripts.fasta"), sample=samples),
       os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', 'final_metatranscriptome.fasta'),
       os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukprot', 'eukprot_annotation.m8'),
       os.path.join(config['output_dir'], 'annotation', 'functional_eggnog', 'functional_annotation.emapper.annotations'),
       expand(os.path.join(config['output_dir'], 'quantification', '{sample}'), sample=samples),
       os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukulele', 'EUKulele_done.txt')