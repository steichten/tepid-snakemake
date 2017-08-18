# TEPID-snakemake workflow
snakemake pipeline for TEPID transposon variant calling

The `Snakefile` lists the input requirements at the top. You'll want:

- YAHA indexed genome
- bowtie2 indexed genome
- folder called `raw` containing your paired-end fastq files (edit as needed to glob the correct sample names in your set)
- annotation of transposable elements per TEPID format (I'm using the TEPID inclued brachypodium v1.2 annotations)
- location to the TEPID scripts
- ID the samplename that contains your reference sample for phasing of variants
