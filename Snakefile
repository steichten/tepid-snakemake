# hard-coded locations of things
SAMPLES = glob_wildcards("raw/{S}_R1_001.fastq.gz").S
REFERENCE = "/home/steve/tepid_test/genome/Bd21Control_SNPincorp_sgr1_genome"
TEDB = "/home/steve/bin/TEPID/Annotation/Brachypodium/Brachy_TE_v2.2.bed.gz"
SCRIPTS = "/home/steve/bin/TEPID/Scripts"
REFLINE = "Bd21-t1-WGS_S1"
# bash safe mode
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")


rule final_output:
  input: 
     "genotyped_insertions.bed",
     "genotyped_deletions.bed"

  
  #setup the enviroment and required dependencies and software

#tepid-map (bowtie2 and yaha) 
rule tepid_map:
    input:
        r1 = 'raw/{sample}_R1_001.fastq.gz',
        r2 = 'raw/{sample}_R2_001.fastq.gz'
    output:
        bam   = 'align/{sample}/{sample}.bam',
        bai   = 'align/{sample}/{sample}.bam.bai',
        split = 'align/{sample}/{sample}.split.bam',
        umap  = 'align/{sample}/{sample}.umap.fastq'
    params:
        ref=REFERENCE,
    shell:
          "tepid-map -x {params.ref} -y {params.ref}.X15_01_65525S -p 12 -s 200 -n {wildcards.sample} -1 {input.r1} -2 {input.r2} &&"
          " mv {wildcards.sample}.bam {output.bam} &&"
          " mv {wildcards.sample}.bam.bai {output.bai} &&"
          " mv {wildcards.sample}.split.bam {output.split} &&"
          " mv {wildcards.sample}.umap.fastq {output.umap}"
  
#tepid-discover
rule tepid_discover:
     input:
        bam   = 'align/{sample}/{sample}.bam',
        split = 'align/{sample}/{sample}.split.bam'
     output:
        inbed    = 'align/{sample}/insertions_{sample}.bed',
        inreads  = 'align/{sample}/insertion_reads_{sample}.txt',
        delbed   = 'align/{sample}/deletions_{sample}.bed',
        delreads = 'align/{sample}/deletion_reads_{sample}.txt',
        logs     = 'logs/tepid_discover_log_{sample}.txt'
     params:
        tedb = TEDB
     shell:
        "tepid-discover -p 12 -n {wildcards.sample} -c {input.bam} -s {input.split} -t {params.tedb} &&"
        "mv insertions_{wildcards.sample}.bed {output.inbed} &&"
        "mv insertion_reads_{wildcards.sample}.txt {output.inreads} &&"
        "mv deletions_{wildcards.sample}.bed {output.delbed} &&"
        "mv deletion_reads_{wildcards.sample}.txt {output.delreads} &&"
        "mv tepid_discover_log_{wildcards.sample}.txt {output.logs}"


#merge
rule sample_bed_merge:
    input:
        expand("align/{sample}/insertions_{sample}.bed",sample=SAMPLES),
        expand("align/{sample}/deletions_{sample}.bed",sample=SAMPLES)
    output:
        inbed  = "align/insertions.bed",
        inpoly = "align/insertions_poly_te.bed",
        delbed  = "align/deletions.bed"
    params:
        scripts = SCRIPTS
    shell:
          "cd align &&"
          "python {params.scripts}/merge_insertions.py -f insertions &&"
          "python {params.scripts}/merge_deletions.py -f deletions &&"
          "cd ../"

#refine
rule tepid_refine:
    input:
        inbed  = "align/insertions.bed",
        inpoly = "align/insertions_poly_te.bed",
        delbed  = "align/deletions.bed",
        bam   = 'align/{sample}/{sample}.bam',
        split = 'align/{sample}/{sample}.split.bam'
    output:
        ambin      = "align/{sample}/ambiguous_insertion_{sample}.bed",
        ambdel     = "align/{sample}/ambiguous_deletion_{sample}.bed",
        spin       = "align/{sample}/second_pass_insertion_{sample}.bed",
        spinreads  = "align/{sample}/second_pass_reads_insertion_{sample}.txt",
        spdel      = "align/{sample}/second_pass_deletion_{sample}.bed",
        spdelreads = "align/{sample}/second_pass_reads_deletion_{sample}.txt"
    params:
        tedb = TEDB,
        samples = SAMPLES,
        inbed  = "../insertions.bed",
        inpoly = "../insertions_poly_te.bed",
        delbed  = "../deletions.bed",
        bam   = '{sample}.bam',
        split = '{sample}.split.bam'
    shell:
          "cd align/{wildcards.sample} &&"
          "tepid-refine "
          " -i {params.inbed}"
          " -d {params.delbed}"
          " -p 12"
          " -t {params.tedb}"
          " -n {wildcards.sample}"
          " -c {params.bam}"
          " -s {params.split}"
          " -a <(echo {params.samples} | tr ' ' '\n')"
          #" mv ambiguous_insertion_{wildcards.sample}.bed {output.ambin} &&"
          #" mv ambiguous_deletion_{wildcards.sample}.bed {output.ambdel} &&"
          #" mv second_pass_insertion_{wildcards.sample}.bed {output.spin} &&"
          #" mv second_pass_reads_insertion_{wildcards.sample}.txt {output.spinreads} &&"
          #" mv second_pass_deletion_{wildcards.sample}.bed {output.spdel} &&"
          #" mv second_pass_reads_deletion_{wildcards.sample}.txt {output.spdelreads}"

#cat first and second pass
rule combine_bed:
    input:
        inbed      = "align/{sample}/insertions_{sample}.bed",
        delbed     = "align/{sample}/deletions_{sample}.bed",
        spin       = "align/{sample}/second_pass_insertion_{sample}.bed",
        spdel      = "align/{sample}/second_pass_deletion_{sample}.bed"
    output:
        refin  = "align/{sample}/refined_insertions_{sample}.bed",
        refdel = "align/{sample}/refined_deletions_{sample}.bed"
    params:
        scripts = SCRIPTS
    shell:
          "cat {input.spin} {input.inbed} > {output.refin} && "
          "cat {input.spdel} {input.delbed} > {output.refdel}"
          
#merge2
rule sample_bed_merge2:
    input:
        expand("align/{sample}/refined_insertions_{sample}.bed",sample=SAMPLES),
        expand("align/{sample}/refined_deletions_{sample}.bed",sample=SAMPLES),
        expand("align/{sample}/ambiguous_insertion_{sample}.bed",sample=SAMPLES),
        expand("align/{sample}/ambiguous_deletion_{sample}.bed",sample=SAMPLES)
    output:
        refin     = "align/refined_insertions.bed",
        refinpoly = "align/refined_insertions_poly_te.bed",
        refdel    = "align/refined_deletions.bed",
        ambin     = "align/ambiguous_insertion.bed",
        ambdel    = "align/ambiguous_deletion.bed",
    params:
        scripts = SCRIPTS
    shell:
          "cd align && "
          "python {params.scripts}/merge_insertions.py -f refined_insertions && "
          "python {params.scripts}/merge_deletions.py -f refined_deletions && "
          "python {params.scripts}/merge_insertions.py -f ambiguous_insertion && "
          "python {params.scripts}/merge_deletions.py -f ambiguous_deletion && "
          "cd ../"

#genotype
rule tepid_genotype:
    input:
        refin     = "align/refined_insertions.bed",
        refdel    = "align/refined_deletions.bed",
        ambin     = "align/ambiguous_insertion.bed",
        ambdel    = "align/ambiguous_deletion.bed",
    output:
        genoin  = "genotyped_insertions.bed",
        genodel = "genotyped_deletions.bed"
    params:
        scripts = SCRIPTS,
        samples = SAMPLES,
        refline = REFLINE
    shell:
        "python {params.scripts}/genotype.py "
        " -d  "
        " -a {input.ambdel} "
        " -m {input.refdel} "
        " -s <(echo {params.samples} | tr ' ' '\n') "
        " -r {params.refline} > {output.genodel} &&"
        " python {params.scripts}/genotype.py "
        " -i  "
        " -a {input.ambin} "
        " -m {input.refin} "
        " -s <(echo {params.samples} | tr ' ' '\n') "
        " -r {params.refline} > {output.genoin}" 
        
#
