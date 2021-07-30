rule all:
  
  input:
    "/home/jovyan/pipeline_data/SRR622461_1_fastp_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_2_fastp_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_sorted.bam.bai",
    "/home/jovyan/pipeline_data/SRR622461_sorted_flagstats.txt",
    "/home/jovyan/pipeline_data/SRR622461_refined.bam.bai",
    "/home/jovyan/pipeline_data/SRR622461_refined_flagstats.txt",
    "/home/jovyan/pipeline_data/SRR622461_CYP2C19.vcf"

rule fastp_filter:
  
  input:
    fwd="{sample}_1.fastq.gz",
    rev="{sample}_2.fastq.gz"
    
  output:
    fwd=temp("{sample}_1_fastp.fastq.gz"),
    rev=temp("{sample}_2_fastp.fastq.gz"),
    json="{sample}_fastp.json",
    html="{sample}_fastp.html",
    out="{sample}_fastp.txt"
    
  shell:
    "fastp --thread 8 -i {input.fwd} -o {output.fwd} \
    -I {input.rev} -O {output.rev} --disable_adapter_trimming \
    --length_required 36 -3 --correction --json {output.json} \
    --html {output.html} 2> {output.out}"

rule fastqc:
  
  input:
    
     "{sample}.fastq.gz"
    
  output:
    
     "{sample}_fastqc.html",
     "{sample}_fastqc.zip"
    
  shell:
     "fastqc -t 8 {input}"

rule bwa_mem:
  
  input:
     ref="/home/jovyan/pipeline_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
     fwd="{sample}_1_fastp.fastq.gz",
     rev="{sample}_2_fastp.fastq.gz"
    
  output:
     temp("{sample}.sam")
  
  shell:
     "bwa mem -t 8 -R '@RG\\tID:1\\tLB:library\\tPL:Illumina\\tPU:lane1\\tSM:NA12878' \
     {input.ref} {input.fwd} {input.rev} > {output}"

rule sam_to_bam:
  
  input:
     "{sample}.sam"
    
  output:
    temp("{sample}.bam")
    
  shell:
     "samtools view -b {input} -o {output} -@ 8"

rule samtools_sort:
  
  input:
     "mapped_reads/{sample}.bam"
      
  output:
     "sorted_reads/{sample}.bam"
      
  shell:
     "samtools sort -T sorted_reads/{wildcards.sample} "
     "-O bam {input} > {output}"
    
rule samtools_index:
  
  input:
     "sorted_reads/{sample}.bam"
      
  output:
     "sorted_reads/{sample}.bam.bai"
      
  shell:
     "samtools index {input}"

rule samtools_flagstat:
  
  input:
    "{sample}.bam"
    
  output:
    "{sample}_flagstats.txt"
    
  shell:
    "samtools flagstat {input} -@ 8 > {output}"

rule picard_remove_duplicates:
  
  input:
    "{sample}_sorted.bam"
    
  output:
    bam="{sample}_refined.bam",
    metrics="{sample}_dupl_metrics.txt"
    
  shell:
    "picard MarkDuplicates I={input} O={output.bam} \
    METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true"

rule gatk_haplotype_caller:
  
  input:
    ref="/home/jovyan/pipeline_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    bam="{sample}_refined.bam"
    
  output:
    vcf="{sample}.vcf",
    idx="{sample}.vcf.idx"
    
  shell:
    "gatk HaplotypeCaller -R {input.ref} \
    -I {input.bam} -O {output.vcf} -L 10 -mbq 20"

rule vcftools_filter:
  
  input:
    "{sample}.vcf"
    
  output:
    temp("{sample}.recode.vcf"),
    "{sample}.log"
    
  shell:
    "vcftools --vcf {input} --out {wildcards.sample} \
    --minDP 3 --minQ 20 --recode --recode-INFO-all"

rule vcftools_exclude:
  
  input:
    "{sample}.recode.vcf"
    
  output:
    temp("{sample}_filtered.recode.vcf"),
    "{sample}_filtered.log"
    
  shell:
    "vcftools --vcf {input} --out {wildcards.sample}_filtered \
    --max-missing 1 --recode --recode-INFO-all"

rule snpeff_annotate:
  
  input:
    "{sample}_filtered.recode.vcf"
    
  output:
    vcf="{sample}_annotated.vcf",
    html="{sample}_snpEff_summary.html",
    txt="{sample}_snpEff_summary.genes.txt"
    
  shell:
    "snpEff -Xmx4G GRCh38.99 {input} -stats {output.html} > {output.vcf}"

rule snpsift_filter:
  
  input:
    "{sample}_annotated.vcf"
    
  output:
    "{sample}_CYP2C19.vcf"
    
  shell:
    "grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X,Y,M]' 'CYP2C19'\" > {output}"
    
rule plot_quals:
  
  input:
     "calls/{output}all.vcf"
      
  output:
     "plots/{sample}quals.svg"
      
  script:
     "scripts/plot-quals.py"
      
      
