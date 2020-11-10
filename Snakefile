import glob

SAMPLES_FILE = ''

genomes = {}

with open(SAMPLES_FILE, 'r') as f:
    for line in f:
        line = line.strip().split("\t")
        
        if line[0] not in genomes:
            genomes[line[0]] = line[1]


file_types = ["tajima", "snpden", "pi"]

rule all:
    input: 
        [f"{key}/{value}/1k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/{value}/10k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/{value}/100k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/{value}/filtered.merged.vcf.gz" for key, value in genomes.items()],
        [f"{key}/{value}/filtered.merged.vcf.gz.tbi" for key, value in genomes.items()],
        [f"{key}/{value}/{cvg}.bb" for key, value in genomes.items() for cvg in ['clean', 'low', 'high']],
        [f"{key}/{value}/trackDb.txt" for key, value in genomes.items()],
rule move:
    input:
        [f"{key}/1k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/10k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/100k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/filtered.merged.vcf.gz" for key, value in genomes.items()],
        [f"{key}/filtered.merged.vcf.gz.tbi" for key, value in genomes.items()],
        [f"{key}/{cvg}.bb" for key, value in genomes.items() for cvg in ['clean', 'low', 'high']]  
    output:
        [f"{key}/{value}/1k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/{value}/10k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/{value}/100k.{f}.bw" for key, value in genomes.items() for f in file_types],
        [f"{key}/{value}/filtered.merged.vcf.gz" for key, value in genomes.items()],
        [f"{key}/{value}/filtered.merged.vcf.gz.tbi" for key, value in genomes.items()],
        [f"{key}/{value}/{cvg}.bb" for key, value in genomes.items() for cvg in ['clean', 'low', 'high']]  

    run:
        comb = zip(input,output)

        for i in comb:
            command = f"mv {i[0]} {i[1]}"
            shell(command)

rule write_hub_files:
    output:
        trackdb = ("{directory}/{genome}/trackDb.txt"),
    shell:
         "python3 scripts/write_hub_files.py {wildcards.directory} {wildcards.genome}"



rule concat_vcfs:
    input:
        ("{directory}")
    output:
        temp("{directory}/merged.vcf.gz")
    log:
        "logs/bcftools/{directory}/mergevcfs.log"
    params:
        output_type="-O z"
    
    shell:
        "bcftools concat {params.output_type} {input}/*.vcf.gz > {output} 2> {log}"

rule filter_vcf:
    input:
        "{directory}/merged.vcf.gz"
    output:
        "{directory}/filtered.merged.vcf.gz"
    log:
        "logs/vcftools/{directory}/filter_vcf.log"
    shell:
        "vcftools --gzvcf {input} --remove-filtered-all --recode --stdout 2>>{log} | bgzip > {output}"

rule tabix:
    input:
        "{directory}/filtered.merged.vcf.gz"
    output:
        "{directory}/filtered.merged.vcf.gz.tbi"
    params:
        p = "-p vcf"
    shell:
        "tabix {params.p} {input}"
rule calc_tajima:
    input:
        "{directory}/filtered.merged.vcf.gz"
    output:
        one_t=temp("{directory}/1k.tajima"),
        ten_t=temp("{directory}/10k.tajima"),
        hundred_t=temp("{directory}/100k.tajima"),
    log:
        "logs/vcftools/{directory}/calc_tajima.log"
    shell:
        """
        vcftools --gzvcf {input} --TajimaD 1000 --stdout > {output.one_t} 2>> {log}
        vcftools --gzvcf {input} --TajimaD 10000 --stdout > {output.ten_t} 2>> {log}
        vcftools --gzvcf {input} --TajimaD 100000 --stdout > {output.hundred_t} 2>> {log}
        """

rule calc_pi:
    input:
        "{directory}/filtered.merged.vcf.gz"
    output:
        one_p=temp("{directory}/1k.pi"),
        ten_p=temp("{directory}/10k.pi"),
        hundred_p=temp("{directory}/100k.pi"),
    log:
        "logs/vcftools/{directory}/calc_pi.log"
    shell:
        """
        vcftools --gzvcf {input} --window-pi 1000 --stdout > {output.one_p} 2>> {log}
        vcftools --gzvcf {input} --window-pi 10000 --stdout > {output.ten_p} 2>> {log}
        vcftools --gzvcf {input} --window-pi 100000 --stdout > {output.hundred_p} 2>> {log}
        """

rule calc_snpden:
    input:
        "{directory}/filtered.merged.vcf.gz"
    output:
        one_s=temp("{directory}/1k.snpden"),
        ten_s=temp("{directory}/10k.snpden"),
        hundred_s=temp("{directory}/100k.snpden")

    log:
        "logs/vcftools/{directory}/calc_snpden.log"
    shell:
        """
        vcftools --gzvcf {input} --SNPdensity 1000 --maf 0.001 --stdout > {output.one_s} 2>> {log}
        vcftools --gzvcf {input} --SNPdensity 10000 --maf 0.001 --stdout > {output.ten_s} 2>> {log}
        vcftools --gzvcf {input} --SNPdensity 100000 --maf 0.001 --stdout > {output.hundred_s} 2>> {log}
        """

rule convert_to_bed:
    input:
        one_t=("{directory}/1k.tajima"),
        one_s=("{directory}/1k.snpden"),
        one_p=("{directory}/1k.pi"),
        ten_t=("{directory}/10k.tajima"),
        ten_s=("{directory}/10k.snpden"),
        ten_p=("{directory}/10k.pi"),
        hund_t=("{directory}/100k.tajima"),
        hund_s=("{directory}/100k.snpden"),
        hund_p=("{directory}/100k.pi"),


    output:
        temp("{directory}/1k.tajima.bg"),
        temp("{directory}/1k.snpden.bg"),
        temp("{directory}/1k.pi.bg"),
        temp("{directory}/10k.tajima.bg"),
        temp("{directory}/10k.snpden.bg"),
        temp("{directory}/10k.pi.bg"),
        temp("{directory}/100k.tajima.bg"),
        temp("{directory}/100k.snpden.bg"),
        temp("{directory}/100k.pi.bg"),

    shell:
        "python3 scripts/makebed.py {input.one_t} {input.one_s} {input.one_p}"


rule bed_sort_1:
    input:
        one_t=("{directory}/1k.tajima.bg"),
        one_s=("{directory}/1k.snpden.bg"),
        one_p=("{directory}/1k.pi.bg"),


    output:
        one_t=temp("{directory}/1k.tajima.sorted.bg"),
        one_s=temp("{directory}/1k.snpden.sorted.bg"),
        one_p=temp("{directory}/1k.pi.sorted.bg"),


    shell:
        """
        scripts/bedSort {input.one_t} {output.one_t}
        scripts/bedSort {input.one_s} {output.one_s}
        scripts/bedSort {input.one_p} {output.one_p}
        """

rule bed_sort_10:
    input:
        ten_t=("{directory}/10k.tajima.bg"),
        ten_s=("{directory}/10k.snpden.bg"),
        ten_p=("{directory}/10k.pi.bg"),

    output:
        ten_t=temp("{directory}/10k.tajima.sorted.bg"),
        ten_s=temp("{directory}/10k.snpden.sorted.bg"),
        ten_p=temp("{directory}/10k.pi.sorted.bg"),

    shell:
        """
        scripts/bedSort {input.ten_t} {output.ten_t}
        scripts/bedSort {input.ten_s} {output.ten_s}
        scripts/bedSort {input.ten_p} {output.ten_p}
        """


rule bed_sort_100:
    input:
        hund_t=("{directory}/100k.tajima.bg"),
        hund_s=("{directory}/100k.snpden.bg"),
        hund_p=("{directory}/100k.pi.bg"),

    output:
        hund_t=temp("{directory}/100k.tajima.sorted.bg"),
        hund_s=temp("{directory}/100k.snpden.sorted.bg"),
        hund_p=temp("{directory}/100k.pi.sorted.bg"),

    shell:
        """
        scripts/bedSort {input.hund_t} {output.hund_t}
        scripts/bedSort {input.hund_s} {output.hund_s}
        scripts/bedSort {input.hund_p} {output.hund_p}
        """


rule bedtobigwig_1:
    input:
        one_t=("{directory}/1k.tajima.sorted.bg"),
        one_s=("{directory}/1k.snpden.sorted.bg"),
        one_p=("{directory}/1k.pi.sorted.bg"),


    output:
        one_t=("{directory}/1k.tajima.bw"),
        one_s=("{directory}/1k.snpden.bw"),
        one_p=("{directory}/1k.pi.bw"),

    params:
        chroms="{directory}/chrom_sizes"
    shell:
        """
        scripts/bedGraphToBigWig {input.one_t} {params.chroms} {output.one_t}
        scripts/bedGraphToBigWig {input.one_s} {params.chroms} {output.one_s}
        scripts/bedGraphToBigWig {input.one_p} {params.chroms} {output.one_p}
        """     

rule bedtobigwig_10:
    input:
        ten_t=("{directory}/10k.tajima.sorted.bg"),
        ten_s=("{directory}/10k.snpden.sorted.bg"),
        ten_p=("{directory}/10k.pi.sorted.bg"),
    output:
        ten_t=("{directory}/10k.tajima.bw"),
        ten_s=("{directory}/10k.snpden.bw"),
        ten_p=("{directory}/10k.pi.bw"),
    params:
        chroms="{directory}/chrom_sizes"
    shell:
        """
        scripts/bedGraphToBigWig {input.ten_t} {params.chroms} {output.ten_t}
        scripts/bedGraphToBigWig {input.ten_s} {params.chroms} {output.ten_s}
        scripts/bedGraphToBigWig {input.ten_p} {params.chroms} {output.ten_p}
        """

rule bedtobigwig_100:
    input:
        hund_t=("{directory}/100k.tajima.sorted.bg"),
        hund_s=("{directory}/100k.snpden.sorted.bg"),
        hund_p=("{directory}/100k.pi.sorted.bg"),
    output:
        hund_t=("{directory}/100k.tajima.bw"),
        hund_s=("{directory}/100k.snpden.bw"),
        hund_p=("{directory}/100k.pi.bw"),
    params:
        chroms="{directory}/chrom_sizes"
    shell:
        """
        scripts/bedGraphToBigWig {input.hund_t} {params.chroms} {output.hund_t}
        scripts/bedGraphToBigWig {input.hund_s} {params.chroms} {output.hund_s}
        scripts/bedGraphToBigWig {input.hund_p} {params.chroms} {output.hund_p}
        """

rule coverage_beds: 
    input:
        d = "{directory}"

    output:
        clean = temp("{directory}/clean.bed"),
        low = temp("{directory}/low.bed"),
        high = temp("{directory}/high.bed")

    run:

        clean_file = glob.glob(f"{input.d}/*_clean_coverage_sites_merged.bed")[0]
        low_file = glob.glob(f"{input.d}/*_low_coverage_sites_merged.bed")[0]
        high_file = glob.glob(f"{input.d}/*_low_coverage_sites_merged.bed")[0]

        commands = [f"scripts/bedSort {clean_file} {output.clean}", f"scripts/bedSort {low_file} {output.low}", f"scripts/bedSort {high_file} {output.high}"]
    
        for c in commands:
            shell(c)

rule bedtobigbed:
    input:
        clean = ("{directory}/clean.bed"),
        low = ("{directory}/low.bed"),
        high = ("{directory}/high.bed")
    params:
        chroms = "{directory}/chrom_sizes"
    output:
        clean = ("{directory}/clean.bb"),
        low = ("{directory}/low.bb"),
        high = ("{directory}/high.bb")

    shell:
        """
        scripts/bedToBigBed {input.clean} {params.chroms} {output.clean}
        scripts/bedToBigBed {input.low} {params.chroms} {output.low}
        scripts/bedToBigBed {input.low} {params.chroms} {output.high}
        """
