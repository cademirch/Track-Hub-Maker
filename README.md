# Track-Hub-Maker
A Snakemake workflow and associated scripts to make population genomic statistics tracks for the UCSC Genome browser.

# Installation and Setup
1. Clone the repository: `git clone https://github.com/cademirch/Track-Hub-Maker`  
2. Create and activate conda environment: ```conda create --name <env> --file <req.txt> && conda activate <env>```
3. Create a tab seperated file `samples.tsv` in which the first column is the sample name and the second column is the accession of the genome on which the tracks are to be displayed. The sample name needs to be the same as the name of the directory that contains the raw data.
4. Ensure your directory structure looks like so: 
```
├── helper_files
│   └── trackDb.txt
├── README.md
├── req.txt
├── sample_1
│   └── s1.vcf.gz
├── sample_2
│   └── s2.vcf.gz
├── samples.tsv
├── scripts
│   ├── bedGraphToBigWig
│   ├── bedSort
│   ├── bedToBigBed
│   ├── makebed.py
│   └── write_hub_files.py
└── Snakefile
```
5. Run the Snakemake pipeline : `snakemake`
