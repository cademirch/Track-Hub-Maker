import glob
import gzip
import sys



def make_chrom_table(directory):
    contig_lengths = {}
    with gzip.open(f"{directory}/filtered.merged.vcf.gz", 'rt') as f:
        with open(f"{directory}/chrom_sizes", 'w') as out:
            for line in f:
                if line.startswith("##contig"):
                    line = line.strip().split("=")
                    contig = line[2].split(',')[0]
                    length = line[3].strip('>')
                    contig_lengths[contig] = int(length)
                    out.write(f"{contig}\t{length}\n")
                else:
                    continue
    return contig_lengths

def tajima(tajima_file, chr_table,bin_size):

    tajima_out = open(f"{tajima_file}.bg", "w")
    with open(tajima_file) as inp:
        next(inp)

        for line in inp:

            line = line.strip().split()

            chrom = line[0]
            if chrom not in chr_table:
                pass
            else:
                start = int(line[1])

                end = start + (bin_size-1)
                if end >= chr_table[chrom]:
                    end = chr_table[chrom]-1
                tajima_int = line[3]

                tajima_out.write(f'{chrom}\t{start}\t{end}\t{tajima_int}\n')
    tajima_out.close()
def snpden(snpden_file, chr_table,bin_size):

    snp_out = open(f"{snpden_file}.bg", "w")
    with open(snpden_file) as inp:
        next(inp)

        for line in inp:

            line = line.strip().split()

            chrom = line[0]
            if chrom not in chr_table:
                pass
            else:
                start = int(line[1])

                end = start + (bin_size-1)
                if end >= chr_table[chrom]:
                    end = chr_table[chrom]-1
                snp_int = line[2]

                snp_out.write(f'{chrom}\t{start}\t{end}\t{snp_int}\n')
    snp_out.close()

def pi(pi_file, chr_table):
    pi_out = open(f"{pi_file}.bg", "w")
    with open(pi_file) as inp:
        next(inp)
        for line in inp:
            line = line.strip().split()
            chrom = line[0]
            if chrom not in chr_table:
                pass
            else:

                start = line[1]
                end = int(line[2])
                if end >= chr_table[chrom]:
                    end = chr_table[chrom]-1

                pi_value = line[4]

                pi_out.write(f'{chrom}\t{start}\t{end}\t{pi_value}\n')
    pi_out.close()

def main():
    
    in_tajima = sys.argv[1]
    in_snpden = sys.argv[2]
    in_pi = sys.argv[3]
    step_size = [1000, 10000, 100000]
    steps = ['1k', '10k', '100k']
    directory = in_tajima.split("/")[0]
    table_file = make_chrom_table(directory)

    for i in range(len(step_size)):
        in_tajima = f"{directory}/{steps[i]}.tajima"
        in_snpden = f"{directory}/{steps[i]}.snpden"
        in_pi = f"{directory}/{steps[i]}.pi"
        tajima(in_tajima,table_file,step_size[i])
        snpden(in_snpden,table_file,step_size[i])
        pi(in_pi,table_file)
if __name__ == '__main__':
    main()

        


