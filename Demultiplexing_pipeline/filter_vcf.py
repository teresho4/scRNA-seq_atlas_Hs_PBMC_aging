intersect_file='output/only_unique_snp.tsv'
list_chr_pos = list()
with open(intersect_file) as f:
    for line in f:
        subline = line.split('\t')[0:2]
        subline[1] = subline[1].replace('\n','')
        list_chr_pos.append(subline)

outF = open('output/only_unique_snp.vcf' , 'w')
with open('souporcell_merged_sorted_vcf.vcf') as f:
    for line in f:
        chr_pos = line.split('\t')[0:2]
        if chr_pos in list_chr_pos:
            outF.write(line)
outF.close()
