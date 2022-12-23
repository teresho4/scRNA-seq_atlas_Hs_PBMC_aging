library(dplyr)
library(tidyr)

`%ni%` <- Negate(`%in%`)
patients <- basename(list.dirs('output/', full.names = TRUE, recursive = TRUE))
patients <- patients[!(patients %in% "output")]

all_vcf_tables <- NULL
for(patient in patients){
  vcf_table <- read.csv(paste('./output/', patient,'/vcf_table.tsv', sep=''))
  all_vcf_tables <- rbind(all_vcf_tables, vcf_table)
}

duplicates <- all_vcf_tables[duplicated(all_vcf_tables$CHROM.POS.unknown.GT),]
all_vcf_tables <- all_vcf_tables %>% filter(CHROM.POS.unknown.GT %ni% duplicates)

table_unique <- all_vcf_tables %>% separate(CHROM.POS.unknown.GT, c('CHROM', 'POS', 'GT'), sep = '\t')
table_unique$Patient <- NULL
table_unique$GT <- NULL
write.table(table_unique, 'output/only_unique_snp.tsv', row.names = F, quote = F, sep='\t')