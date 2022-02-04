library(tidyverse)
library(microseq)



#################################################
### Downloading the assembly_summary_refseq.txt
###
ok <- download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt",
                    destfile = "assembly_summary_archaea.txt")
ok <- download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
                    destfile = "assembly_summary_bacteria.txt")



######################
### Reading the table
###
read_delim("assembly_summary_archaea.txt", delim = "\t", skip = 1, quote = "") %>%
  rename(assembly_accession = `# assembly_accession`) -> arc.tbl
read_delim("assembly_summary_bacteria.txt", delim = "\t", skip = 1, quote = "") %>%
  rename(assembly_accession = `# assembly_accession`) -> bac.tbl
bac.tbl %>%
  bind_rows(arc.tbl) %>%
  filter(refseq_category == "reference genome") %>%
  mutate(file_prefix = basename(ftp_path)) %>%
  select(-wgs_master, -isolate, -excluded_from_refseq) -> rfq.tbl
save(rfq.tbl, file = "rfq.tbl.RData")



################
### Downloading
###
for(i in 1:nrow(rfq.tbl)){
  fna.file <- str_c(rfq.tbl$file_prefix[i], "_genomic.fna.gz")
  download.file(file.path(rfq.tbl$ftp_path[i], fna.file),
                destfile = file.path("fna", fna.file))

  faa.file <- str_c(rfq.tbl$file_prefix[i], "_protein.faa.gz")
  download.file(file.path(rfq.tbl$ftp_path[i], faa.file),
                destfile = file.path("faa", faa.file))

  gff.file <- str_c(rfq.tbl$file_prefix[i], "_genomic.gff.gz")
  download.file(file.path(rfq.tbl$ftp_path[i], gff.file),
                destfile = file.path("gff", gff.file))
  R.utils::gunzip(file.path("gff", gff.file))
}



#####################
### Updating rfq.tbl
###
load("rfq.tbl.RData")
rfq.tbl %>%
   mutate(GC = 0) -> rfq.tbl
for(i in 1:nrow(rfq.tbl)){
  readFasta(file.path("fna", str_c(rfq.tbl$file_prefix[i], "_genomic.fna.gz"))) %>%
    mutate(N = str_length(str_remove_all(Sequence, "[^ACTG]"))) %>%
    mutate(N_GC = str_length(str_remove_all(Sequence, "[AT]"))) -> gnm
  rfq.tbl$GC[i] <- sum(gnm$N_GC) / sum(gnm$N)

}
save(rfq.tbl, file = "rfq.tbl.RData")



########################
### Cleaning GFF-tables
###
load("rfq.tbl.RData")
rfq.tbl %>% 
  mutate(trans_table = "") -> rfq.tbl
for(i in 1:nrow(rfq.tbl)){
  cat("Genome", i, "\n")
  readFasta(file.path("faa", str_c(rfq.tbl$file_prefix[i], "_protein.faa.gz"))) %>% 
    mutate(Tag = word(Header, 1, 1)) -> faa
  
  suppressWarnings(readGFF(file.path("gff", str_c(rfq.tbl$file_prefix[i], "_genomic.gff")))) %>% 
    filter(Type == "CDS") %>% 
    mutate(Tag = str_extract(Attributes, "ID=.+?;")) %>% 
    mutate(Tag = str_remove_all(Tag, "ID=cds-|;")) %>% 
    arrange(Tag) -> gff
  
  rfq.tbl$trans_table[i] <- str_c(unique(str_remove(str_extract(gff$Attributes, "transl_table=.+"), "transl_table=")), collapse = ",")
  
  gff %>% 
    full_join(faa, by = "Tag") %>% 
    select(-Tag) -> gffaa.tbl
  save(gffaa.tbl, file = file.path("gff", str_c(rfq.tbl$file_prefix[i], ".gff.faa.RData")))
  cat("   found", sum(is.na(gffaa.tbl$Header)), "annotations without protein sequence\n")
  cat("   found", sum(is.na(gffaa.tbl$Seqid)), "protein sequences without annotations\n")
}
save(rfq.tbl, file = "rfq.tbl.RData")
