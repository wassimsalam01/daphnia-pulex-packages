if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("AnnotationDbi", "AnnotationHub", "Biostrings",
                       "BSgenome", "GenomicFeatures", "rtracklayer"))
# Set your path
setwd("/path/to/seed-file")

# BSgenome package
## Download fna.gz sequence and export as 2bit file
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/134/715/GCF_021134715.1_ASM2113471v1/GCF_021134715.1_ASM2113471v1_genomic.fna.gz",
              destfile = "daphnia_pulex.fasta.gz")
pulex_fasta = Biostrings::readDNAStringSet("daphnia_pulex.fasta.gz")
pulex_2bit = file.path(getwd(), "daphnia_pulex.2bit")
rtracklayer::export.2bit(pulex_fasta, pulex_2bit)

## Creating BSgenome.Dpulex.NCBI.ASM2113471v1 library
BSgenome::forgeBSgenomeDataPkg("BSgenome.Dpulex.NCBI.ASM2113471v1-seed", verbose = TRUE)
## Quit R, and build the source package (tarball) with...
system('R CMD build BSgenome.Dpulex.NCBI.ASM2113471v1/')
system('R CMD check BSgenome.Dpulex.NCBI.ASM2113471v1_1.0.0.tar.gz')
system('R CMD INSTALL BSgenome.Dpulex.NCBI.ASM2113471v1_1.0.0.tar.gz')

library(BSgenome.Dpulex.NCBI.ASM2113471v1)

# TxDb package
## Gene Transfer Format (GTF) file for building Transcript Database (TxDb)
download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/134/715/GCF_021134715.1_ASM2113471v1/GCF_021134715.1_ASM2113471v1_genomic.gtf.gz",
              destfile = "daphnia_pulex.gtf.gz")

## Chromosome sizes info file
download.file(url = "https://hgdownload.soe.ucsc.edu/hubs/GCF/021/134/715/GCF_021134715.1/GCF_021134715.1.chrom.sizes.txt",
              destfile = "daphnia_pulex.chrom.sizes.txt")
chrom_info = read.csv(file = "daphnia_pulex.chrom.sizes.txt", header = FALSE, sep = "\t", col.names = c("chr","size"))
chrom_info = rbind(chrom_info[order(chrom_info$chr[1:12]),],chrom_info[13,])

## Creating Transcript Database (TxDb)
seqinfo_Dpulex = GenomeInfoDb::Seqinfo(seqnames = chrom_info$chr, seqlengths = chrom_info$size, isCircular = logical(13), genome = "Dpulex")
## Build metadata dataframe
name = c("Resource URL", "Type of Gene ID", "exon_nrow", "cds_nrow")
value = c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/134/715/GCF_021134715.1_ASM2113471v1/", "Entrez Gene ID", "159649", "113453")

pulex_metadata = data.frame(name, value)

txdb_pulex = GenomicFeatures::makeTxDbFromGFF(file = "databases/daphnia_pulex.gtf.gz",
                                      format = "auto",
                                      dataSource = "NCBI",
                                      organism = "Daphnia pulex",
                                      taxonomyId = 6669,
                                      chrominfo = seqinfo_Dpulex,
                                      metadata = pulex_metadata)

GenomicFeatures::makeTxDbPackage(txdb_pulex,
                                 version = "1.0",
                                 maintainer = "Wassim Salam <wassimsalam49@gmail.com>",
                                 author = "Wassim Salam",
                                 destDir = ".",
                                 pkgname = "TxDb.Dpulex.NCBI.ASM2113471v1.knownGene")

install.packages("TxDb.Dpulex.NCBI.ASM2113471v1.knownGene", repos = NULL, type = "source")
library(TxDb.Dpulex.NCBI.ASM2113471v1.knownGene)

# org.Dp.eg.db package
ah = AnnotationHub::AnnotationHub()
ah = subset(ah, species == "Daphnia pulex")

org.Dp.eg.db = AnnotationHub::AnnotationHub()[["AH115573"]]
GIDkeys = keys(org.Dp.eg.db,"GID")

DpSym = AnnotationDbi::select(org.Dp.eg.db, keys = GIDkeys, columns = c("GID","ENTREZID","SYMBOL","GENENAME"))
DpChr = na.omit(AnnotationDbi::select(org.Dp.eg.db, keys = GIDkeys, columns = c("GID","CHR")))
DpGO = na.omit(AnnotationDbi::select(org.Dp.eg.db, keys = GIDkeys, columns = c("GID","GO","EVIDENCE")))

AnnotationForge::makeOrgPackage(gene_info = DpSym, chromosome = DpChr, go = DpGO,
                                version="1.0",
                                maintainer="Wassim Salam <wassimsalam49@gmail.com>",
                                author="Wassim Salam <wassimsalam49@gmail.com>",
                                outputDir = ".",
                                tax_id="6669",
                                genus="Daphnia",
                                species="pulex",
                                goTable = "go")

install.packages("org.Dpulex.eg.db/", repos = NULL, type = "source")
library(org.Dpulex.eg.db)
############################
