N <- function(n) prettyNum(n, big.mark = ",")

# top dir for the project
topDir <- "../glycotreg/"

# Public HTML for file downloads
public_html <- "http://www.compbio.dundee.ac.uk/user/mgierlinski/glycotreg/"

makeDirs <- function(topDir) {
  lapply(subDirs, function(d) paste0(topDir, d))
}

# Sub-directories
subDirs <- list(
  top = "",
  fastq = "fastq/",
  qc = "qc/",
  multiqc = "multiqc/",
  genome = "genome/",
  starmap = "starmap/",
  bam = "bam/",
  bedgraph = "bedgraph/",
  readcount = "readcount/",
  salmon = "salmon/",
  download = "download/",
  data = "data/"
)

# All directories
dirs <- makeDirs(topDir)

genomeFile <- paste0(dirs$genome, "Mus_musculus.GRCm38.dna_rm.primary_assembly.fa")
transFile <- paste0(dirs$genome, "Mus_musculus.GRCm38.cds.all.fa")
gtfFile <- paste0(dirs$genome, "Mus_musculus.GRCm38.93.gtf")

