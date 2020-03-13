
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('filtered/SRS015072_R1.fastq.gz', 'filtered/SRS015072_R2.fastq.gz', "R1" = c('filtered/SRS015072_R1.fastq.gz'), "R2" = c('filtered/SRS015072_R2.fastq.gz')),
    output = list('model/ErrorRates_R1.rds', 'model/ErrorRates_R2.rds', 'figures/ErrorRates_R1.pdf', 'figures/ErrorRates_R2.pdf', "errR1" = 'model/ErrorRates_R1.rds', "errR2" = 'model/ErrorRates_R2.rds', "plotErr1" = 'figures/ErrorRates_R1.pdf', "plotErr2" = 'figures/ErrorRates_R2.pdf'),
    params = list(),
    wildcards = list(),
    threads = 4,
    log = list('logs/dada2/learnErrorRates.txt'),
    resources = list(),
    config = list("java_mem" = 30, "threads" = 4, "truncLen" = c(0, 0), "maxEE" = c(2, 2), "truncQ" = 2, "learn_nbases" = '100e6', "chimera_method" = 'consensus', "max_length_variation" = 7, "idtaxa_dbs" = list("Silva" = '/mnt/d/github/amplicon-seq-dada2/silva_138/SILVA_SSU_r138_2019.RData', "GTDB" = '/mnt/d/github/amplicon-seq-dada2/GTDB/GTDB_r89-mod_June2019.RData', "RDP" = 'http://www2.decipher.codes/Classification/TrainingSets/RDP_v16_March2018.RData')),
    rule = 'learnErrorRates'
)
######## Original script #########
library(dada2)
library(ggplot2)
sink(snakemake@log[[1]])


errF <- learnErrors(snakemake@input[['R1']], nbases=snakemake@config[["learn_nbases"]], multithread=snakemake@threads,randomize = TRUE)
errR <- learnErrors(snakemake@input[['R2']], nbases=snakemake@config[["learn_nbases"]], multithread=snakemake@threads,randomize = TRUE)

save(errF,file=snakemake@output[['errR1']])
save(errR,file=snakemake@output[['errR2']])



## ---- plot-rates ----
plotErrors(errF,nominalQ=TRUE)
ggsave(snakemake@output[['plotErr1']])
plotErrors(errR,nominalQ=TRUE)
ggsave(snakemake@output[['plotErr2']])




