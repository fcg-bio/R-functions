## Original file = /scicore/home/ivanek/GROUP/ivanek/Rscripts/runSTAR.R
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @return 
##' @author Robert Ivanek
##' @param genome 
##' @param projDir 
##' @param storeDir 
##' @param chromOrder 
##' @param sjdbGTFfile 
##' @param sjdbOverhang 
##' @param sjdbFileChrStartEnd 
##' @param sjdbGTFchrPrefix 
##' @param sjdbGTFfeatureExon 
##' @param sjdbGTFtagExonParentTranscript 
##' @param sjdbGTFtagExonParentGene 
##' @param sjdbInsertSave 
##' @param genomeSAindexNbases 
##' @param genomeChrBinNbits 
##' @param genomeSAsparseD 
generateSTARindex <- function(genome,
                              projDir=".",
                              storeDir=NULL,
                              chromOrder=NULL,
                              sjdbGTFfile=NULL,
                              sjdbOverhang=100,
                              sjdbFileChrStartEnd=NULL,
                              sjdbGTFchrPrefix=NULL,
                              sjdbGTFfeatureExon="exon",
                              sjdbGTFtagExonParentTranscript="transcript_id",
                              sjdbGTFtagExonParentGene="gene_id",
                              sjdbInsertSave="Basic",
                              genomeSAindexNbases=14,
                              genomeChrBinNbits=18,
                              genomeSAsparseD=1) {
    ## 
    indDir <- sprintf("%s/genome_index", projDir)
    if(!dir.exists(indDir)) dir.create(indDir, recursive=TRUE)
    ## 
    if (!is.null(storeDir) & file.exists(sprintf("%s/genomeParameters.txt", storeDir))) {   # missing proper test for STAR index
        ## rsync files to index dir (/scratch)
        system(sprintf("rsync -a %s/* %s/", storeDir, indDir))
    } else {
        genDir <- sprintf("%s/genome", projDir)
        if(!dir.exists(genDir)) dir.create(genDir, recursive=TRUE)
        ## create directory for FASTA sequences, download chromFa.tar.gz from UCSC
        if (!file.exists(sprintf("%s/md5sum.txt", genDir))) {
            download.file(url=sprintf("http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/md5sum.txt", genome),
                          destfile=sprintf("%s/md5sum.txt", genDir))
        }
        if (!file.exists(sprintf("%s/chromFa.tar.gz", genDir))) {
            download.file(url=sprintf("http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/chromFa.tar.gz", genome),
                          destfile=sprintf("%s/chromFa.tar.gz", genDir))
        }
        ##
        md5old <- read.table(file=sprintf("%s/md5sum.txt", genDir), row.names=2, stringsAsFactors=FALSE)
        md5old <- md5old["chromFa.tar.gz",1]
        md5new <- tools::md5sum(sprintf("%s/chromFa.tar.gz", genDir))
        if (md5new != md5old)
            stop(sprintf("md5sum of downloaded %s/chromFa.tar.gz file does not match the original one", genDir))
        ## extract FASTA
        untar(tarfile=sprintf("%s/chromFa.tar.gz", genDir), exdir=genDir, compressed="gzip", extras="--overwrite")
        ## run STAR genomeGenerate
        faFiles <- list.files(genDir, full.names=TRUE, pattern="\\.fa$")
        if (!is.null(chromOrder)) {
            names(faFiles) <- sub("\\.fas?t?a?", "", basename(faFiles))
            faFiles <- faFiles[chromOrder]
        }
        faFiles <- paste(faFiles, collapse=" ")
        ## parse the STAR parameters
        addParams <- list(
            genomeSAindexNbases = genomeSAindexNbases,
            genomeChrBinNbits = genomeChrBinNbits,
            genomeSAsparseD = genomeSAsparseD,
            sjdbOverhang = sjdbOverhang,
            sjdbFileChrStartEnd = sjdbFileChrStartEnd,
            sjdbGTFfile = sjdbGTFfile,
            sjdbGTFchrPrefix = sjdbGTFchrPrefix,
            sjdbGTFfeatureExon = sjdbGTFfeatureExon,
            sjdbGTFtagExonParentTranscript = sjdbGTFtagExonParentTranscript,
            sjdbGTFtagExonParentGene = sjdbGTFtagExonParentGene,
            sjdbInsertSave = sjdbInsertSave
        )
        addParams <- lapply(addParams, function(x) if (is.null(x)) "-" else x)
        addParams <- paste(paste("--", names(addParams)," ", unlist(addParams), sep=""), collapse=" ")
        addParams <- ""
        ##
        cmd <- sprintf("STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s %s",
                       indDir, faFiles, addParams)
        system(cmd)
        if (!is.null(storeDir)) {
            ## rsync index files back to storeDir (home)
            if(!dir.exists(storeDir)) dir.create(storeDir, recursive=TRUE)
            system(sprintf("rsync -a %s/* %s/", indDir, storeDir))
        }
    }
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param sampleFile 
##' @param genome 
##' @param projDir 
##' @param runThreadN 
##' @param bamStoreDir 
##' @param outFilterType 
##' @param outFilterMultimapNmax 
##' @param outSAMunmapped 
##' @param outSAMtype 
##' @param limitBAMsortRAM 
##' @param addParams 
##' @return 
##' @author Robert Ivanek
runSTARmapping1stPass <- function(sampleFile,
                                  genome,
                                  projDir=".",
                                  runThreadN=1,
                                  bamStoreDir=NULL,
                                  outFilterType="BySJout",
                                  outFilterMultimapNmax=1,
                                  outSAMunmapped="Within",
                                  outSAMtype="BAM SortedByCoordinate",
                                  limitBAMsortRAM=5e10,
                                  addParams="") {
    
    if(!is.null(bamStoreDir) && !dir.exists(bamStoreDir))
        dir.create(bamStoreDir, recursive=TRUE)
    ## 
    bamDir <- sprintf("%s/common_bam", projDir)
    if(!dir.exists(bamDir)) dir.create(bamDir, recursive=TRUE)
    ## 
    indDir <- sprintf("%s/genome_index", projDir)
    if(!dir.exists(indDir))
        stop("Missing genome index, please run the function \"generateSTARindex\" !")
    ## 
    proj1 <- QuasR:::createQProject(
        sampleFile=sampleFile,
        genome=grep(sprintf("%s$", genome), BSgenome::installed.genomes(), value=T),
        paired=NULL,
        splicedAlignment=TRUE,
        alignmentParameter=NULL,
        maxHits=1,
        aligner="Rbowtie",
        snpFile=NULL,
        auxiliaryFile=NULL,
        bisulfite="no",
        projectName="QProject",
        cacheDir="/scratch",
        lib.loc=installed.packages()[grep(sprintf("%s$", genome), BSgenome::installed.genomes(), value=T),"LibPath"],
        alignmentsDir=bamDir)
    ## load index in shared memory if there is more than one FASTQ file (or one pair of files)
    if (length(proj1) > 1) {
        system(sprintf("STAR --runMode alignReads --genomeDir %s --genomeLoad LoadAndExit",  indDir))
        on.exit(system(sprintf("STAR --runMode alignReads --genomeDir %s --genomeLoad Remove",  indDir)))
        genomeLoad <- "LoadAndKeep"
    } else {
        genomeLoad <- "NoSharedMemory"
    }
    ## setting parameters
    ## are the files compressed? assuming all files are the same
    if (length(grep("\\.gz$", proj1@reads[1,1]))) {
        readFilesCommand <- "zcat"
    } else if (length(grep("\\.bz$", proj1@reads[1,1]))) {
        readFilesCommand <- "bzcat"
    } else {
        readFilesCommand <- "-"
    }
    ## mapping, currently using for loop
    for (i in 1:length(proj1)) {
        ## fastq file names
        if (proj1@paired == "no") {
            readFilesIn <-  proj1@reads[i,1]
        } else if (proj1@paired == "fr") {
            readFilesIn <-  paste(proj1@reads[i,1], proj1@reads[i,2], sep=" ")
        } else {
            stop("Only paired end mode \"fr\" is currently suported!")
        }
        ## combine directory and file prefix
        outFileNamePrefix <- basename(tools::file_path_sans_ext(proj1@reads[i,1], compression=TRUE))
        outFileNamePrefix <- sprintf("%s/%s_", bamDir, outFileNamePrefix)
        ## command
        cmd <- sprintf("STAR --runMode alignReads --genomeDir %s --genomeLoad %s --runThreadN %s",
                       indDir, genomeLoad, runThreadN)
        cmd <- sprintf("%s --readFilesIn %s --readFilesCommand %s", cmd, readFilesIn , readFilesCommand)
        cmd <- sprintf("%s --outSAMattrRGline ID:%s SM:%s PL:ILLUMINA", cmd, proj1@reads$SampleName[i], proj1@reads$SampleName[i])
        cmd <- sprintf("%s --outFileNamePrefix %s", cmd, outFileNamePrefix)
        cmd <- sprintf("%s --outFilterType %s --outFilterMultimapNmax %s", cmd, outFilterType, outFilterMultimapNmax)
        cmd <- sprintf("%s --outSAMtype %s --outSAMunmapped  %s ", cmd, outSAMtype, outSAMunmapped)
        cmd <- sprintf("%s --limitBAMsortRAM %.0f", cmd, limitBAMsortRAM)
        cmd <- sprintf("%s %s", cmd, addParams)
        system(cmd)
        proj1@alignments[i, 1] <- sprintf("%sAligned.sortedByCoord.out.bam", outFileNamePrefix)
    }
    ## index files 
    system(sprintf("parallel -j %s samtools index {} ::: %s", runThreadN,
           paste(proj1@alignments[,1], collapse=" ")))
    
    ## rsync bam files back to bamStoreDir (home)
    if (!is.null(bamStoreDir)) {
        if (!dir.exists(bamStoreDir))
            dir.create(bamStoreDir, recursive=TRUE)
        system(sprintf("rsync -a %s/* %s/", bamDir, bamStoreDir))
        ## change the path to bam files
        proj1@alignmentsDir <- bamStoreDir
        proj1@alignments[,1] <- sprintf("%s/%s", bamStoreDir, basename(proj1@alignments[,1]))
    }
    return(proj1)
}

