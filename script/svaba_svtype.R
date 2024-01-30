# // modification from https://github.com/walaj/svaba/issues/4
#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
	stop("Rscript svaba_svtype.R input.vcf output.vcf", call.=FALSE)
}

# Read in svaba vcf file
cols <- colnames(read.table(pipe(paste0('grep -v "##" ', args[1], ' | sed s/#//')), header=T))
#cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))
cols[10] <- unlist(strsplit(cols[10], '\\.'))[10]
cols[11] <- unlist(strsplit(cols[11], '\\.'))[10]
svaba_uniq <- read.table(pipe(paste0('grep -v "##" ', args[1], ' | sed s/#//')), header=T, col.names = cols, stringsAsFactors = FALSE)

for ( i in 1:length(svaba_uniq[,1]) ) {
    # Find mate pair
    root <- gsub(":[12]", "", svaba_uniq[i,]$ID)
    mate1 <- paste0(root, ":1")
    mate2 <- paste0(root, ":2")
    alt1 <- svaba_uniq[svaba_uniq$ID == mate1, ]$ALT
    alt2 <- svaba_uniq[svaba_uniq$ID == mate2, ]$ALT
    chr1 = gsub(':[0-9]+(\\[|\\])[A-Z]*', "", alt1);	chr1 = gsub('[A-Z]*(\\[|\\])', "", chr1);
    chr2 = gsub(':[0-9]+(\\[|\\])[A-Z]*', "", alt2);	chr2 = gsub('[A-Z]*(\\[|\\])', "", chr2);

    # Determine sv type based on breakpoint orientation
    sv_type = ""
    if ( chr1 != chr2 ) {
	sv_type <- "=BND"
    } else {
    	if ((grepl("\\[", alt1) & grepl("\\[", alt2)) | (grepl("\\]", alt1) & grepl("\\]", alt2))) {
        	sv_type <- "=INV"        
    	} else if (grepl("[A-Z]+\\[", alt1) & grepl("^\\]", alt2)) {
        	sv_type <- "=DEL"
    	} else if (grepl("^\\]", alt1) & grepl("[A-Z]+\\[", alt2)) {
        	sv_type <- "=DUP"
    	} else { 
        	sv_type <- "=UNKNOWN"
   	}
    }
    svaba_uniq[i,]$INFO = gsub('=BND', sv_type, svaba_uniq[i,]$INFO)
}

write.table(svaba_uniq, file=args[2], quote=F, row.names=F, col.names=F, sep='\t')

