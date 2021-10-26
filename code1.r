# Creating and setting up working directory
dir.home <- Sys.getenv("HOME")
dir.results <- file.path(dir.home, "Documents", "cmb", "stat1", "tp1")
dir.create(path = dir.results, showWarnings = FALSE, recursive = TRUE)
setwd(dir.results)

# Downloading and extracting required files
library(R.utils)
download.file("http://ftp.ensembl.org/pub/release-104/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.104.gtf.gz", "genome.gtf.gz")
genome <- gunzip("genome.gtf.gz")
download.file("https://raw.githubusercontent.com/gsgautham98/probastat1/main/chrom_sizes.tsv", "chrom_sizes.tsv")
download.file("https://raw.githubusercontent.com/gsgautham98/probastat1/main/codon_frequencies.tab", "codon_frequencies.tab")

# Importing main dataset
feature.table <- read.table("genome.gtf", comment.char = "#", sep = "\t", header = FALSE, row.names = NULL)
names(feature.table) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Computing coding genes length
feature.table$length <- feature.table$end - feature.table$start
knitr::kable(table(feature.table$feature), col.names = c("Feature", "Freq"))
cds <- subset(feature.table, feature == "CDS")
cds.count <- table(cds$seqname)
knitr::kable(cds.count, col.names = c("Chromosome", "CDS Freq"))
chromosomes <- read.delim("chrom_sizes.tsv", sep = "\t")
chromosomes <- chromosomes[order(chromosomes$chrom),]
genes <- subset(feature.table, feature == "gene")
genes.count <- table(genes$seqname)
knitr::kable(genes.count, col.names = c("Chromosome", "Gene Freq"))
genes.per.mb <- genes.count / (chromosomes$size * 0.00001)
barplot(genes.per.mb, xlab = "Chromosome", ylab = "Gene density")

# Histogram of CDS lengths
cds.length.hist <- hist(cds$length, breaks = 74, xlab = "CDS length", main = "Histogram of CDS length")
print(cds.length.hist)

# Polygon of frequency of genes of median classes
genes.length.break <- table(findInterval(genes$length, cds.length.hist$breaks))
plot(genes.length.break, type = "l", xlab = "Class index", ylab = "Frequency", main = "Frequency of genes of each median class")

# CDS length visualisation using a boxplot
boxplot(cds$length ~ cds$seqname, xlab = "Chromosome", ylab = "CDS lengths", main = "Boxplot of CDS length in each chromosome")

# Function to find the mode
mode.finder <- function(a) {
  a.unique <- unique(a)
  a.unique[which.max(tabulate(match(a, a.unique)))]
}

# Calculating descriptive parameters
chrom3 <- subset(genes, seqname == "III")
chrom3.mean <- mean(chrom3$length)
chrom3.median <- median(chrom3$length)
chrom3.mode <- mode.finder(chrom3$length)
sprintf("For chromosome III gene lengths, mean = %f, median = %f, mode = %f", chrom3.mean, chrom3.median, chrom3.mode)
chrom3.var <- var(chrom3$length)
chrom3.sd <- sd(chrom3$length) * sqrt((length(chrom3$length) -1) / length(chrom3$length))
chrom3.iqr <- IQR(chrom3$length)
sprintf("Variance = %f, std deviation = %f, interquartile range = %f", chrom3.var, chrom3.sd, chrom3.iqr)
genes.mean <- mean(genes$length)
genes.median <- median(genes$length)
genes.mode <- mode.finder(genes$length)
sprintf("For all gene lengths, mean = %f, median = %f, mode = %f", genes.mean, genes.median, genes.mode)
genes.var <- var(genes$length)
genes.sd <- sd(genes$length)
genes.iqr <- IQR(genes$length)
sprintf("Variance = %f, std deviation = %f, interquartile range = %f", genes.var, genes.sd, genes.iqr)

# Visualising some descriptive parameters
hist(cds$length, breaks = 74, xlab = "CDS length", main = "Histogram of CDS length")
arrows(genes.median, 10, genes.median - (genes.iqr / 2), 10, 0.05, col = "red")
arrows(genes.median, 10, genes.median + (genes.iqr / 2), 10, 0.05, col = "red")
arrows(genes.mean, 0, genes.mean, 800, 0.005, col = "blue")

# Confidence interval computation
n <- length(chrom3$length)
margin <- qt(0.975, df=n-1) * chrom3.sd / sqrt(n)
lower.end <- chrom3.mean - margin
upper.end <- chrom3.mean + margin
sprintf("The confidence interval for chromosome III around the mean is %f to %f", lower.end, upper.end)

# Empirical distribution function visualisation
lengths.mids <- data.frame(cds.length.hist$mids, cds.length.hist$counts)
names(lengths.mids) <- c("mids", "counts")
relative.freq <- lengths.mids$counts / sum(lengths.mids$counts)
lengths.mids <- cbind(lengths.mids, relative.freq)
lengths.mids <- cbind(lengths.mids, cumsum(lengths.mids$counts))
empirical.freq <- ecdf(lengths.mids$counts)
plot(lengths.mids$mids, lengths.mids$counts, type = "h", xlab = "Classes", ylab = "Frequency", main = "Frequency of each median class")
plot(empirical.freq, xlab = "Median class size", ylab = "Relative frequency", main = "Empirical cumulative distribution")


# Computing and visualising expected frequency distribution of different genes based on length
prob.3nt <- c(0.0182536851094, 0.0223781252573, 0.0128875284527, 0.0201230463324)
gene.lengths <- c()
freq.genes <- c()
for (gene.length in 1:200) {
  freq.gene <- 12156679 * (prob.3nt[1] * (1 - sum(prob.3nt[2:4])) ^ gene.length)
  freq.genes <- append(freq.genes, freq.gene)
  gene.lengths <- append(gene.lengths, gene.length)
}
plot(gene.lengths, freq.genes, type = "l", xlab = "Gene length", ylab = "Frequency")

