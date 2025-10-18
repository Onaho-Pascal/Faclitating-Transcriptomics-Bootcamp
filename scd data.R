# install.packages("BiocManager")
install.packages("edgeR")
install.packages("ggplot2")



library(edgeR)
library(ggplot2)


scd_data <- read.delim("C:/Users/user/Documents/Bioinformatics Projects/Sickle Cell Project/Sickle-Cell-Research/scd_raw_data.txt", row.names =1, header =TRUE)

#setwd("C:/Users/user/Documents/Bioinformatics Projects/Sickle Cell Project/Sickle-Cell-Research")
#rownames(counts) <- counts[, 1]   # assign first column as row names
#counts <- counts[, -1]            # remove that first column

BiocManager::install("edgeR")
library(edgeR)
boxplot(scd_data,
        las = 2,
        col = "red",
        main = "Raw Count Data (Before log normalization)",
        ylab = "Raw Counts")

log_scd_data <- log2(scd_data + 1)

boxplot(log_scd_data,
        las = 2,
        col = "red",
        main = "Log Transformed Data (After log normalization)",
        ylab = "Log Counts")

dge <- DGEList(counts = log_scd_data)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)


boxplot(logCPM,
        las = 2,
        col = "red",
        main = "CPM Transformed Data (After log normalization)",
        ylab = "CPM Counts")
