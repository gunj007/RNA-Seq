# FlagsTat Alignment visualise 
#Provide alignment statistics and its visualization report, 
#including the percentage of aligned reads, 
#mapped reads, and potential issues with multi-mapping.


#if (!require(stringr)) install.packages("stringr", dependencies = TRUE)
#if (!require(dplyr)) install.packages("dplyr", dependencies = TRUE)
#if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)

# Load necessary libraries
library(dplyr)
library(stringr)
library(ggplot2)


file_path <- "7bamalignment_summary.txt"
alignment_data <- readLines(file_path)

samples <- c()
total_reads <- c()
mapped_reads <- c()
percentage_mapped <- c()


for (line in alignment_data) {
  
  # Detect total reads line for a new sample
  if (grepl("in total", line)) {
    # Extract sample ID and total reads
    sample_name <- str_extract(line, "\\d+ \\+ 0")  # Extracts the first number before "+ 0"
    total <- as.numeric(unlist(str_extract_all(line, "\\d+"))[1])
    
    # Add to the lists
    samples <- c(samples, sample_name)
    total_reads <- c(total_reads, total)
  }
  
  
  if (grepl("mapped", line) && grepl("N/A", line)) {
    mapped <- as.numeric(unlist(str_extract_all(line, "\\d+"))[1])
    mapped_reads <- c(mapped_reads, mapped)
    
    percent <- round((mapped / total) * 100, 2)
    percentage_mapped <- c(percentage_mapped, percent)
  }
}



  alignment_summary <- data.frame(
    Samplename = c("Heart_ZT0_1","Heart_ZT0_2","Heart_ZT12_1","Heart_ZT12_2","Liver_ZT0_1","Liver_ZT0_2","Liver_ZT12_1","Liver_ZT12_2"),
    Sample = samples,
    Total_Reads = total_reads,
    Mapped_Reads = mapped_reads,
    Percentage_Mapped = percentage_mapped  )
  
  # Display the summary
  print(alignment_summary)
  write.csv(alignment_summary, file = "1flagstatAlignmentSummary.csv")
  # Plot alignment statistics
  ggplot(alignment_summary, aes(x = Samplename, y = Percentage_Mapped)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    xlab("Sample") +
    ylab("Percentage Mapped (%)") +
    ggtitle("Alignment Statistics - Percentage of Mapped Reads per Sample") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

