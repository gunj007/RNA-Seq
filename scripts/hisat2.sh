#!/bin/bash

# Check if the required directories are provided
if [ $# -ne 3 ]; then
  echo "Usage: $0 <input_directory> <genome_directory> <gtf_directory>"
  exit 1
fi

# Set input arguments
input_dir="$1"
genome_dir="$2"
gtf_dir="$3"
bin_dir="${input_dir}/bin"
# Create output dir for Hisat2
sam_dir="${input_dir}/sam"
bam_dir="${input_dir}/bam"
samsum="${bin_dir}/samsum6"
sammet="${bin_dir}/sammetrics6"
bamstat="${bin_dir}/bamstat7"

mkdir -p "$bam_dir" "$sam_dir" "$bamstat" "$samsum" "$sammet"

echo "Bin-folders created..."
# Step 1: Count .gz files in the input directory (inside fastp folder)
fastp_dir="${input_dir}/fastp"
if [ ! -d "$fastp_dir" ]; then
  echo "Directory $fastp_dir not found. Please check the path."
  exit 1
fi

# Count .gz files
gz_count=$(ls "$fastp_dir"/*.gz 2>/dev/null | wc -l)
if [ "$gz_count" -eq 0 ]; then
  echo "Total .gz files in the folder: $fastp_dir"
  exit 1
else
  echo "Number of trimmed files in $fastp_dir: $gz_count"
fi

# Step 2: Check if genome folder exists
# If genome is not built, build it using: 1st download genome.fa.gz file
# 2nd gunzip it "gunzip mgiGenome/GRCm39.primary_assembly.genome.fa.gz"
# 3rd build genome "hisat2-build mgiGenome/GRCm39.primary_assembly.genome.fa genome"

if [ -d "$genome_dir/genome" ]; then
  echo "Genome folder available: genome"
else
  echo "Genome folder not found in $genome_dir"
  exit 1
fi


# Step 3: Perform hisat2 alignment for each paired-end read in fastp folder
echo "Running HISAT2..."

for r1 in "$fastp_dir"/*R1.fastq.gz; do
  # Extract base name (without _R1.fastq.gz)
  base_name=$(basename "$r1" R1.fastq.gz)
  
  # Corresponding R2 file
  r2="${fastp_dir}/${base_name}R2.fastq.gz"
  
  # Output SAM and BAM filenames
  output_sam="${sam_dir}/${base_name}.sam"
  output_bam="${bam_dir}/${base_name}.bam"
  sorted_bam="${bam_dir}/${base_name}_sorted.bam"
  output_sum="${samsum}/${base_name}_alignmentSummary.txt"
  output_flag="${bamstat}/${base_name}_flagstat.txt"
  output_met="${sammet}/${base_name}_alignmentmetrics.txt"

  # Check if R2 file exists
  if [ -f "$r2" ]; then
    echo "Processing ${base_name}..."
    
    # Run hisat2 alignment you can add -p 10 threads crashing ERR137OOM 
    hisat2 -p 3 -x "${genome_dir}/genome/genome" -1 "$r1" -2 "$r2" -S "$output_sam" --time --summary-file "${output_sum}" --met-file "${output_met}" > ${bin_dir}/S6_hisat2.log 2>&1
    
    echo "SAM Summary file: $output_sum"
    echo "Alignment completed: ${base_name}"
    
    # Convert SAM to BAM and then sort Bam
    samtools view -bS "$output_sam" > "$output_bam"
    samtools sort "$output_bam" -o "$sorted_bam" 
    echo "BAM file: $sorted_bam"
    # Remove SAM file
    du -h "$output_sam" && rm "$output_sam"
    du -h "$output_bam" && rm "$output_bam"
    echo "SAM file removed"
    echo -e  "\n"
    
  else
    echo "R2 file not found for ${base_name}, skipping..."
  fi
done

echo "BAM flagsat..."
	for bamFile in "$bam_dir"/*d.bam; do
	    base_name=$(basename "$bamFile" .bam)
	    echo "Flagstat processes: ${base_name}."

	    samtools flagstat "$bamFile" >> "$output_flag" &
	done

	# Wait for all background processes to complete
	wait


echo "flagstat completed..."  
echo -e  "\n"

# Filter for total and mapped lines from all text files and save to alignment_summary.txt
grep -E "total|mapped" ${bamstat}/*.txt > ${bin_dir}/7bamalignment_summary.txt

echo " Aligmnment completed!"
echo -e  "\n"

