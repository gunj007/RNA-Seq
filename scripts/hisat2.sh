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

# Create output dir for Hisat2
sam_dir="${input_dir}/sam"
bam_dir="${input_dir}/bam"

mkdir -p "$bam_dir"
mkdir -p "$sam_dir"

# Step 1: Count .gz files in the input directory (inside fastp folder)
fastp_dir="${input_dir}/fastp"
if [ ! -d "$fastp_dir" ]; then
  echo "Directory $fastp_dir not found. Please check the path."
  exit 1
fi

# Count .gz files
gz_count=$(ls "$fastp_dir"/*.gz 2>/dev/null | wc -l)
if [ "$gz_count" -eq 0 ]; then
  echo "No .gz files found in the folder: $fastp_dir"
  exit 1
else
  echo "Number of .gz files in $fastp_dir: $gz_count"
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
echo "Running HISAT2 for all paired-end reads..."

for r1 in "$fastp_dir"/*R1.fastq.gz; do
  # Extract base name (without _R1.fastq.gz)
  base_name=$(basename "$r1" R1.fastq.gz)
  
  # Corresponding R2 file
  r2="${fastp_dir}/${base_name}R2.fastq.gz"
  
  # Output SAM and BAM filenames
  output_sam="${sam_dir}/${base_name}.sam"
  output_bam="${input_dir}/bam/${base_name}.bam"
  output_sum="${input_dir}/sam/${base_name}.txt"

  # Check if R2 file exists
  if [ -f "$r2" ]; then
    echo "Processing ${base_name}..."
    
    # Run hisat2 alignment you can add -p 10 threads crashing ERR137OOM 
    hisat2 -p 4 -x "${genome_dir}/genome/genome" -1 "$r1" -2 "$r2" -S "$output_sam" --summary-file $output_sum --time --quiet

    echo "HISAT2 alignment completed for ${base_name}: $output_sam"
    
    # Convert SAM to BAM
    samtools view -bS "$output_sam" > "$output_bam"
    
    # Remove SAM file
    rm "$output_sam"
    echo "BAM file created and SAM file removed: $output_bam"
  else
    echo "R2 file not found for ${base_name}, skipping..."
  fi
done

echo "All alignments completed successfully!"

