#!/bin/bash

# Check if input directory is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_directory>"
  exit 1
fi

# Set input directory
input_dir="$1"
bin_dir="${input_dir}/bin"
# Create the qc reports directory
qc_dir="${input_dir}/qcreports"
fastp_dir="${input_dir}/fastp"
fastprepo="${fastp_dir}/fastprepo"
fastplog="${bin_dir}/fastplog3"
mkdir -p "$fastplog"
mkdir -p "$qc_dir"
mkdir -p "$fastp_dir"
mkdir -p "$bin_dir"
mkdir -p "$fastprepo"
# Step 1: Run FastQC on all fastq.gz files
echo "Running FastQC..."
fastqc -o "$qc_dir" "$input_dir"/*.fastq.gz --threads 3 > ${bin_dir}/S1_fastqc1.log 2>&1
echo "Completed FastQC..."

# Step 2: Run MultiQC to aggregate FastQC reports
echo "Running MultiQC..."
multiqc -o "$qc_dir" "$qc_dir"/*.zip > ${bin_dir}/S2_multiqc1.log 2>&1
echo "Completed MultiQC..."

echo "QC1 process completed!"
echo -e "\n"

# Step 3: Run fastp for paired-end reads (adjust file names as needed)
echo "Running fastp..."

for file in "$input_dir"/*_R1.fastq.gz
do
  # Extract the base name
  base_name=$(basename "$file" _R1.fastq.gz)
  
  # Set the input and output filenames for R1 and R2
  r1="$input_dir/${base_name}_R1.fastq.gz"
  r2="$input_dir/${base_name}_R2.fastq.gz"
  out_r1="$fastp_dir/${base_name}R1.fastq.gz"
  out_r2="$fastp_dir/${base_name}R2.fastq.gz"
  
  # Define output for fastp reports (HTML and JSON)
  html_report="$fastprepo/${base_name}_fastp.html"
  json_report="$fastprepo/${base_name}_fastp.json"
  
  # Define error log specific for this sample
  log_file="$fastplog/${base_name}_fastp_error.log"

  # Run fastp
  fastp -i "$r1" -I "$r2" -o "$out_r1" -O "$out_r2" --detect_adapter_for_pe --thread 8 -h "$html_report" -j "$json_report" > "$log_file" 2>&1 

  # Notify that this sample is done (optional). will make dir for logs
  echo "Processed $base_name "
done
echo "Completed trimming..."
echo -e "\n"

qc2_dir="${fastp_dir}/qcreportstrim"
mkdir -p "$qc2_dir"


# Step 4: Re-run FastQC on all fastq.gz files which are trimmed rm -t as OOM
echo "Running FastQC for trimmed fastq files..."
fastqc -o "$qc2_dir" "$input_dir"/fastp/*.fastq.gz -t 3 > ${bin_dir}/S4_fastqc2.log 2>&1
echo "Completed FastQC..."

# Step 5: Run MultiQC to aggregate Trimmed FastQC reports
echo "Running MultiQC..."
multiqc -o "$qc2_dir" "$qc2_dir"/*.zip > ${bin_dir}/S5_multiqc2.log 2>&1
echo "Completed MultiQC..."

echo "QC2 process completed!"
echo -e "\n"
