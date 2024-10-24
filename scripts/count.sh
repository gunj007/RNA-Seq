#!/bin/bash
# *NOTE: Use nohup while running the script
# Check if the required directories are provided

if [ $# -ne 3 ]; then
  echo "Usage: $0 <input_directory> <genome_directory> <gtf_directory>"
  exit 1
fi

# Set input arguments
input_dir="$1"
genome_dir="$2"
gtf_dir="$3"

# Step 1: Run qc.sh for quality control
echo "Running quality control (qc.sh)..."
bash ../biostateai/scripts/qc.sh "$input_dir"

# Check if qc.sh ran successfully
if [ $? -ne 0 ]; then
  echo "QC step failed. Exiting..."
  exit 1
fi

echo "######### QC step completed successfully #########"
echo -e "\n"
# Step 2: Run hisat.sh for alignment
echo "Running alignment (hisat2.sh)..."
bash ../biostateai/scripts/hisat2.sh "$input_dir" "$genome_dir" "$gtf_dir"

# Check if hisat.sh ran successfully
if [ $? -ne 0 ]; then
  echo "HISAT2 alignment step failed. Exiting..."
  exit 1
fi

echo "######### Alignment step completed successfully #########"
echo -e "\n"
# Step 3: Perform feature counting with featureCounts
echo "Performing feature counting with featureCounts..."

# Prepare the list of BAM files to count
bam_files=("$input_dir/bam"/*.bam)

# Check if there are any BAM files
bam_count=${#bam_files[@]}
if [ "$bam_count" -eq 0 ]; then
  echo "No BAM files found in $input_dir/bam. Exiting..."
  exit 1
else
  echo "Number of BAM files found: $bam_count"
fi

# Run featureCounts
featureCounts -p -t gene --extraAttributes gene_name,gene_type --primary -a "$gtf_dir" -o "$input_dir/featurecounts.txt" "${bam_files[@]}"

# Check if featureCounts ran successfully
if [ $? -ne 0 ]; then
  echo "featureCounts step failed. Exiting..."
  exit 1
fi

echo "Feature counting completed successfully."
echo -e "\n"
# Step 4: Process counts file
echo "Processing counts file..."
sed '1d' "$input_dir/bam/featurecounts.txt" > "$input_dir/allfeaturecounts.tsv"

# Check if sed command was successful
if [ $? -ne 0 ]; then
  echo "Error processing counts file. Exiting..."
  exit 1
fi

echo "Counts file processed successfully: counts.tsv"
echo -e "\n"

# Display the completion message
center_text "***RNA-Seq Pipeline Analysis Completed!***"
