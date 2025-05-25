# MLST
Eimeria maxima typing using MLST. The gene sequenced were dnaj, rcc, and hsp.

#activate qiime2 envt

module load qiime2/2024.10
source activate qiime2-amplicon-2024.10

#copy necessary folders/files

#create manifest files- for forward and reverse reads separately (SE) and for PE.
# nano manifest.py

import os

# Define the base directory (hsp folder)
base_dir = "/scratch/aubsxf001/data/qiime/hsp"

# Define the output directory (qiime folder)
output_dir = "/scratch/aubsxf001/data/qiime"

# Initialize manifest content for forward, reverse, and paired-end reads
forward_manifest_content = "sample-id\tabsolute-filepath\n"  # Corrected header for single-end forward reads
reverse_manifest_content = "sample-id\tabsolute-filepath\n"  # Corrected header for single-end reverse reads
pe_manifest_content = "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n"  # Header for paired-end reads

# Iterate through each folder in the base directory
for folder in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder)
    
    # Check if it's a directory
    if os.path.isdir(folder_path):
        # Find all FASTQ files in the folder
        fastq_files = [f for f in os.listdir(folder_path) if f.endswith('.fastq.gz')]
        
        # Separate R1 and R2 files
        r1_files = [f for f in fastq_files if '_R1_' in f]
        r2_files = [f for f in fastq_files if '_R2_' in f]
        
        # Ensure there is exactly one R1 and one R2 file
        if len(r1_files) == 1 and len(r2_files) == 1:
            r1_path = os.path.join(folder_path, r1_files[0])
            r2_path = os.path.join(folder_path, r2_files[0])
            
            # Extract the sample ID from the file name
            sample_id = r1_files[0].split('_')[0]  # Extract the sample ID (e.g., RH052)
            
            # Add to forward manifest
            forward_manifest_content += f"{sample_id}\t{r1_path}\n"
            
            # Add to reverse manifest
            reverse_manifest_content += f"{sample_id}\t{r2_path}\n"
            
            # Add to paired-end manifest
            pe_manifest_content += f"{sample_id}\t{r1_path}\t{r2_path}\n"

# Write the forward manifest file
forward_manifest_path = os.path.join(output_dir, "forward_manifest.csv")
with open(forward_manifest_path, 'w') as f:
    f.write(forward_manifest_content)

# Write the reverse manifest file
reverse_manifest_path = os.path.join(output_dir, "reverse_manifest.csv")
with open(reverse_manifest_path, 'w') as f:
    f.write(reverse_manifest_content)

# Write the paired-end manifest file
pe_manifest_path = os.path.join(output_dir, "paired_end_manifest.csv")
with open(pe_manifest_path, 'w') as f:
    f.write(pe_manifest_content)

print(f"Forward manifest file created at: {forward_manifest_path}")
print(f"Reverse manifest file created at: {reverse_manifest_path}")
print(f"Paired-end manifest file created at: {pe_manifest_path}")


#run the script
python manifest_hsp.py

#Removing the primer sequences (20 bp length) and denoise data

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs hsp-paired-end-demux.qza \
  --p-trim-left-f 20 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table table-hsp-paired-end.qza \
  --o-representative-sequences rep-seqs-hsp-paired-end.qza \
  --o-denoising-stats stats-hsp-paired-end.qza

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-hsp-forward.qza \
  --p-trim-left 20 \
  --p-trunc-len 0 \
  --o-table table-hsp-forward.qza \
  --o-representative-sequences rep-seqs-hsp-forward.qza \
  --o-denoising-stats stats-hsp-forward.qza

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-hsp-reverse.qza \
  --p-trim-left 20 \
  --p-trunc-len 0 \
  --o-table table-hsp-reverse.qza \
  --o-representative-sequences rep-seqs-hsp-reverse.qza \
  --o-denoising-stats stats-hsp-reverse.qza# Visualize the feature table

#Generate visualization

qiime feature-table summarize \
  --i-table table-hsp-forward.qza \
  --o-visualization table-hsp-forward.qzv

qiime metadata tabulate \
  --m-input-file stats-hsp-forward.qza \
  --o-visualization stats-hsp-forward.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-hsp-forward.qza \
  --o-visualization rep-seqs-hsp-forward.qzv

qiime feature-table summarize \
  --i-table table-hsp-reverse.qza \
  --o-visualization table-hsp-reverse.qzv

qiime metadata tabulate \
  --m-input-file stats-hsp-reverse.qza \
  --o-visualization stats-hsp-reverse.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-hsp-reverse.qza \
  --o-visualization rep-seqs-hsp-reverse.qzv

qiime feature-table summarize \
  --i-table table-hsp-paired-end.qza \
  --o-visualization table-hsp-paired-end.qzv

qiime metadata tabulate \
  --m-input-file stats-hsp-paired-end.qza \
  --o-visualization stats-hsp-paired-end.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-hsp-paired-end.qza \
  --o-visualization rep-seqs-hsp-paired-end.qzv

#De-novo clustering of the rep-seqs at 99% identity (because we set the OTU percent identity to 97%, most clusters consist of only one sequence).

qiime vsearch cluster-features-de-novo \
  --i-table table-hsp-forward.qza \
  --i-sequences rep-seqs-hsp-forward.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-hsp-forward-99.qza \
  --o-clustered-sequences rep-seqs-hsp-forward-99.qza

#Generate visualization of the clustered sequences

qiime feature-table summarize \
  --i-table table-hsp-forward-97.qza \
  --o-visualization table-hsp-forward-97.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-hsp-forward-97.qza \
  --o-visualization rep-seqs-hsp-forward-97.qzv

#Export representative sequences in FASTA format for BLAST analysis

qiime tools export \
  --input-path rep-seqs-hsp-forward-97.qza \
  --output-path rep-seqs-hsp-forward-97-export
mv rep-seqs-hsp-forward-97-export/dna-sequences.fasta rep-seqs-hsp-forward-97.fasta

#Blast the sequences with query Eimeria/ reference sequence saved in fasta format

blastn \
  -query rep-seqs-hsp-forward-99.fasta \
  -db nt \
  -remote \
  -entrez_query "Eimeria [ORGN]" \
  -outfmt "6 qseqid sseqid stitle pident qcovs evalue" \
  -max_target_seqs 10 \
  -out blast-results-hsp-forward-99.tsv

# Exports and converts the clustered feature table to a TSV file for downstream analysis
qiime tools export \
  --input-path table-hsp-forward-97.qza \
  --output-path exported-table-hsp-97

biom convert \
  -i exported-table-hsp-97/feature-table.biom \
  -o exported-table-hsp-97/feature-table-hsp-99.tsv \
  --to-tsv


#Convert the original feature table to TSV format to map the clusters to original sequence IDs

qiime tools export \
  --input-path table-hsp-forward.qza \
  --output-path exported-table-original

biom convert \
  -i exported-table-original/feature-table.biom \
  -o exported-table-original/feature-table.tsv \
  --to-tsv
  
setwd("C:/Users/shahn/Downloads/rcc")

# Load required packages
library(readr)
library(Biostrings)
library(dplyr)

# Import the TSV file without column titles
blasted_rep_seqs <- read_tsv("blast-rcc-220-200-98.tsv", col_names = FALSE)

# Assign column names to the data frame
colnames(blasted_rep_seqs) <- c("qseqid", "sseqid", "stitle", "pident", "qcovs", "evalue", "qseq", "sseq")

#get the unique sequences
unique_rep_seqs <- unique(blasted_rep_seqs$qseqid)

#get sequence lengths
blasted_rep_seqs$query_length <- nchar(blasted_rep_seqs$qseq)
blasted_rep_seqs$subject_length <- nchar(blasted_rep_seqs$qseq)
blasted_rep_seqs$differenz <- nchar(blasted_rep_seqs$qseq) - nchar(blasted_rep_seqs$qseq) #checking: it should be 0 for all

# Collapse the data frame
blasted_rep_seqs_collapsed <- blasted_rep_seqs %>%
  group_by(qseqid) %>%
  summarise(
    qseq = paste(qseq, collapse = " "),
    sseq = paste(sseq, collapse = " "),
    query_length = sum(query_length),
    row_count = n(),
    .groups = "drop"  # To ungroup after summarising
  )

write_tsv(blasted_rep_seqs_collapsed, "blasted_rep_seqs_collapsed_98.tsv")


# Filter sequences with length > 300
blasted_rep_seqs_long <- blasted_rep_seqs_collapsed %>% 
  filter(query_length > 300)

# Get the unique sequence IDs that meet our criteria
selected_ids <- unique(blasted_rep_seqs_long$qseqid)
write_tsv(blasted_rep_seqs_long, "blasted_rep_seqs_filtered_98.tsv")

# Read the original FASTA file
original_fasta <- readDNAStringSet("rep-seqs-220-200-98.fasta")

# Extract sequences that match our selected IDs
filtered_sequences <- original_fasta[names(original_fasta) %in% selected_ids]

# Create histogram of sequence lengths
sequence_lengths <- width(filtered_sequences)
hist(sequence_lengths, 
     main = "Distribution of Sequence Lengths (>300 bp)",
     xlab = "Sequence Length (bp)",
     ylab = "Frequency",
     col = "lightblue")

# Save the filtered sequences to a new FASTA file
writeXStringSet(filtered_sequences, 
                filepath = "rep-seqs-rcc-98-filtered.fasta")

# Print summary
cat("Original number of sequences:", length(original_fasta), "\n")
cat("Number of sequences after filtering:", length(filtered_sequences), "\n")
cat("Filtered sequences saved to: rep-seqs-paired-98-filtered.fasta\n")
