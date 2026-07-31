#in the terminal before starting (tried to put into main code but so far it isn't working)

minimap2 -ax map-ont NAME_OF_THE_FILE.fasta NAME_OF_THE_FILE.plus.rawreads.fastq.gz | samtools sort -o aligned_reads.bam

samtools index aligned_reads.bam

#later python script, need to fix the first part to run it as one 

import pysam
import csv
import pandas as pd
from Bio import SeqIO
import argparse
import subprocess
import os

#def run_alignment_pipeline(fasta_file, fastq_file, bam_file, threads=4):

    #map_sort_cmd = f"minimap2 -t {threads} -ax map-ont {fasta_file} {fastq_file} | samtools sort -@ {threads} -o {bam_file}"
    #subprocess.run(map_sort_cmd, shell=True, check=True)
    
    #print(f"Indexing BAM file: {bam_file}...")
    #index_cmd = f"samtools index {bam_file}"
    #subprocess.run(index_cmd, shell=True, check=True)
    
    #print("Alignment and indexing complete")

def generate_direct_transposed_table(fasta_path, bam_path, mapping_path, output_csv):
    consensus_record = next(SeqIO.parse(fasta_path, "fasta"))
    consensus_seq = str(consensus_record.seq)
    consensus_length = len(consensus_seq)
    
    print(f"Consensus length: {consensus_length} bp")
    match_pct = []
    with open(mapping_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            match_pct.append(row.get('Match[%]', ''))

    matrix = []
    for i in range(consensus_length):
        mp = match_pct[i] if i < len(match_pct) else ''
        matrix.append([consensus_seq[i], mp])
    print(f"matrix")

    read_count = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        for read in bam_file.fetch():
            if read.is_unmapped:
                continue
          
            read_col = [''] * consensus_length
            
            for q_pos, r_pos in read.get_aligned_pairs():
                if r_pos is not None and q_pos is not None:
                    if r_pos < consensus_length:
                        read_col[r_pos] = read.query_sequence[q_pos]
            
            for i in range(consensus_length):
                matrix[i].append(read_col[i])

            read_count += 1
            if read_count % 1000 == 0:
                print(f"Processed {read_count} reads...")

    print(f"almost there")

    with open(output_csv, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(matrix)
   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Directly generate transposed read alignment matrix.")
    parser.add_argument("-c", "--consensus", default="pDB005-1.fasta", help="Path to input FASTA file")
    parser.add_argument("-b", "--bam", default="aligned_reads.bam", help="Path to input BAM file")
    parser.add_argument("-m", "--mapping", default="pDB005-1.mapping.tsv", help="Path to input mapping TSV file")
    parser.add_argument("-o", "--output", default="read_alignment_table_transposed.csv", help="Path to output CSV file")
    parser.add_argument("-f", "--fastq", default="pDB005-1.plus.rawreads.fastq.gz", help="Optional: Path to FASTQ to run minimap2 automatically first")
    
    args = parser.parse_args()

   #if args.fastq:
        #run_alignment_pipeline(args.consensus, args.fastq, args.bam)

    generate_direct_transposed_table(
        fasta_path=args.consensus,
        bam_path=args.bam,
        mapping_path=args.mapping,
        output_csv=args.output
    )

print(f"done, your result is in {args.output}")
