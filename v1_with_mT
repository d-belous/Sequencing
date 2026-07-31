#no file flexibility, but worked 

import pysam
import csv
import pandas as pd
from Bio import SeqIO

def generate_alignment_table(fasta_path, bam_path, output_csv):
  
    consensus_record = next(SeqIO.parse(fasta_path, "fasta"))
    consensus_seq = str(consensus_record.seq)
    consensus_length = len(consensus_seq)
    
    print(f"Consensus length: {consensus_length} bp")
    print(f"Parsing alignments and writing to {output_csv}...")


    with pysam.AlignmentFile(bam_path, "rb") as bam_file, open(output_csv, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(list(consensus_seq))
        
        read_count = 0

        for read in bam_file.fetch():
            if read.is_unmapped:
                continue

            row = [''] * consensus_length

            aligned_pairs = read.get_aligned_pairs()
            
            for q_pos, r_pos in aligned_pairs:

                if r_pos is not None and q_pos is not None:
                    if r_pos < consensus_length: 
                        base = read.query_sequence[q_pos]
                        row[r_pos] = base # Place the base in the exact consensus column
            
            writer.writerow(row)
            read_count += 1

    print(f"Done! Successfully assigned {read_count} reads to the consensus.")

if __name__ == "__main__":
    CONSENSUS_FILE = "pDB005-1.fasta"
    BAM_FILE = "aligned_reads.bam"
    OUTPUT_CSV = "read_alignment_table.csv"
    
    generate_alignment_table(CONSENSUS_FILE, BAM_FILE, OUTPUT_CSV)

  
  pd.read_csv("read_alignment_table.csv", header=None).T.to_csv("read_alignment_table_transposed.csv", index=False, header=False)
