import gzip
import sys

def parse_length_range(range_str):
    """Parses a string like '[3000; 5000]' into a tuple of ints."""
    cleaned = range_str.strip().strip("[]").replace(";", ",")
    parts = [int(p.strip()) for p in cleaned.split(",")]
    if len(parts) != 2:
        raise ValueError(f"Invalid range format: {range_str}. Expected format like [3000; 5000]")
    return sorted(parts)

def process_fastq(input_file, range1_str, range2_str, out1_file, out2_file, basecall_model=None):
    min_l1, max_l1 = parse_length_range(range1_str)
    min_l2, max_l2 = parse_length_range(range2_str)

    print(f"Filtering reads and preserving/injecting headers:")
    print(f"  Range 1: [{min_l1}, {max_l1}] -> {out1_file}")
    print(f"  Range 2: [{min_l2}, {max_l2}] -> {out2_file}")
    if basecall_model:
        print(f"  Injecting basecall model tag: {basecall_model}")

    count1 = 0
    count2 = 0
    skipped = 0

    open_in = gzip.open if input_file.endswith(".gz") else open
    open_out = gzip.open if out1_file.endswith(".gz") else open

    with open_in(input_file, "rt") as fin, \
         open_out(out1_file, "wt") as fout1, \
         open_out(out2_file, "wt") as fout2:
        
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            if not qual:
                break  # Malformed file safety check

            read_len = len(seq.strip())

            # Format header, optionally injecting the basecall model tag if missing
            cleaned_header = header.strip()
            if basecall_model and "basecall_model_version_id" not in cleaned_header and "DS=" not in cleaned_header:
                cleaned_header = f"{cleaned_header} basecall_model_version_id={basecall_model}"
            
            formatted_header = cleaned_header + "\n"

            if min_l1 <= read_len <= max_l1:
                fout1.write(f"{formatted_header}{seq}{plus}{qual}")
                count1 += 1
            elif min_l2 <= read_len <= max_l2:
                fout2.write(f"{formatted_header}{seq}{plus}{qual}")
                count2 += 1
            else:
                skipped += 1

    print("\nProcessing complete:")
    print(f"  Saved {count1} reads to {out1_file}")
    print(f"  Saved {count2} reads to {out2_file}")
    print(f"  Skipped {skipped} reads outside both ranges.")

if __name__ == "__main__":
    input_fastq = "Evol-kan-day1.minus.rawreads copy.fastq.gz"
    
    # Define your two target read length boundaries
    range_one = "[3000; 4000]"
    range_two = "[6500; 8000]"
    
    output_one = "rawreads_range1.fastq.gz"
    output_two = "rawreads_range2.fastq.gz"

    # Set your sequencing basecaller model here (e.g., from your MinKNOW run report)
    # Leave as None if your original headers already contain model tags and you just want strict preservation.
    target_basecall_model = "dna_r10.4.1_e8.2_400bps_hac@v4.2.0"

    try:
        process_fastq(input_fastq, range_one, range_two, output_one, output_two, basecall_model=target_basecall_model)
    except FileNotFoundError:
        print(f"Error: Could not find input file '{input_fastq}'. Please ensure it is in the working directory.", file=sys.stderr)
        sys.exit(1)
