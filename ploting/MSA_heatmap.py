import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os

def parse_msa_file(filepath):
    """
    Parses an MSA file. Handles standard FASTA format or 
    block-interleaved text formats (like the one in image_11f4a8.png).
    """
    sequences = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()

    if not lines:
        return sequences

    # Check if the file is standard FASTA
    if lines[0].startswith(">"):
        name = None
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                if name not in sequences:
                    sequences[name] = ""
            elif name and line:
                sequences[name] += line
    else:
        # Parse block/interleaved text format
        for line in lines:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # Expecting at least two columns: [Name] ... [Sequence Chunk]
            if len(parts) >= 2:
                name = parts[0]
                seq = parts[-1]
                if name not in sequences:
                    sequences[name] = ""
                sequences[name] += seq
                
    return sequences

def plot_msa_heatmap(sequences):
    """
    Generates a discrete 10-color heatmap from the parsed sequences 
    and saves it as a 4:3 ratio PNG file.
    """
    names = list(sequences.keys())
    if not names:
        print("No sequences found. Please check your file.")
        return

    # Determine maximum sequence length
    seq_len = max(len(seq) for seq in sequences.values())

    # Initialize matrices: one for numeric scores, one for characters
    scores = np.full((len(names), seq_len), np.nan)
    chars = np.full((len(names), seq_len), "", dtype=object)

    for i, name in enumerate(names):
        seq = sequences[name]
        for j, char in enumerate(seq):
            chars[i, j] = char
            # Convert digits to integers for the heatmap
            if char.isdigit():
                scores[i, j] = int(char)
            # "-" and letters remain np.nan (will be plotted as white)

    # Define Custom Discrete Colormap (N=10 for values 0 through 9)
    color_stops = ["#00008B", "#ADD8E6", "#FFB6C1", "#800020"]
    cmap = LinearSegmentedColormap.from_list("discrete_gradient", color_stops, N=10)
    cmap.set_bad(color='white') # Handles np.nan (gaps and letters)

    # Set figure size to a 4:3 width-to-height ratio (e.g., 16x12 inches)
    width = 16
    height = 12
    fig, ax = plt.subplots(figsize=(width, height))
    
    # Plot the matrix using discrete color-stepping (vmax/vmin shifted for integers 0-9)
    cax = ax.imshow(scores, cmap=cmap, aspect='auto', vmin=-0.5, vmax=9.5)

    # Overlay text for "-" and letters
    rows, cols = scores.shape
    for i in range(rows):
        for j in range(cols):
            val = chars[i, j]
            if val == '-':
                ax.text(j, i, '-', ha='center', va='center', color='black', fontsize=8)
            elif val.isalpha():
                ax.text(j, i, val, ha='center', va='center', color='black', fontsize=7)

    # Axis formatting
    ax.set_yticks(np.arange(len(names)))
    ax.set_yticklabels(names, family='monospace', fontsize=9)
    ax.set_xlabel("Alignment Position", fontsize=12)
    ax.set_title("Protein Alignment Score Heatmap", fontsize=14)

    # Add a discrete colorbar with ticks for 0-9
    cbar = plt.colorbar(cax, ax=ax, pad=0.01, ticks=np.arange(0, 10))
    cbar.set_label("Score", fontsize=12)

    plt.tight_layout()
    
    # Save directly to a 4:3 PNG file
    output_filename = "msa_heatmap_sortase.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved successfully as {output_filename}")
    plt.close()

if __name__ == "__main__":
    # Replace with the path to your txt or fasta file
    file_path = "sortase_al.txt" 
    
    if os.path.exists(file_path):
        parsed_seqs = parse_msa_file(file_path)
        plot_msa_heatmap(parsed_seqs)
    else:
        print(f"File not found: {file_path}. Please update the file_path variable.")
