import argparse
import os
import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from Levenshtein import distance as levenshtein_distance

def read_txt_file(file_path):
    data = pd.read_csv(file_path, sep='\t', header=None)
    return data

def extract_subsequences(sequence, length):
    return [sequence[i:i+length] for i in range(len(sequence) - length + 1)]

def find_best_match(sequence, barcodes, max_distance):
    best_distance = max_distance
    best_barcode = None
    for barcode in barcodes:
        distance = levenshtein_distance(sequence, barcode)
        if distance < best_distance:
            best_distance = distance
            best_barcode = barcode
    return best_barcode, best_distance

def main(args):
    # Read input files
    fastq_data = read_txt_file(args.txt_file)
    barcodes_data = read_txt_file(args.barcodes)
    
    barcodes = barcodes_data[0].tolist()
    max_distance = args.max_distance
    
    output_data = []

    for index, row in fastq_data.iterrows():
        seq_name = row[0]
        sequence = row[5]
        
        if len(sequence) < args.read_length:
            continue

        sequence = sequence[:args.read_length]
        subsequences = extract_subsequences(sequence, args.barcode_length)

        best_barcode = None
        best_distance = max_distance

        for subseq in subsequences:
            barcode, distance = find_best_match(subseq, barcodes, best_distance)
            if distance < best_distance:
                best_barcode = barcode
                best_distance = distance
        print([seq_name, best_barcode, best_distance])
        output_data.append([seq_name, best_barcode, best_distance])

    output_df = pd.DataFrame(output_data, columns=['Sequence Name', 'Best Matched Barcode', 'Edit Distance'])
    output_df.to_csv('matched_barcodes.tsv', sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--txt_file", type=str, required=True, help="Path to the txt file containing sequences")
    parser.add_argument("--barcodes", type=str, required=True, help="Path to the file containing known barcodes")
    parser.add_argument("--read_length", type=int, default=50, help="Length of the sequence to be considered")
    parser.add_argument("--max_distance", type=int, default=10, help="Maximum allowed edit distance")
    parser.add_argument("--barcode_length", type=int, default=26, help="Length of the barcode")
    
    args = parser.parse_args()
    main(args)
