#!/usr/bin/env python3
################################################################################
#  Contigs sequences extraction                                                #
#  Sacc2Omics - FungiALab                                                      #
#                                                                              #
#  This script extracts gene sequences from strain proteomes based on BLASTP   #
#  hit coordinates. It automates:                                              #
#    1. Parsing TSV hit files per strain.                                      #
#    2. Mapping coordinates to contigs in proteome FASTA files.                #
#    3. Extracting sequences (with start/stop positions).                      #
#    4. Validating sequences (start codon = M).                                #
#    5. Writing multi-FASTA files per gene with extracted sequences.           #
#                                                                              #
#  Author: Laura Varón Pozuelo                                                 #
#  Date:   2025-08-27                                                          #
################################################################################


import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


base_dir = "."

hits_dir = os.path.join(base_dir, "hits")
proteomes_dir = os.path.join(base_dir, "proteomes")
multifasta_dir = os.path.join(base_dir, "multifasta")
os.makedirs(multifasta_dir, exist_ok=True)


def parse_tsv_hits(hits_dir: str) -> dict:
    """
    Parse TSV files containing BLAST hits and build a dictionary of strains.

    attr:
        hits_dir (str): Path to the directory containing TSV files.

    returns:
        dict: Dictionary with strain as key and list of [
              gene, contig, start, end] as values.
    """
    strains = {}
    for file in os.listdir(hits_dir):
        if file.endswith('.tsv') and not file.startswith('compilation'):
            with open(os.path.join(hits_dir, file), 'r') as f:
                strain = file.replace('.tsv', '')
                for line in f:
                    if not line.startswith('gene'):
                        line = line.strip().split('\t')
                        if strain not in strains:
                            strains[strain] = []
                        strains[strain].append(
                            [line[0], line[1], int(line[8]), int(line[9])])
    return strains


def extract_sequences(strains: dict, proteomes_dir: str) -> dict:
    """
    Extract sequences from proteome FASTA files based on TSV coordinates.

    attr:
        strains (dict): Dictionary with strain as key and 
                        [gene, contig, start, end] as values.
        proteomes_dir (str): Path to directory containing proteome FASTA files.

    returns:
        dict: Dictionary with gene as key and list of (contig_name, sequence) 
              tuples as values.
    """
    coord_chr = {}

    for strain, contigs in strains.items():
        for contig_info in contigs:
            gene, contig_id, start, end = contig_info

            # Fix inverted coordinates
            if start > end:
                start, end = end, start

            contig_name = 'chr' + contig_id
            fasta_file = os.path.join(proteomes_dir, f"{strain}.fasta")

            if os.path.exists(fasta_file):
                for seq_record in SeqIO.parse(fasta_file, "fasta"):
                    if seq_record.id == contig_name:
                        seq_record.id = strain
                        seq = seq_record.seq[start - 1:end]

                        if gene not in coord_chr:
                            coord_chr[gene] = []

                        # Only keep if sequence starts with 'M'
                        if seq.startswith('M'):
                            if (f'{strain}_{contig_name}', seq) not in coord_chr[gene]:
                                coord_chr[gene].append((f'{strain}_{contig_name}', seq))
                        else:
                            print(
                                f"Sequence for {gene} in {strain}_{contig_name} doesn't start with 'M'")
                        break
    return coord_chr


def write_multifasta(gene: str, sequences: list, multifasta_dir: str) -> None:
    """
    Write extracted sequences to a multifasta file for each gene.

    attr:
        gene (str): Gene name used for the output filename.
        sequences (list): List of (contig_name, sequence) tuples.
        multifasta_dir (str): Directory where multifasta files will be saved.

    returns:
        None: Writes one FASTA file per gene.
    """
    filename = f"{gene}_multi.fasta"
    filepath = os.path.join(multifasta_dir, filename)
    unique_sequences = {contig: seq for contig, seq in sequences}
    seq_records = [
        SeqRecord(
            seq, id=contig, description="") 
        for contig, seq in unique_sequences.items()]
    with open(filepath, "w") as f:
        if len(seq_records) > 1:
            SeqIO.write(seq_records, f, "fasta")


if __name__ == "__main__":
    strains = parse_tsv_hits(hits_dir)
    coord_chr = extract_sequences(strains, proteomes_dir)

    for gene, sequences in coord_chr.items():
        write_multifasta(gene, sequences, multifasta_dir)

    print("Multi-FASTA files created in ", multifasta_dir)
