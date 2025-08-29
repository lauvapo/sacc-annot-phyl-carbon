#!/usr/bin/env python3
################################################################################
#  BLASTP Pipeline                                                             #
#  Sacc2Omics - FungiALab                                                      #
#                                                                              #
#  This script runs BLASTP to identify target genes in the proteomes of        #
#  multiple strains. It automates:                                             #
#    1. Building a multi-FASTA file from strain proteomes.                     #
#    2. Creating a BLASTP database.                                            #
#    3. Running BLASTP with query sequences.                                   #
#    4. Filtering results by identity and coverage.                            #
#    5. Extracting sequences corresponding to filtered hits.                   #
#                                                                              #
#  Author: Laura Varón Pozuelo                                                 #
#  Date:   2025-08-27                                                          #
################################################################################

import os
import subprocess
from Bio import SeqIO

base_dir = "."
proteomes_dir = os.path.join(base_dir, "proteomes")
blast_dir = os.path.join(base_dir, "blast")
db_dir = os.path.join(blast_dir, "db")
genes_dir = os.path.join(base_dir, "general")
hits_dir = os.path.join(base_dir, "hits")
faa_dir = os.path.join(proteomes_dir, "aa_hits")

compilation_fasta = os.path.join(proteomes_dir, "compilation.fasta")
fasta_basename = os.path.splitext(os.path.basename(compilation_fasta))[0]
blast_db_path = os.path.join(db_dir, fasta_basename)
query_files = [os.path.join(genes_dir, f) 
               for f in os.listdir(genes_dir) 
               if f.endswith((".fa", ".fasta"))]
blast_output_file = os.path.join(blast_dir, "blast_results.tsv")
filtered_output_dir = os.path.join(hits_dir, "by_strain")


def create_multifasta(input_dir, output_fasta, strain_list=None):
    """
    Combine proteome FASTA files into a single multi-FASTA file for BLASTP.
    
    Args:
        input_dir (str): Directory containing proteome FASTA files.
        output_fasta (str): Path to the output multi-FASTA file.
        strain_list (list, optional): List of strain names to include 
                                      (default: None).
    
    Returns:
        None
    """
    n = 0
    with open(output_fasta, 'w') as out_f:
        for file in os.listdir(input_dir):
            if file.endswith(('.fasta', '.fa')):
                strain_name = os.path.splitext(file)[0]
                if strain_list and strain_name not in strain_list:
                    continue
                with open(os.path.join(input_dir, file), 'r') as f:
                    print(f"Processing {file}")
                    n += 1
                    for line in f:
                        if line.startswith('>'):
                            out_f.write(f'>{strain_name}_{line[1:]}')
                        else:
                            out_f.write(line)
    print(f"\nMulti-FASTA created with {n} proteome files.\n")


def create_blast_db(fasta_path, db_dir):
    """
    Create a BLASTP database from a FASTA file.
    
    Args:
        fasta_path (str): Path to the input FASTA file.
        db_dir (str): Directory where the BLAST database will be stored.
    
    Returns:
        None
    """
    print("\nCreating BLASTP database...\n")
    fasta_basename = os.path.splitext(os.path.basename(fasta_path))[0]
    output_path = os.path.join(db_dir, fasta_basename)
    cmd = f"makeblastdb -in {fasta_path} -dbtype prot -out {output_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print("BLASTP database created successfully.")
    else:
        print("Error creating BLASTP database:")
        print(result.stderr)


def run_blastp(query, db, output):
    """
    Run BLASTP search with a query against a protein database.
    
    Args:
        query (str): Path to the query FASTA file.
        db (str): Path to the BLAST database (without extension).
        output (str): Path to the BLAST results file (TSV).
    
    Returns:
        None
    """
    print("Running BLASTP...")

    custom_header = [
        "gene", "accession", "%match", "align_length", "mismatch",
        "identity", "gapopen", "query_length", "qstart", "qend", 
        "start", "end", "e-value", "bitscore", "coverage"
    ]
    blast_fields = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "nident", "gapopen", "qlen", "qstart", "qend", "sstart", 
        "send", "evalue", "bitscore", "qcovs"
    ]

    with open(output, 'w') as out:
        out.write("\t".join(custom_header) + "\n")

    cmd = [
        'blastp',
        '-query', query,
        '-db', db,
        '-evalue', '1e-5',
        '-outfmt', '6 ' + ' '.join(blast_fields),
        '-num_threads', '4',
        '-word_size', '3',
        '-gapopen', '11',
        '-gapextend', '1',
        '-matrix', 'BLOSUM62',
        '-threshold', '11',
        '-comp_based_stats', '2',
        '-seg', 'no'
    ]

    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode == 0:
        with open(output, 'a') as out:
            out.write(result.stdout)
        print(f"BLASTP completed. Results saved in {output}")
    else:
        print("Error running BLASTP:")
        print(result.stderr)


def filter_blast_results(blast_file, output_dir, identity_cutoff=90.0):
    """
    Filter BLASTP results by identity and split into per-strain files.
    
    Args:
        blast_file (str): Path to the BLAST results file (TSV).
        output_dir (str): Directory where filtered results will be saved.
        identity_cutoff (float, optional): Identity threshold (default: 90.0).
    
    Returns:
        None
    """
    print(f"\nFiltering results with identity ≥ {identity_cutoff}%...\n")
    os.makedirs(output_dir, exist_ok=True)
    strain_files = {}
    header = [
        "gene", "accession", "%match", "align_length", "mismatch",
        "identity", "gapopen", "query_length", "qstart", "qend", "start", 
        "end", "e-value", "bitscore", "coverage"
    ]

    with open(blast_file, 'r') as f:
        for line in f:
            if line.startswith('gene'):
                continue
            fields = line.strip().split('\t')
            try:
                identity = float(fields[2])
            except ValueError:
                continue

            if identity >= identity_cutoff:
                strain = fields[1].split('_')[0]
                if strain not in strain_files:
                    file_path = os.path.join(output_dir, f'{strain}.tsv')
                    strain_files[strain] = open(file_path, 'w')
                    strain_files[strain].write('\t'.join(header) + '\n')
                strain_files[strain].write('\t'.join(fields) + '\n')

    for f in strain_files.values():
        f.close()
    print(f"Filtered results saved in {output_dir}")


def extract_hit_sequences(
        blast_dir: str, multifasta: str, output_dir: str) -> None:
    """
    Extract sequences of BLASTP hits into strain-specific FASTA files.
    
    Args:
        blast_dir (str): Directory containing filtered per-strain BLAST result 
                         TSVs.
        multifasta (str): Path to the combined multi-FASTA file.
        output_dir (str): Directory to save strain-specific FASTA files.
    
    Returns:
        None
    """
    print("\nExtracting protein sequences from filtered hits...\n")
    os.makedirs(output_dir, exist_ok=True)

    seq_dict = SeqIO.to_dict(SeqIO.parse(multifasta, "fasta"))

    for file in os.listdir(blast_dir):
        if not file.endswith(".tsv"):
            continue
        strain_name = file.replace('.tsv', '')
        protein_ids = set()

        with open(os.path.join(blast_dir, file), 'r') as f:
            next(f)
            for line in f:
                sseqid = line.strip().split('\t')[1]
                protein_ids.add(sseqid)

        output_faa = os.path.join(output_dir, f"{strain_name}.faa")
        with open(output_faa, 'w') as out_fasta:
            found = 0
            for pid in protein_ids:
                if pid in seq_dict:
                    SeqIO.write(seq_dict[pid], out_fasta, "fasta")
                    found += 1
            print(f"{strain_name}: {found} sequences extracted.")

    print(f"\nFASTA files written to {output_dir}")


if __name__ == '__main__':
    print("\nStarting BLASTP pipeline...\n")
    create_multifasta(proteomes_dir, compilation_fasta)
    create_blast_db(compilation_fasta, db_dir)
    for query_file in query_files:
        query_name = os.path.splitext(os.path.basename(query_file))[0]
    
        blast_output_file = os.path.join(
            blast_dir, f"{query_name}_blast_results.tsv")
        filtered_output_dir = os.path.join(hits_dir, "by_strain", query_name)
        faa_out_dir = os.path.join(faa_dir, query_name)
    
        run_blastp(query_file, blast_db_path, blast_output_file)
        filter_blast_results(blast_output_file, filtered_output_dir)
        extract_hit_sequences(
            filtered_output_dir, compilation_fasta, faa_out_dir)
        
    print("\nPipeline completed.\n")
