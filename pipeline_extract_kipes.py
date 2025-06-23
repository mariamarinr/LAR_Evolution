#!/usr/bin/env python3
"""
Extract high-quality LAR/DFR/ANR protein sequences from KIPES output
and map them to corresponding CDS sequences using a CSV table.

Steps:
1. Parse KIPES output folders and extract AA sequences with 100% conserved residues.
2. Parse AA FASTA headers (format: species@sequence_id).
3. Read CDS mapping table (CSV with columns 'Species', 'CDS_file').
4. Match protein sequences to CDS and write final CDS FASTA.

Usage:
    python3 pipeline_extract_kipes_cds.py \
        --main-folder PATH_TO_KIPES_OUTPUT \
        --cds-table PATH_TO_CDS_TABLE_CSV \
        --output-cds OUTPUT_CDS_FASTA \
        [--tmp-aa TEMP_AA_FASTA] \
        [--skipped-headers SKIPPED_HEADERS_FILE] \
        [--dry-run] \
        [--clean-temp]

Arguments:
    --main-folder          Path to KIPES output folder (required)
    --cds-table            CSV file with 'Species', 'CDS_file' columns (required)
    --output-cds           Output CDS FASTA file (required)
    --tmp-aa               Temporary AA FASTA file (default: temp-aa-sequences.fasta)
    --skipped-headers      Output file for skipped AA headers (default: skipped_headers.txt)
    --dry-run              Dry run mode: parse and match sequences, but do not write output
    --clean-temp           Delete temporary AA FASTA after run
"""

import os
import csv
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_aa_sequences(main_folder, tmp_aa_output, dry_run=False):
    print("[STEP 1] Extracting AA sequences from KIPES output...")

    total_written = 0
    with open(tmp_aa_output, "w") as outfile:
        for subfolder in os.listdir(main_folder):
            subfolder_path = os.path.join(main_folder, subfolder)
            if not os.path.isdir(subfolder_path):
                continue

            summary_file = os.path.join(subfolder_path, "summary.txt")
            if not os.path.exists(summary_file):
                continue

            ids_to_extract = set()
            with open(summary_file, "r") as summary:
                for line in summary:
                    columns = line.strip().split()
                    if len(columns) >= 4:
                        seq_id = columns[0]
                        gene_type = columns[1]
                        conserved = columns[3]
                        if gene_type.startswith(("LAR", "DFR", "ANR")) and conserved == "100.0":
                            ids_to_extract.add(seq_id)

            if ids_to_extract:
                for gene_type in ["LAR", "DFR", "ANR"]:
                    fasta_file = os.path.join(subfolder_path, "final_pep_files", f"{gene_type}.fasta")
                    if os.path.exists(fasta_file):
                        for record in SeqIO.parse(fasta_file, "fasta"):
                            if record.id in ids_to_extract:
                                if not dry_run:
                                    SeqIO.write(record, outfile, "fasta")
                                total_written += 1

    print(f"[INFO] Extracted {total_written} AA sequences to {tmp_aa_output}")
    return total_written

def parse_aa_fasta(tmp_aa_output):
    print("[STEP 2] Parsing extracted AA FASTA...")

    aa_records = list(SeqIO.parse(tmp_aa_output, "fasta"))
    species_to_ids = {}
    id_to_header = {}
    skipped_headers = []

    for record in aa_records:
        header = record.id
        if "@" not in header:
            print(f"[SKIP] Malformed header (no @): {header}")
            skipped_headers.append(header)
            continue

        species_key, rest = header.split("@", 1)
        gene_id = rest
        species_to_ids.setdefault(species_key, set()).add(gene_id)
        id_to_header[(species_key, gene_id)] = header

    return aa_records, species_to_ids, id_to_header, skipped_headers

def read_cds_table(cds_table_file):
    print("[STEP 3] Reading CDS table...")

    species_to_cdsfiles = {}
    with open(cds_table_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            species_name = row.get("Species", "").strip()
            cds_path = row.get("CDS_file", "").strip()
            if species_name and cds_path and Path(cds_path).is_file():
                species_to_cdsfiles.setdefault(species_name, []).append(cds_path)

    return species_to_cdsfiles

def extract_matching_cds(species_to_ids, id_to_header, species_to_cdsfiles, final_cds_output, dry_run=False):
    print("[STEP 4] Extracting matching CDS sequences...")

    output_records = []
    total_found = 0
    total_missing = 0

    for species_key, gene_ids in species_to_ids.items():
        matched_species = None
        for species_name in species_to_cdsfiles:
            if species_key.lower() in species_name.lower():
                matched_species = species_name
                break

        if not matched_species:
            print(f"[WARNING] No matching CDS entry for species '{species_key}'")
            total_missing += len(gene_ids)
            continue

        cds_files = species_to_cdsfiles[matched_species]
        for cds_file in cds_files:
            try:
                cds_dict = SeqIO.to_dict(SeqIO.parse(cds_file, "fasta"))
            except Exception as e:
                print(f"[ERROR] Failed to parse CDS file {cds_file}: {e}")
                continue

            for gene_id in gene_ids:
                if gene_id in cds_dict:
                    seq = cds_dict[gene_id]
                else:
                    matches = [seq for key, seq in cds_dict.items() if gene_id in key]
                    seq = matches[0] if matches else None

                if seq:
                    new_header = id_to_header[(species_key, gene_id)]
                    new_record = SeqRecord(seq.seq, id=new_header, description="")
                    output_records.append(new_record)
                    total_found += 1
                else:
                    print(f"[NOT FOUND] Gene ID '{gene_id}' not found in {Path(cds_file).name}")
                    total_missing += 1

    if output_records and not dry_run:
        SeqIO.write(output_records, final_cds_output, "fasta")
        print(f"[DONE] Extracted {len(output_records)} CDS sequences to {final_cds_output}")

    return total_found, total_missing

def main():
    parser = argparse.ArgumentParser(description="Pipeline: Extract KIPES LAR/DFR/ANR AA → Map to CDS → Write CDS FASTA")
    parser.add_argument("--main-folder", required=True, help="Path to KIPES output folder")
    parser.add_argument("--cds-table", required=True, help="CSV with columns 'Species', 'CDS_file'")
    parser.add_argument("--output-cds", required=True, help="Output FASTA for CDS sequences")
    parser.add_argument("--tmp-aa", default="temp-aa-sequences.fasta", help="Temporary AA FASTA path")
    parser.add_argument("--skipped-headers", default="skipped_headers.txt", help="Output file for skipped AA headers")
    parser.add_argument("--dry-run", action="store_true", help="Dry run (no files written)")
    parser.add_argument("--clean-temp", action="store_true", help="Delete temp AA file after run")

    args = parser.parse_args()

    aa_written = extract_aa_sequences(args.main_folder, args.tmp_aa, dry_run=args.dry_run)
    aa_records, species_to_ids, id_to_header, skipped_headers = parse_aa_fasta(args.tmp_aa)
    species_to_cdsfiles = read_cds_table(args.cds_table)
    total_found, total_missing = extract_matching_cds(
        species_to_ids, id_to_header, species_to_cdsfiles, args.output_cds, dry_run=args.dry_run
    )

    print("\n=== SUMMARY ===")
    print(f"AA sequences extracted: {aa_written}")
    print(f"AA records parsed: {len(aa_records)}")
    print(f"CDS sequences found: {total_found}")
    print(f"CDS sequences missing: {total_missing}")
    print(f"Skipped malformed AA headers: {len(skipped_headers)}")
    print("================")

    if skipped_headers and not args.dry_run:
        with open(args.skipped_headers, "w") as f:
            for header in skipped_headers:
                f.write(f"{header}\n")
        print(f"[INFO] Skipped headers saved to {args.skipped_headers}")

    if args.clean_temp and Path(args.tmp_aa).exists():
        os.remove(args.tmp_aa)
        print(f"[INFO] Temporary file {args.tmp_aa} deleted.")

if __name__ == "__main__":
    main()
