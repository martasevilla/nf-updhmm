#!/usr/bin/env python3
import sys, csv, re, argparse

HEADER = [
    "fam_id","proband_id","mother_id","father_id",
    "path_vcf_proband","path_vcf_mother","path_vcf_father",
    "path_sv_proband","path_sv_mother","path_sv_father"
]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate nf-core/updhmm samplesheet",
        epilog="Usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"
    )
    parser.add_argument("FILE_IN", help="Input samplesheet CSV")
    parser.add_argument("FILE_OUT", help="Validated samplesheet CSV")
    return parser.parse_args()

def print_error(msg, line=""):
    sys.stderr.write(f"ERROR: {msg}\n{line}\n")
    sys.exit(1)

def check_file(path, allowed_exts, allow_dash=False):
    if allow_dash and path == "-":
        return
    if not any(path.endswith(ext) for ext in allowed_exts):
        print_error(f"File {path} does not end with one of {allowed_exts}")

def check_samplesheet(file_in, file_out):
    seen_fams = set()
    seen_ids = set()

    with open(file_in, newline="") as fin, open(file_out, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        if reader.fieldnames != HEADER:
            print_error(f"Header mismatch! Found {reader.fieldnames}, expected {HEADER}")

        writer = csv.DictWriter(fout, fieldnames=HEADER)
        writer.writeheader()

        for i, row in enumerate(reader, start=2):
            fam = row["fam_id"].strip()
            if not fam:
                print_error("fam_id is empty!", f"Line {i}")
            if fam in seen_fams:
                print_error("Duplicate fam_id!", f"Line {i}")
            seen_fams.add(fam)

            for role in ["proband_id","mother_id","father_id"]:
                sid = row[role].strip()
                if not re.match(r"^[A-Za-z0-9_.-]+$", sid):
                    print_error(f"Invalid characters in {role}: {sid}", f"Line {i}")
                if sid in seen_ids:
                    print_error(f"Duplicate sample ID: {sid}", f"Line {i}")
                seen_ids.add(sid)

            check_file(row["path_vcf_proband"], [".vcf.gz"])
            check_file(row["path_vcf_mother"], [".vcf.gz"])
            check_file(row["path_vcf_father"], [".vcf.gz"])
            check_file(row["path_sv_proband"], [".bed"], allow_dash=True)
            check_file(row["path_sv_mother"], [".bed"], allow_dash=True)
            check_file(row["path_sv_father"], [".bed"], allow_dash=True)

            writer.writerow(row)

def main():
    args = parse_args()
    check_samplesheet(args.FILE_IN, args.FILE_OUT)

if __name__ == "__main__":
    main()
