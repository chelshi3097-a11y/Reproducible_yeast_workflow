import gzip

files = [
    r"C:\Users\Maria\Downloads\Bioinformatics_projects\rnaseq_yeast\data_raw_reads\ERR458500_yeast_MUT1.fastq.gz",
    r"C:\Users\Maria\Downloads\Bioinformatics_projects\rnaseq_yeast\data_raw_reads\ERR458501_yeast_MUT2.fastq.gz"
]

for f in files:
    with gzip.open(f, "rt") as fh:
        next(fh)  # skip header
        seq = next(fh).strip()  # sequence line
        print(f"{f}: Read length = {len(seq)}")
