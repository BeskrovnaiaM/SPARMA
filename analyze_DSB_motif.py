import pysam
from collections import defaultdict
import bisect
import subprocess

BAM_FILE = "example.bam"
GENOME_FASTA = "example.fasta"
WINDOW_SIZE = 10
TARGET_DISTANCE = 0
OUTPUT_HIST = "output_histogram.tsv"
OUTPUT_BED = "output_pairs.bed"
OUTPUT_FASTA = "output_sequences.fa"
OUTPUT_FILTERED_FASTA = "output_filtr_sequences.fa"
WEBLOGO_OUTPUT = "output_logo.png"


def extract_5prime_coords(bam_file):
    plus_ends = defaultdict(list)
    minus_ends = defaultdict(list)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            chrom = read.reference_name
            name = read.query_name
            if not read.is_reverse:
                plus_ends[chrom].append((read.reference_start, name))
            else:
                minus_ends[chrom].append((read.reference_end - 1, name))
    return plus_ends, minus_ends


def find_pairs_with_distance(plus_ends, minus_ends, window_size, target_dist):
    selected_reads = set()
    distance_reads = defaultdict(list)
    for chrom in set(plus_ends) & set(minus_ends):
        plus_list = sorted(plus_ends[chrom])
        minus_list = sorted(minus_ends[chrom])
        minus_positions = [m[0] for m in minus_list]
        for p_pos, p_name in plus_list:
            left = bisect.bisect_left(minus_positions, p_pos - window_size - 1)
            right = bisect.bisect_right(minus_positions, p_pos + window_size)
            for m_pos, m_name in minus_list[left:right]:
                dist = (m_pos - p_pos) + 1
                if -window_size <= dist <= window_size:
                    distance_reads[dist].append((p_name, m_name))
                    if dist == target_dist:
                        selected_reads.add(p_name)
                        selected_reads.add(m_name)
    return selected_reads, distance_reads


def create_bed_for_plus_reads(bam_file, selected_reads, output_bed):
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_bed, 'w') as bed:
        for read in bam:
            if read.is_reverse or read.is_secondary or read.is_supplementary:
                continue
            if read.query_name not in selected_reads:
                continue
            end_pos = read.reference_start
            start = max(0, end_pos - 5)
            end = end_pos + 5
            bed.write(f"{read.reference_name}\t{start}\t{end}\t{read.query_name}\t.\t+\n")


def run_bedtools(genome_fasta, input_bed, output_fasta):
    cmd = f"bedtools getfasta -fi {genome_fasta} -bed {input_bed} -s -name -fo {output_fasta}"
    subprocess.run(cmd, shell=True, check=True)


def filter_fasta_by_length(input_fasta, output_fasta, length=10):
    with open(input_fasta) as fin, open(output_fasta, 'w') as fout:
        header, seq = "", ""
        for line in fin:
            if line.startswith('>'):
                if header and len(seq) == length:
                    fout.write(header + seq + '\n')
                header, seq = line, ""
            else:
                seq += line.strip()
        if header and len(seq) == length:
            fout.write(header + seq + '\n')


def run_weblogo(input_fasta, output_logo):
    cmd = (f"weblogo -f {input_fasta} -o {output_logo} -A dna -F png -c classic -U probability --annotate '-5,-4,-3,-2,-1,+1,+2,+3,+4,+5' --size large --resolution 1200")
    subprocess.run(cmd, shell=True, check=True)


def plot_histogram(distance_reads, window_size):
    histogram = {dist: len(pairs) for dist, pairs in distance_reads.items()}
    full_hist = {d: histogram.get(d, 0) for d in range(-window_size, window_size + 1)}
    with open(OUTPUT_HIST, 'w') as f:
        for d in sorted(full_hist):
            f.write(f"{d}\t{full_hist[d]}\n")


def main():
    plus_ends, minus_ends = extract_5prime_coords(BAM_FILE)
    selected_reads, distance_reads = find_pairs_with_distance(plus_ends, minus_ends, WINDOW_SIZE, TARGET_DISTANCE)
    print(f"reads at target distance {TARGET_DISTANCE}: {len(selected_reads)}")
    create_bed_for_plus_reads(BAM_FILE, selected_reads, OUTPUT_BED)
    run_bedtools(GENOME_FASTA, OUTPUT_BED, OUTPUT_FASTA)
    filter_fasta_by_length(OUTPUT_FASTA, OUTPUT_FILTERED_FASTA)
    run_weblogo(OUTPUT_FILTERED_FASTA, WEBLOGO_OUTPUT)
    plot_histogram(distance_reads, WINDOW_SIZE)
    print(f"histogram: {OUTPUT_HIST}")
    print(f"bed: {OUTPUT_BED}")
    print(f"sequences: {OUTPUT_FILTERED_FASTA}")
    print(f"weblogo: {WEBLOGO_OUTPUT}")


if __name__ == "__main__":
    main()
