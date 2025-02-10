#!/bin/bash
set -xeuo pipefail

if [[ $OSTYPE = linux-gnu ]]; then
    color="--color=always"
else
    color=""
fi

function strobealign() {
    echo "Testing '${@}'" >&2
    target/debug/rstrobes "${@}"
#     if ! target/debug/rstrobes "${@}" 2> testlog.txt; then
#         cat testlog.txt
#         echo "Failure"
#         exit 1
#     fi
}

function diff() {
    if ! env diff -u ${color} "$1" "$2"; then
      echo "Failure running 'diff $1 $2'"
      exit 1
    fi
}

cargo build

# Unit tests
cargo test

# Ensure the binary is available
samtools --version > /dev/null

# Single-end SAM, M CIGAR operators
strobealign --no-PG tests/phix.fasta tests/phix.1.fastq > phix.se.m.sam
if samtools view phix.se.m.sam | cut -f6 | grep -q '[X=]'; then false; fi
rm phix.se.m.sam

# Single-end PAF
strobealign -x tests/phix.fasta tests/phix.1.fastq | tail -n 11 > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Single-end PAF (stdin input)
cat tests/phix.1.fastq | strobealign -x tests/phix.fasta - | tail -n 11 > phix.se.paf
diff tests/phix.se.paf phix.se.paf
rm phix.se.paf

# Paired-end PAF
strobealign -x tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq | tail -n 11 > phix.pe.paf
diff tests/phix.pe.paf phix.pe.paf
rm phix.pe.paf

# Single-end abundance estimation
strobealign --aemb tests/phix.fasta tests/phix.1.fastq > phix.abun.se.txt
diff tests/phix.abun.se.txt phix.abun.se.txt
rm phix.abun.se.txt

# Paired-end abundance estimation
strobealign --aemb tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq > phix.abun.pe.txt
diff tests/phix.abun.pe.txt phix.abun.pe.txt
rm phix.abun.pe.txt

# Build a separate index
strobealign --no-PG -r 150 tests/phix.fasta tests/phix.1.fastq > without-sti.sam
strobealign -r 150 -i tests/phix.fasta
strobealign --no-PG -r 150 --use-index tests/phix.fasta tests/phix.1.fastq > with-sti.sam
diff without-sti.sam with-sti.sam
rm without-sti.sam with-sti.sam

# Create index requires -r or reads file
if strobealign --create-index tests/phix.fasta > /dev/null 2> /dev/null; then false; fi

# --details output is proper SAM
strobealign --details tests/phix.fasta tests/phix.1.fastq tests/phix.2.fastq 2> /dev/null | samtools view -o /dev/null
strobealign --details tests/phix.fasta tests/phix.1.fastq 2> /dev/null | samtools view -o /dev/null

strobealign -C --no-PG --rg-id=1 --rg=SM:sample --rg=LB:library tests/phix.fasta tests/phix.tags.fastq > with-tags.sam
diff tests/phix.tags.sam with-tags.sam
rm with-tags.sam

# Secondary alignments

# No secondary alignments on phix
strobealign --no-PG tests/phix.fasta tests/phix.1.fastq > no-secondary.sam
strobealign --no-PG -N 5 tests/phix.fasta tests/phix.1.fastq > with-secondary.sam
test $(samtools view -f 0x100 -c with-secondary.sam) -eq 0
rm no-secondary.sam with-secondary.sam

# Secondary alignments for repeated phiX
cp tests/phix.fasta repeated-phix.fasta
echo ">repeated_NC_001422" >> repeated-phix.fasta
sed 1d tests/phix.fasta >> repeated-phix.fasta
strobealign --no-PG repeated-phix.fasta tests/phix.1.fastq > no-secondary.sam
strobealign --no-PG -N 5 repeated-phix.fasta tests/phix.1.fastq > with-secondary.sam
test $(samtools view -f 0x100 -c with-secondary.sam) -gt 0

# Removing secondary alignments gives same result as not producing them in the first place
samtools view -h --no-PG -F 0x100 with-secondary.sam > with-secondary-only-primary.sam
diff no-secondary.sam with-secondary-only-primary.sam
rm no-secondary.sam with-secondary.sam with-secondary-only-primary.sam repeated-phix.fasta

echo "Success"
