# snp-id
Tools for identifying SNP markers for iPLEX and other genotyping assays

## Prerequisites
- BioPython
- Pysam
- Faidx

## Introduction
The **snp-id** scripts can be used to identify SNPs suitable for iPLEX or other genotyping assays.  Low coverage (<10x) WGS data with multiple biological replicates for each species are required.  High quality diagnostic SNPs in low-polymorphic regions are chosen and output in a format suitable for downstream assay design, e.g. as input for the MassARRAY Typer Assay Designer.

## Workflow Example
Assume that you are developing an iPLEX assay to distinguish three closely related species – _species_A_, _species_B_ and _species_C_, and have sequenced 5 biological replaces for each species.  The following inputs must be prepared or generated:

### Reference Genomes
A reference genome for each species in FASTA format is required.  If a reference genome is not readily available for your species, you can simply choose to perform a de novo assembly on one of the replicates and use that as the reference.  Relatively low coverage draft genome sequences should be sufficient.

### Multi-genome Alignment
A genome alignment in XFMA format must be provided.  You can use [progressiveMauve](http://darlinglab.org/mauve/user-guide/progressivemauve.html) to produce the alignment.

### Replicate Alignments
Each replicate needs to be aligned to its respective reference, e.g. with [BWA](https://github.com/lh3/bwa).  The BAM file for each replicate is required as part of the input.

### Input JSON File
Once the above data have been generated and assuming that they are all in the same directory, create a JSON file as below:

```
{
  "reference": {
    "spA": "species_a.fasta",
    "spB": "species_b.fasta ",
    "spC": "species_c.fasta"
        },
  "bams": {
    "spA": [
      "species_a_repl1.bam",
      "species_a_repl2.bam",
      "species_a_repl3.bam",
      "species_a_repl4.bam",
      "species_a_repl5.bam"
    ],
    "spB": [
      "species_b_repl1.bam",
      "species_b_repl2.bam",
      "species_b_repl3.bam",
      "species_b_repl4.bam",
      "species_b_repl5.bam"
    ],
    "spC": [
      "species_c_repl1.bam",
      "species_c_repl2.bam",
      "species_c_repl3.bam",
      "species_c_repl4.bam",
      "species_c_repl5.bam"
    ]
  },
  "mauve_alignment": "spA-spB-spC.xmfa"
}
```

### Running the program
The only argument for the main program, **search_iplex.py**, is the JSON input file.  Diagnostic SNPs and flanking sequences are written to standard output while logs are written to standard error, so you should capture them in separate files, e.g.
```
$ <path_to_snp-id>/search_iplex.py input.json >out.txt 2>log
```
Depending on the size of the genomes, this may take a few hours.

### Identifying contaminants (optional)
You may want to identify possible contaminants in the output and exclude them from downstream analysis.  The **blast_iplex.py** script can be used to blast the output of **search_iplex.py** against the NCBI nt database, e.g.

```
$ BLASTDB=<path_to_blastdb>/nt <path_to_snp-id>/blast_iplex.py out.txt > out.blast.txt
```

## Availability
The **snp-id** scripts are released under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) and are freely available at [Github](https://github.com/ClockLabX/snp-id/).
