# the_validator

Repository for a snakemake pipeline for validating snakemake bioinformatics pipelines.

## background

Input reference sequence, generate...

<img src="https://github.com/aineniamh/artic-noro/blob/master/primer-schemes/noro2kb/V2/noro2kb.amplicons.png">

This pipeline starts with the basecalled reads, gathers them together and applies a read length filter, discarding reads less than 1000 bases or greater than 3000 bases. It then demultiplexes the reads using ``qcat`` and bins them into ``amplicons`` by BLASTing them against a reference database of amplicon nucleotide sequences. The best genotype is assigned per amplicon based on ``blastn`` results. Individual bed files and reference files are created using custom scripts and an instance of Nick Loman's ``minion`` pipeline is run for each amplicon for each sample. The output is then tidied up into respective directories.

## pipeline steps

<img src="https://github.com/aineniamh/artic-noro/blob/master/dag_one_sample.svg">

## setup

Although not a requirement, an install of conda will make the setup of this pipeline on your local machine much more streamlined. I have created an ``artic-noro`` conda environment which will allow you to access all the software required for the pipeline to run. To install conda, visit here https://conda.io/docs/user-guide/install/ in a browser. 

> *Recommendation:* Install the `64-bit Python 3.6` version of Miniconda

Once you have a version of conda installed on your machine, clone this repository by typing into the command line:

```bash
git clone https://github.com/aineniamh/artic-noro.git
```

Build the conda environment by typing:

```bash
conda env create -f artic-noro/envs/artic-noro.yaml
```

To activate the environment, type:

```bash
source activate artic-noro
```

To deactivate the environment, enter:

```bash
conda deactivate
```

## customising the pipeline

To run the analysis using snakemake, you will want to customise the ``config.yaml`` file.

Inside the ``config.yaml`` file, you can change your MinION ``run_name`` and the barcode names. Ensure ```path_to_fast5``` and ```path_to_fastq``` point to where your data is.

## running the pipeline

To start the pipeline, in a terminal window in the artic-noro directory, simply enter:

```bash
snakemake
```

If you wish to run your pipeline using more than one core (**recommended**), enter:

```bash
snakemake --cores X
```

where X is the number of threads you wish to run.

## pipeline description

1. setup ``artic fieldbioinformatics`` package \
Automatic setup of this on startup of the pipeline. Gives the user access to ``artic minion`` script for step below.
2. gather \
Parses all of the basecalled fastq files from ``guppy``, applies a length filter that can be customised in the ``config.yaml`` file and writes the reads to a single file ``run_name_all.fastq``. This script also searches the fastq directories for ``sequencing_summary`` files and combines them into a single file: ``run_name_sequencing_summary.txt``. These files will be output in the ``pipeline_output`` directory.
3. demultiplex_qcat \
For each read in the ``run_name_all.fastq`` file, identifies barcodes and outputs reads into respective files, binned by barcode. These files appear in the ``demultiplexed`` directory, in ``pipeline_output``.
4. make_amplicon_files \
Using the primer-scheme bed file, amplicon summary files are generated using a custom python script. The script produces a figure showing the amplicon span over the reference genomes, a csv file with amplicon information and a fasta file with extracted amplicon sequences.
5. fastq_to_fasta \
Fastq files are converted to fasta files using ``seqtk`` for the purposes of blasting.
6. blastn \
For each ``barcode.fastq`` file, each read is blasted against a sequence database containing all amplicon sequences.
7. bin \
This step parses each blast output and assesses for each read what the best blast hit is. The reads are then binned by amplicon and, for each amplicon, the best reference sequence is calculated. This determines which reference is most suited to take forward into nanopolish for each amplicon for each barcode. It also creates the respective bed file for the ``artic minion`` pipeline to use.

<!-- 4. minimap2_index \
Indexes a panel of reference sequences for minimap2.
5. minimap2 \
For each barcode, maps the reads against the panel of reference sequences and produces a ``.paf`` file.
6. find_top_reference \
For each ``barcode``, identifies the reference with the greatest number of reads mapping to it and creates a new reference file ``primer-schemes/noro2kb/V_barcode/barcode.reference.fasta`` and a new bed file ``primer-schemes/noro2kb/V_barcode/barcode.scheme.bed``.
7. minimap_to_top_reference \
Re-maps the reads for each demultiplexed file against their respective top reference and outputs a sam file in ``pipeline_output/best_ref_mapped_reads/``.
8. Quick reference generation \
``samtools`` is used to sort the reads and ``bcftools`` is then used to call variants, normalise for indels and call a quick consensus sequence ``pipeline_output/consensus/{barcode}.cns.fasta``. This consensus sequence is then renamed and saved in ``primer-schemes/noro2kb/V_barcode/`` as ``barcode.reference.fasta``. -->
8. nanopolish_index \
Creates the nanopolish index necessary for running nanopolish in the next step. It accesses the gathered fastq and sequencing summary files from step 2 and also the signal-level fast5 data.
9. artic_minion \
The ``artic minion`` pipeline, written by Nick Loman, is then run for each barcode in order to generate a high-quality consensus sequence, using an approach informed by signal-level data. This pipeline performs the following steps:
    * Maps against a given reference and sorts reads using ``bwa`` and ``samtools`` respectively.
    * Runs the ``artic align_trim`` script. This script takes in a bed file and your alignment and assesses whether the primers are correctly paired according to the bed file, discarding reads that are not, and normalises the read coverage across the genome. It is run twice, first to trim off the barcodes and the primers and second to just trim off the barcodes.
    * Loads the ``nanopolish index`` created in step 8.
    * Runs ``nanopolish variants`` twice, on the barcode-and-primer-trimmed bam and on the barcode-trimmed bam.
    * Generates a variant frequency plot.
    * Runs ``margin_cons``, a custom script that filters the variants, masking sites that do not reach the depth threshold of 20 and do not reach a quality threshold of 200, and produces a consensus sequence with 'N' masking on the relevant sites. It uses the vcf from nanopolish without primer-trimming but the primer-trimmed bam file so that primer sequences do not count towards depth calculation. A report is also generated.
10. organise_minion_output \
Moves artic_minion output files into respective ``pipeline_output/minion_output/barcode`` directories on completion of the pipeline.
11. generate_genome (Not done yet) \
Overlays the consensus sequences for each amplicon and creates a whole-genome reference sequence. 