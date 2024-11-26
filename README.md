# BKC
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/bkc/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/bkc/releases)

BKC is a program for <i>k</i>-mer counting in barcoded FASTQ/FASTA files obtained from technologies like 10x, Visium.
BKC is one of many projects developed by [REFRESH Bioinformatics Group](https://refresh-bio.github.io).
If your data are not barcoded and you need to count <i>k</i>-mers, please visit [KMC](https://github.com/refresh-bio/KMC) project.

## Quick start
### Getting the executable
The simplest way to get the BKC executable is to download the newest release for an appropriate operating system from [BKC releases]().

### Downloading sample data
The example requires [SRAtools](https://github.com/ncbi/sra-tools) to be on path for downloading one run from [SRA](https://www.ncbi.nlm.nih.gov/sra).
```
prefetch SRR6211483
fastq-dump --split-3 SRR6211483
```
You should obtain two FASTQ (suffixed by `_1` and `_2`) of sizes approx. 1GB and 2.2GB.

# Prapare input file with a description of data
Then, it is necessary to prepare a single CSV file describing the input.
Each run should be given as a single line and FASTQs should be comma-separated (`_1` file with CBC and UMI) must be given as the first.
In the example, there is a single run only, so it suffices to:
```
echo "SRR6211483_1.fastq,SRR6211483_2.fastq" > fl.txt
```

### Counting and dumping <i>k</i>-mers
In the barcoded data, the first (`_1`) file usually stores CBC and UMI. 
CBCs are used to separate data from different cells.
You should provide CBC and UMI lengths according to the technology used for sequencing.
For 10x v3, CBCs are of length 16, and UMIs are of length 12. 
For the older v2 chemistry, CBCs are still of length 16, but UMIs are shorter (10).
In our example, the data are from 10x v3 technology.
The length of <i>k</i>-mers is specified with `--leader_len` option.

```
./bkc --input_name fl.txt --cbc_len 16 --umi_len 12 --leader_len 27 --output_name SRR6211483_single.bkc
./bkc_dump --input_name SRR6211483_single.bkc --output_name SRR6211483_single.txt
```

### Counting and dumping <i>k</i>-mer pairs
The origins of BKC are in our [SPLASH](https://github.com/refresh-bio/SPLASH) project, where it was used as a highly-specialized counting tool for pairs of <i>k</i>-mers.
Then, BKC grew and was offered as a more general tool.
That's why BKC offers a mode for counting <i>k</i>-mer pairs (with optional gap).
To use it in this way, you have to specify the `--mode pair` explicitly and provide the length of the 2nd <i>k</i>-mer (follower).
```
./bkc --mode pair --input_name fl.txt --cbc_len 16 --umi_len 12 --leader_len 27 --follower_len 27 --output_name SRR6211483_pair.bkc
./bkc_dump --input_name SRR6211483_pair.bkc --output_name SRR6211483_pair.txt
```

## Compiling from source
BKC can also be cloned and compiled from sources. 
For Linux you need GCC/G++ 11 or newer.
```
git clone --recurse-submodules https://github.com/refresh-bio/bkc.git
cd bkc
make -j
```

For MacOS you need GCC/G++ 11 or newer.
Please be aware of setting CC and CXX in the right way. 
```
git clone --recurse-submodules https://github.com/refresh-bio/bkc.git
cd bkc
make -j CC=gcc-12 CXX=g++12   # usualy gcc is an alias for clang, so you need to specify gcc/g++ version
```

In both cases, the `bkc` and `bkc_dump` executable files will be placed in the `bin\` subdirectory.

## bkc details
### Main options
* `--mode <single|pair>` &ndash; selects the mode:
  * `single` &ndash; counting of $k$-mers (default),
  * `pair` &ndash; counting of $k$-mer pairs,
  * `filter` &ndash; CBC filtering and UMI deduplication.
* `--cbc_len <int>` &ndash; CBC len (default: 16, min: 10, max: 16).
* `--umi_len <int>` &ndash; UMI len (default: 12, min: 8, max: 16).
* `--leader_len <int>` &ndash; length of $k$-mer in `single` mode or length of the 1st $k$-mer of a pair in `pair` mode (default: 27, min: 1, max: 31).
* `--follower_len <int>` &ndash; length of the 2nd $k$-mer of a pair in `pair` mode (default: 0, min: 0, max: 31).
* `--gap_len <int>` &ndash; in `pair` mode, the leader and follower $k$-mers can be separated by some gap (default: 0, min: 0, max: 4294967295).
* `--n_threads <int>` &ndash; no. threads (default: 8, min: 0, max: 256).
* `--canonical` &ndash; turns on canonical k-mers (default: false); works only in single mode.
* `--verbose <int>` &ndash; verbosity level (default: 0, min: 0, max: 2).

### Input data specification
* `--input_format <fasta|fastq>` &ndash; select input format (default: fastq).
* `--input_name <file_name>` &ndash; file name with a list of pairs (comma separated) of barcoded files; 1st contains CBC+UMI.
* `--technology <10x|visium>` &ndash; sequencing technology (default: 10x).
* `--soft_cbc_umi_len_limit <int>` &ndash; tolerance of CBC+UMI len (default: 0, min: 0, max: 1000000000). It happens that `_1` reads are longer than CBC_len+UMI_len. With this option, you can specify how much longer they can be. BKC will, however, use only a prefix of such reads.
* `--cbc_filtering_thr <int>` &ndash; [UMItools](https://github.com/CGATOxford/UMI-tools) applies CBC filtering (by removing rare CBCs). BKC follows the same strategy if you specify the threshold as 0 (default). Nevertheless, you can also specify the number of reads the CBC must contain to prevent it from filtering out. (default: 0, min: 0, max: 4294967295)
* `--allow_strange_cbc_umi_reads` &ndash; use this option to prevent the application from crashing when the CBC+UMI read length is outside the acceptable range (either shorter than CBC_len+UMI_len or longer than CBC_len+UMI_len+soft_cbc_umi_len_limit). Use with care as such strange reads highly suggest that there is something wrong with the data.
* `--apply_cbc_correction` &ndash; apply CBC correction (similar to UMI tools).

### Output options
* `--output_format <bkc|splash>` &ndash; allows to specify the output format (default: bkc). As said above, BKC originated in the SPLASH project, in which we use a slightly different output format.
* `--output_name <file_name>` &ndash;- output file name (default: ./results.bkc).
* `--sample_id <int>` &ndash; sample id (default: 0). Usually you do not need to care about this id. Nevertheless, if you want to merge outputs from several runs of BKC for various data and want to be able to distinguish between them, specifying `sample_id` can be helpful.
* `--n_splits <int>` &ndash; if you want to have the output to be split into several files, you can use this switch (default: 1, min: 1, max: 256). The $k$-mers are split according to some hash value. This can be helpful if you want to process BKC output in parallel by many threads/processes.
* `--log_name <file_name>` &ndash; path to cbc log files (default: ); if not provided, log will not be produced. This can be helpful if you want to see how the threshold between the trusted and non-trusted CBCs was selected.
* `--filtered_input_path <string>` &ndash; path to filtered input files (default: ). BKC can work as a filter, providing CBC filtering and UMI deduplication of reads. If you want to use it in this way, you need to provide also `--export_filtered_input_mode` option.
* `--export_filtered_input_mode <none|first|second|both>` &ndash;- specifies which reads will be outputted (default: none).
* `--max_count <int>` &ndash; sometimes it is unimportant how big the counters are if we know they are big enough. Here, you can specify the max. counter value (default: 65535, min: 1, max: 4294967295). This allows some space to be saved in the output file.
* `--zstd_level <int>` &ndash; BKC output files are internally compressed. We use some custom encoding followed by [ZSTD](https://github.com/facebook/zstd) compression. With this option, you can redefine the default compression level(default: 6, min: 0, max: 19).

### Filtering options
* `--predefined_cbc <file_name>` &ndash; sometimes, it is helpful to provide the list of trusted CBCs rather than looking for them in the reads (default: ). The format of this file depends on the `--technology` parameter:
  * For `10x`, it should be a plain list of CBCs (one per line).
  * For `visium`, it should be in the format defining the tissue CBCs, i.e., each line should much the regex `([ACGT]+)-(.+),([0-9]+),[0-9]+,[0-9]+,[0-9]+,[0-9]+`, where the 1st block is for CBC, and the 4th should be `1`.
* `--poly_ACGT_len <int>` &ndash; all leaders containing polyACGT of this length will be filtered out (0 means no filtering) (default: 0, min: 0, max: 31).
* `--artifacts <file_name>` &ndash; a path to artifacts, each leader containing artifact will be filtered out. This can be useful if you want to remove some patterns from the reads.
* `--apply_filter_illumina_adapters` &ndash; if used, leaders containing Illumina adapters will be filtered out (the adapters are hardcoded in BKC).
* `--leader_sample_counts_threshold <int>` &ndash; keep only leaders with counts > leader_sample_counts_threshold (default: 5, min: 0, max: 255). Use the small values with care as you can obtain huge output files.

## bkc_dump details
### Options
* `--input_name <file_name>` &ndash; name of BKC file.
* `--output_name <file_name>` &ndash; output file name (default: ./results.txt).
* `--n_splits <int>` &ndash; no. splits (default: 1, min: 1, max: 256). Must be the same as `--n_splits` used in `bkc` run.

### Dump format
The format is TSV (tab-separated) and each line is composed of:
* sample id &ndash; specified with `--sample_id` option of bkc,
* CBC,
* leading $k$-mer,
* following $k$-mer &ndash; this is present only in `pair` mode of bkc,
* counter.

## More examples
Now let's use 10x sample data.
First, we need to download and unpack the files
```
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar
```

This dataset contains two runs for the same sample.
To produce BKC input file we need:
```
echo "pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz" > fl2.txt
echo "pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz,pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz" >> fl2.txt
```

Then, we can count $k$-mers:
```
./bkc --input_name fl2.txt --cbc_len 16 --umi_len 12 --leader_len 27 --n_threads 16 --output_name pbmc_1k_v3_single.bkc
./bkc_dump --input_name pbmc_1k_v3_single.bkc --output_name pbmc_1k_v3_single.txt
```
or pairs of $k$-mers:
```
./bkc --mode pair --input_name fl2.txt --cbc_len 16 --umi_len 12 --leader_len 27 --follower_len 27 --n_threads 16 --output_name pbmc_1k_v3_pair.bkc
./bkc_dump --input_name pbmc_1k_v3_pair.bkc --output_name pbmc_1k_v3_pair.txt
```

To just do CBC filtering and UMI deduplication:
```
./bkc --mode filter --input_name fl2.txt --cbc_len 16 --umi_len 12 --export_filtered_input_mode both --n_threads 16
```
