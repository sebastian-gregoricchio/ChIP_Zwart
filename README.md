[![Snakemake](https://img.shields.io/badge/snakemake-≥7.24.0-brightgreen.svg)](https://snakemake.github.io)
![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/ChIP_Zwart)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/sebastian-gregoricchio/ChIP_Zwart/LICENSE.md/LICENSE.md)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/ChIP_Zwart?style=social)](https://github.com/sebastian-gregoricchio/ChIP_Zwart/fork)
<!-- ![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/ChIP_Zwart)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/ChIP_Zwart) --->

# ChIP-seq for Zwartlab
## Introduction
In this repository is provide a snakemake-based pipeline for the analyses of ChIP-seq data. It allows the mapping of fastq files (both paired- and single-end) as well as other downstream analyses sturting from the bam files (i.e., mapping filtering, sample correlation, peak calling).

## Installation an dependencies
To install the pipeline it is required to download this repository and the installation of a conda environment is strongly recommended.

### Conda initialization
If never done before proceed to the following steps:
* initialize conda in your environment by: `/opt/miniconda3/bin/conda init`
* actualize the terminal with: `source ~/.bashrc`

On your terminal, it should appear `(base) your.name@harris:~$`. If this is not the case, run `conda activate base`.

Ensure you have the proper conda path (/opt/miniconda3/bin/conda) by running: `which conda`


### Environment installation
To avoid packages version incompatibility a yam file with fixed packages versions is provided in this repository.

For the installation, follow the following steps:
* Place yourself in the directory where the repository should be downloaded with `cd </target/folder>`
* download the GitHub repository with `git clone https://github.com/sebastian-gregoricchio/ChIP_Zwart`, or click on *Code > Download ZIP* on the [GitHub page](https://github.com/sebastian-gregoricchio/ChIP_Zwart)
* install the conda environment from the yaml environment file contained in the repository:<br>
`conda env create -f </target/folder>/ChIP_Zwart/workflow/envs/chip_zwart_condaEnv_stable.yaml`
* activate the environment: `conda activate chip_zwart` (if the env is not activated the pipeline won't work!)

<br/><br/>



## How to run the pipeline
The snakemake pipeline requires at least two files: a) the `.snakefile`, containing all the rules that will be run; b) the `configuration.yaml` file, in which the user can define and customize all the parameters for the different pipeline steps. <br>
Hereafter, the running commands for DNA-mapping and ChIP-seq peak calling will be described for both single and paired-end data.


### DNA-mapping
This short pipeline performs the mapping into a reference genome upon trimming of the raw fastq reads by [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Further, a filter on the mapping quality (MAPQ) is applied and duplicated reads are marked. Notice that in the case of paired-end reads, when present, UMIs (Unique Molecule Identifiers) sequence is added to the indexes ones in the read name. This is allows the marking of the duplicated reads in a UMI-aware manner (reads/fragments that have exactly the same sequence but different UMI-sequence are not marked as duplicates).

Few information must be provided to the pipeline:
* the source fastq directory
* the output directory where you want your results to be stored (if not already available, the pipeline will make it for you)
* whether your data are paired- or single-end
* the genome to use

All the other parameters are already available in the `configfile_DNAmapping.yaml` file or hard-coded in the snakemake file.


To partially avoid unexpected errors during the execution of the pipeline, a so called 'dry-run' is strongly recommended. Indeed, adding a `-n` flag at the end of the snakemake running command will allow snakemake to check that all links and file/parameters dependencies are satisfied before to run the "real" processes. This command will therefore help the debugging process. <br>
*Always activate your environment, otherwise the pipeline won't be able to find the packages required for the analyses.*


**Paired-end** (the `\` must be used every time you go to a new line)
```shell
snakemake \
--cores 20 \
-s </target/folder>/ChIP_Zwart/workflow//DNAmapping.snakefile \
--configfile </target/folder>/ChIP_Zwart/config/configfile_DNAmapping.yaml \
--config \
fastq_directory="/path/to/pairedEnd/fastq_data" \
output_directory="/path/to/results/directory/" \
paired_end="True" \
genome="hg38" \
-n
```

**Single-end** (the `\` must be used every time you go to a new line)
```shell
snakemake \
--cores 20 \
-s </target/folder>/ChIP_Zwart/workflow//DNAmapping.snakefile \
--configfile </target/folder>/ChIP_Zwart/config/configfile_DNAmapping.yaml \
--config \
fastq_directory="/path/to/singleEnd/fastq_data" \
output_directory="/path/to/results/directory/" \
paired_end="False" \
genome="hg19" \
-n
```

If no errors occur, the pipeline can be run with the same command but without the final `-n` flag:

Notice that the absence of errors does not mean that the pipeline will run without any issues; the "dry-run" is only checking whether all the resources are available. <br>

<br/><br/>


#### DNA-mapping workflow
Here after you can see the full potential workflow of the single-end and paired-end DNA-mapping pipeline:

![DNA mapping workflow](https://raw.githubusercontent.com/sebastian-gregoricchio/ChIP_Zwart/main/resources/DNAmapping_workflow_SE.PE.png)


<br/><br/>

#### DNA-mapping results
The results structure is the following:
* *01_trimmed_fastq* -> fastq.gz files that underwent trimming by cutadapt
* *02_BAM* -> mapped reads (bam) filtered for MAPQ, mate-fixed, duplicates marked and eventually UMI-fixed
* *03_quality_controls* -> here you can find the fastQC on the trimmed fastq (if required), with the corresponding multiQC report, as well the multiQC report (flagstat + MarkDuplicates) for the filtered bams


Here an example directory tree (paired-end run):
<pre>
<b><em>output_folder</em></b>
├── <b>01_trimmed_fastq</b>
│   ├── <em>sample</em>_R1_trimmed.fastq.gz
│   └── <em>sample</em>_R2_trimmed.fastq.gz
│
├── <b>02_BAM</b>
│   ├── <em>sample</em>_mapq20_mdup_sorted.bam
│   ├── <em>sample</em>_mapq20_mdup_sorted.bai
│   ├── <b>BWA_summary</b>
│   │   └── <em>sample</em>.BWA_summary.txt
│   ├── <b>flagstat</b>
│   │   └── <em>sample</em>_mapq20_mdup_sorted_flagstat.txt
|   ├── <b>MarkDuplicates_metrics</b>
│   │   └── <em>sample</em>_MarkDuplicates_metrics.txt
│   └── <b>umi_metrics</b>  ### (if UMI present) ##
│       └── <em>sample</em>_UMI_metrics.txt
|
└── <b>03_quality_controls</b>
    ├── <b>multiQC_bam_filtered</b>
    │   └── multiQC_bam_filtered.html
    ├── <b>trimmed_fastq_fastqc</b>
    │   ├── <em>sample</em>_R1_trimmed_fastqc.html
    │   ├── <em>sample</em>_R1_trimmed_fastqc.zip
    │   ├── <em>sample</em>_R2_trimmed_fastqc.html
    │   └── <em>sample</em>_R2_trimmed_fastqc.zip
    └── <b>05_Quality_controls_and_statistics</b>
        └── multiQC_report_trimmed_fastq.html
</pre>

<br/><br/>

#### DNA-mapping config file

| **Parameter**   |  **Description**   |
|------------:|:----------------|
| *umi_present* | Default: `True`. True/False to indicate whether the data contain UMIs (ignored for single-end data). |
| *fastq_suffix* | Default: `".fastq.gz"`. String with the suffix of the source fastq files. |
| *read_suffix* | Default: `['_R1', '_R2']`. A python-formatted list with two strings containing the suffix used to indicate read1 and read2 respectively. In the case of single end reads, only the first value will be read. If your single data do not have any read-prefix set this parameter to: `['', '']` (blank). |
| *cutadapt_trimm_options* | Default: `''` (blank). String indicating additional user-specific values to pass to cutadapt. |
| *fw_adapter_sequence* | Default: `"AGATCGGAAGAGC"`. Sequence of the adapter1 (flag `-a` of cutadapt). |
| *rv_adapter_sequence* | Default: `"AGATCGGAAGAGC"`. Sequence of the adapter2 (flag `-A` of cutadapt). |
| *run_fastq_qc* | Default: `False`. True/False to indicate whether to run the fastQC on the trimmed fastq file. |
| *use_bwamem2* | Default: `False`. True/False to define whether to run [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) instead of [bwa](https://bio-bwa.sourceforge.net/). |
| *bwa_options* | Default: `''` (blank). String indicating additional user-specific values to pass to bwa. |
| *remove_duplicates* | Default: `False`. True/False to define whether remove the duplicates from the bam files (if true the tag in the bams will be *_dedup* instead of *_mdup*). |
| *MAPQ_threshold* | Default: `20`. All reads with a mapping quality (MAPQ) score lower than this value will be filtered out from the bam files. |


<br/><br/>

<hr style="border:2px solid blue">


### CRAM-TO-BAM
It may happen that you start with `.cram` files provided by the GCF. However, the peak-calling pipeline works only with `.bam` files. Therefore, to convert the crams to bams you can use the mini-pipeline provided in this repository. You can also ask to the pipeline to rename your files with the "wz number" indicating `rename_zwart="True"` (otherwise use `rename_zwart="False"`):

```shell
snakemake \
--cores 10 \
-s </target/folder>/ChIP_Zwart/workflow/cramToBam.snakefile \
--config \
cram_directory="/path/to/input/cram_folder" \
bam_out_directory="/path/to/output/bam_folder" \
genome="hg38" \
rename_zwart="True"
```
(the `\` must be used every time you go to a new line)

**NOTE**: remember that the genome used for the conversion must match with the one used to generate the crams. Hence, if the crams do not come from the NKI-GCF please ask for assistance.

<br/><br/>

<hr style="border:2px solid blue">



<br/><br/>

### ChIP-seq peak calling
To facilitate the analyses of the ChIP-seq analyses in the Zwart lab, it is strongly recommended to rename your files so that the files contain the wz number. To do that refer to the section [*Renaming files*](https://github.com/csijcs/snakepipes#renaming-files) of the [*Joe's GitHub*](https://github.com/csijcs/snakepipes#renaming-files).

The pipeline requires a sample configuration file which provides information about ChIP-Input pairs and the type of peak calling to perform (broad or narrow). <br>
This configuration file must be in a tab-delimited txt file format (with column names) containing the following information (respect the column order):

| **target_id**   |   **input_id**   |   **broad**   |
|:----------------|:-----------------|:--------------|
| sample_A        |   input_A-B      |    false      |
| sample_B        |   input_A-B      |    false      |
| sample_C        |   input_C        |    true       |


Few additional information must be provided to the pipeline:
* the source bam directory (e.g. *rename* folder)
* the output directory where you want your results to be stored (if not already available, the pipeline will make it for you)
* whether your data contain UMIs
* whether your data are paired- or single-end
* the genome to use
* the path to the sample configuration table

All the other parameters are already available in the `configfile_peakcalling.yaml` file or hard-coded in the snakemake file.


To partially avoid unexpected errors during the execution of the pipeline, a so called 'dry-run' is strongly recommended. Indeed, adding a `-n` flag at the end of the snakemake running command will allow snakemake to check that all links and file/parameters dependencies are satisfied before to run the "real" processes. This command will therefore help the debugging process. <br>
*Always activate your environment, otherwise the pipeline won't be able to find the packages required for the analyses.*


**Paired-end** (the `\` must be used every time you go to a new line)
```shell
snakemake \
--cores 10 \
-s </target/folder>/ChIP_Zwart/workflow/peakcalling.snakefile \
--configfile </target/folder>/ChIP_Zwart/config/configfile_peakCalling.yaml \
--config \
runs_directory="/path/to/rename" \
output_directory="/path/to/results/directory/" \
sample_config_table="/path/to/sample_configuration.txt" \
paired_end="True" \
umi_present="True" \
genome="hg38" \
-n
```

**Single-end** (the `\` must be used every time you go to a new line)
```shell
snakemake \
--cores 10 \
-s </target/folder>/ChIP_Zwart/workflow/peakcalling.snakefile \
--configfile </target/folder>/ChIP_Zwart/config/configfile_peakCalling.yaml \
--config \
runs_directory="/path/to/rename" \
output_directory="/path/to/results/directory/" \
sample_config_table="/path/to/sample_configuration.txt" \
paired_end="False" \
genome="hg19" \
-n
```

If no errors occur, the pipeline can be run with the same command but without the final `-n` flag:

Notice that the absence of errors does not mean that the pipeline will run without any issues; the "dry-run" is only checking whether all the resources are available. <br>


#### Peak calling workflow
Here after you can see the full potential workflow of the paired-end and single-end ChIP-seq pipeline:

**a) PAIRED-END**
![PE workflow](https://raw.githubusercontent.com/sebastian-gregoricchio/ChIP_Zwart/main/resources/peakcalling_workflow_PE.png)

<br/><br/>

**b) SINGLE-END**
![PE workflow](https://raw.githubusercontent.com/sebastian-gregoricchio/ChIP_Zwart/main/resources/peakcalling_workflow_SE.png)

<br/><br/>



#### Peak calling results
The results structure is the following:
* *01_BAM_filtered* -> bams filtered for mapping quality (MAPQ) and with the duplicates marked/removed
* *02_fastQC_on_BAM_filtered* -> individual fastQC for each filtered bam
* *03_bigWig_bamCoverage* -> bigWig of the bam coverage normalized ([RCPG](https://deeptools.readthedocs.io/en/develop/content/help_glossary.html?highlight=RPGC#abbreviations) = Reads Per Genomic Content) or not (raw_coverage) depending on the sequencing depth
* *04_Called_peaks* -> peaks called with macs (de-blacklisted). If single-end also a folder with the output of [`phantompeakqualtools`](https://www.encodeproject.org/software/phantompeakqualtools/). If the bed file does not contain already the 'chr' for the "canonical" chromosomes, it will be added in a separated file
* *05_Quality_controls_and_statistics* -> this folder contains sample correlations heatmaps and PCAs, a multiQC report containing multiple info (number of reads, duplicates, peak counts and fragmenth lenght, phantom results), statistics on the called peaks (FRiP, number, etc.)


Here an example directory tree:
<pre>
<b><em>output_folder</em></b>
├── <b>01_BAM_filtered</b>
│   ├── <em>sample</em>_mapq20_mdup_sorted.bam
│   ├── <em>sample</em>_mapq20_mdup_sorted.bai
│   ├── <b>flagstat</b>
│   │   └── <em>sample</em>_mapq20_mdup_sorted_flagstat.txt
|   ├── <b>MarkDuplicates_metrics</b>
│   │   └── <em>sample</em>_MarkDuplicates_metrics.txt
│   └── <b>umi_metrics</b>  ### (if UMI present) ##
│       └── <em>sample</em>_UMI_metrics.txt
|
├── <b>02_fastQC_on_BAM_filtered</b>
│   ├── <em>sample</em>_sorted_woMT_dedup_fastqc.html
│   └── <em>sample</em>_sorted_woMT_dedup_fastqc.zip
|
├── <b>03_bigWig_bamCoverage</b>
│   ├── <b>raw_coverage</b>
│   │   └── <em>sample</em>_mapq20_mdup_raw.coverage_bs10.bw
│   └── <b>RPGC_normalized</b>
│       └── <em>sample</em>_mapq20_mdup_RPGC.normalized_bs10.bw
│
├── <b>04_Called_peaks</b>   ### BAM if Single-End ###
│   ├──<b>phantom</b>
│   │   └── <em>sample</em>.phantom.ssp.out
│   ├── <em>sample</em>.filtered.BAMPE_peaks_chr.narrowPeak
│   ├── <em>sample</em>.filtered.BAMPE_peaks.narrowPeak
│   └── <em>sample</em>.filtered.BAMPE_peaks.xls
|
└── <b>05_Quality_controls_and_statistics</b>
    ├── <b>multiQC</b>
    │   └── multiQC_report.html
    ├── <b>peaks_stats</b>
    │    └── all_samples_FRiP_report.tsv
    ├── <b>plotFingerprint</b>   ### optional ###
    |   ├── <b>quality_metrics</b>
    |   │   └── em>sample</em>_fingerPrinting_quality_metrics.txt
    |   └── em>sample</em>_fingerPrinting_plot.pdf
    ├── <b>sample_comparisons_atPeaks</b>
    │   ├── all_peaks_merged_sorted.bed
    │   ├── multiBigWigSummary_matrix_atPeaks.npz
    │   ├── sample_pearson.correlation_heatmap_atPeaks.pdf
    │   ├── sample_spearman.correlation_heatmap_atPeaks.pdf
    │   ├── sample_correlation_PCA.1-2_heatmap_atPeaks.pdf
    │   └── sample_correlation_PCA.2-3_heatmap_atPeaks.pdf
    └── <b>sample_comparisons_wholeGenome</b>
        ├── all_peaks_merged_sorted.bed
        ├── multiBigWigSummary_matrix_atPeaks.npz
        ├── sample_pearson.correlation_heatmap_wholeGenome.pdf
        ├── sample_spearman.correlation_heatmap_wholeGenome.pdf
        ├── sample_correlation_PCA.1-2_heatmap_wholeGenome.pdf
        └── sample_correlation_PCA.2-3_heatmap_wholeGenome.pdf
</pre>


<br/><br/>

#### Peak calling config file

| **Parameter**   |  **Description**   |
|------------:|:----------------|
| *bam_suffix* | Default: `".bam"`. String with the suffix of the source bam files. |
| *umi_present* | Default: `True`. True/False to indicate whether the data contain UMIs (ignored for single-end data). |
| *remove_duplicates* | Default: `False`. True/False to define whether remove the duplicates from the bam files (if true the tag in the bams will be *_dedup* instead of *_mdup*). |
| *MAPQ_threshold* | Default: `20`. All reads with a mapping quality (MAPQ) score lower than this value will be filtered out from the bam files. |
| *bigWig_binSize* | Default: `10`. Size, in bp, of the bins used to compute the normalized bigWig files. |
| *use_macs3* | Default: `False`. True/False to define whether to run [macs3](https://github.com/macs3-project) instead of [macs2](https://pypi.org/project/MACS2/). |
| *macs_qValue_cutoff* | Default: `0.01`. False Discovery Ratio (FDR) (q-value) cutoff used by [MACS](https://github.com/macs3-project/MACS) to filter the significant peaks. |
| *perform_plotFingerprint* | Default: `False`. True/False to define whether perform the finger printing (Lorenz curve). |
| *perform_fragmentSizeDistribution* | Default: `False`. True/False to define whether to plot the fragment size distribution (Paired-end only). |
| *fragment_length* | Default: `200`. Size in bp of the virtual fragment length at which each read should be extended in order to perform the plotFingerprint.  |
| *correlation_heatmap_colorMap* | Default: `'PuBuGn'`. A string indicating the color gradient pattern to use for the correlation heatmaps. This value is passed to matplotlib/seaborn. Therefore, available options (see [matplotlib page](https://matplotlib.org/stable/tutorials/colors/colormaps.html) for examples) are the following: 'Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'Dark2', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired', 'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cividis', 'cool', 'coolwarm', 'copper', 'cubehelix', 'flag', 'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'icefire', 'inferno', 'jet', 'magma', 'mako', 'nipy_spectral', 'ocean', 'pink', 'plasma', 'prism', 'rainbow', 'rocket', 'seismic', 'spring', 'summer', 'tab10', 'tab20', 'tab20b', 'tab20c', 'terrain', 'twilight', 'twilight_shifted', 'viridis', 'vlag', 'winter'. |
|*zScore_heatmap_color*| Default: `"seismic"`. A string indicating the color gradient pattern to use for the peak score heatmaps. This value is passed to matplotlib/seaborn. Therefore, available options (see [matplotlib page](https://matplotlib.org/stable/tutorials/colors/colormaps.html) for examples) are the following: 'Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'Dark2', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired', 'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cividis', 'cool', 'coolwarm', 'copper', 'cubehelix', 'flag', 'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'icefire', 'inferno', 'jet', 'magma', 'mako', 'nipy_spectral', 'ocean', 'pink', 'plasma', 'prism', 'rainbow', 'rocket', 'seismic', 'spring', 'summer', 'tab10', 'tab20', 'tab20b', 'tab20c', 'terrain', 'twilight', 'twilight_shifted', 'viridis', 'vlag', 'winter'. |


<br/><br/>

-----------------
## Package history and releases
A list of all releases and respective description of changes applied could be found [here](https://github.com/sebastian-gregoricchio/ChIP_Zwart/blob/main/NEWS.md).

## Contact
For any suggestion, bug fixing, commentary please report it in the [issues](https://github.com/sebastian-gregoricchio/ChIP_Zwart/issues)/[request](https://github.com/sebastian-gregoricchio/ChIP_Zwart/pulls) tab of this repository.

## License
This repository is under a [GNU General Public License (version 3)](https://github.com/sebastian-gregoricchio/ChIP_Zwart/blob/main/LICENSE.md/LICENSE.md).

<br/>

#### Contributors
![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/ChIP_Zwart?size=50&padding=5&bots=true)
