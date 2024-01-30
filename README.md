 # *SV_standard* - aggregate and convert raw SV calls from mutiple samples

*SV_standard* is a perl script for aggregating SV calls from multiple samples and converting into an expected format for [**FuSViz**](https://github.com/senzhaocode/FuSViz) input. It merges raw SV calls from a range of tools ([*Manta*](https://github.com/Illumina/manta), [*Svaba*](https://github.com/walaj/svaba), [*Delly*](https://github.com/dellytools/delly) and [*Lumpy*](https://github.com/arq5x/lumpy-sv) for DNA-seq data; [*Dragen*](https://www.illumina.com/products/by-type/informatics-products/dragen-secondary-analysis.html), [*STAR-fusion*](https://github.com/STAR-Fusion/STAR-Fusion), [*Arriba*](https://github.com/suhrig/arriba), [*Fusioncatcher*](https://github.com/ndaniel/fusioncatcher) and [*deFuse*](https://github.com/amcpherson/defuse) for RNA-seq data) and convert them into a tab-separated values (TSV) format.

## Quickstart

### 1. Prerequisites

-   [*SURVIVOR*](https://github.com/fritzsedlazeck/SURVIVOR)

    An installation of the tool *SURVIVOR* (version \>=*1.0.7*) is required to merge SVs (i.e., VCF files) from DNA-seq data to generate a consensus or multi-sample file. See [here](https://github.com/fritzsedlazeck/SURVIVOR) for an instruction of installation and usage.

    **NOTE:** *SV_standard* should be installed on a MacOS or Linux/UNIX operating system


### 2. Installation

-   Download

    ```
    wget wget https://github.com/senzhaocode/SV_standard/releases/download/v1.0/SV_standard_v1.0.tar.gz`
    unzip SV_standard_v1.0.tar.gz && mv SV_standard_v1.0 SV_standard
    ```

-   Set environmental variables

    ```
    PERL5LIB="$PERL5LIB:/where_is_path/SV_standard/lib"
    export PERL5LIB
    ```
    
    **NOTE:** we recommend that users add the path of perl libraries to .bashrc, then `source .bashrc`

### 3. Run an example to aggregate SVs called from DNA-seq data

    perl SV_standard.pl --genome hg38 --type DNA --filter PASS \
        --anno anno \
        --input example/DNA/input \
        --output example/DNA/output

The *example/DNA/input* folder contains raw SVs VCF files called from one or several tools (e.g., *Manta*, *Svaba*, *Delly* and *Lumpy*) per sample. Users have to prepare for input files following the folder organization belew:

    example/DNA/input
               |--- T001 # sample name
                  |--- Manta.vcf 
                  |--- Svaba.vcf
               |--- T002 # sample name
                  |--- Delly.vcf
               |--- T003 # sample name
                  |--- Lumpy.vcf

**NOTE:** the raw VCF file (the compressed and indexed ones using bgzip and tabix are recommended) should be named using caller nomenclature. In terms of raw SVs called from *Svaba*, no SV types (e.g., BND, INV, DEL and DUP) are available. Before running SV_standard.pl, we provide an in-house R script (at script folder) to assign SV type to each call and convert original vcf file following FuSViz requirement. For the usage - `Rscript script/svaba_svtype.R svaba_raw.vcf svaba_new.vcf`

An example of *example/DNA/output* folder contains the results:

 - `Final_DNA_SVs.txt` (a tab-separated format file with aggregated DNA-seq raw SV calls used for FuSViz input)
 - `TXX_list` (an intermediate file lists callers used for raw SV calling)
 - `TXX_merge.bedpe` (an intermediate bedpe file contains a merge of raw SVs from multiple callers)
 - `TXX_merge.vcf` (an intermediate vcf file contains a merge of raw SVs from multiple callers)

### 4. Run an example to aggregate SVs called from RNA-seq data

    perl SV_standard.pl --genome hg38 --type RNA \
        --anno anno \
        --input example/RNA/input \
        --output RNA_output

The *example/RNA/input* folder contains raw SVs called from one or several tools (e.g., *Dragen*, *deFuse*, *STAR-fusion*, *Arriba* and *Fusioncatcher*) per sample. Users have to prepare for input files as the following organization below:

     example/RNA/input
                |--- T001 # sample name
                   |--- Arriba.tsv
                   |--- STAR-fusion.tsv
                |--- T002 # sample name
                   |--- Dragen.txt
                   |--- Fusioncatcher.txt

  **NOTE:** the raw input file should be a tab-separated format (TSV or TXT) file that is named using caller nomenclature.

  An example of *example/RNA/output* folder contains the results:

   - `Final_RNA_SVs.txt` (a tab-separated format file with aggregated RNA-seq raw SV calls used for FuSViz input)




    
