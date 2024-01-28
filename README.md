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

### 3. Run an example to aggregate SVs from DNA-seq data

    perl SV_standard.pl --genome hg38 --type DNA --filter PASS \
        --anno anno \
        --input example/DNA/input \
        --output DNA_output

The *example/DNA/input* folder contains raw SVs VCF files called from a range of tools per sample. Users have to prepare for their input files following the folder structure belew:

    example/DNA/input
        |--- T001 (sample name)
            |--- Manta.vcf 
            |--- Svaba.vcf
        |--- T002 (sample name)
            |--- Delly.vcf
        |--- T003 (sample name)
            |--- Lumpy.vcf

**NOTE:** the raw VCF file (the compressed and indexed ones using bgzip and tabix are recommended) should be named using caller nomenclature.
The *DNA_input* folder contains the results:

### 4. Run an example to aggregate SVs from RNA-seq data

    perl SV_standard.pl --genome hg38 --type RNA --filter PASS \
        --anno anno \
        --input example/RNA/input \
        --output RNA_output


    
