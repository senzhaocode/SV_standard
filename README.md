# *SV_standard* - aggregate and convert raw SV calls from mutiple samples

*SV_standard* is a perl script for aggregating SV calls from multiple samples and converting into an expected format for [**FuSViz**](https://github.com/senzhaocode/FuSViz) input. It merges raw SV calls from a range of tools ([*Manta*](https://github.com/Illumina/manta), [*Svaba*](https://github.com/walaj/svaba), [*Delly*](https://github.com/dellytools/delly) and [*Lumpy*](https://github.com/arq5x/lumpy-sv) for DNA-seq data; [*Dragen*](https://www.illumina.com/products/by-type/informatics-products/dragen-secondary-analysis.html), [*STAR-fusion*](https://github.com/STAR-Fusion/STAR-Fusion), [*Arriba*](https://github.com/suhrig/arriba), [*Fusioncatcher*](https://github.com/ndaniel/fusioncatcher) and [*deFuse*](https://github.com/amcpherson/defuse) for RNA-seq data) and convert them into a tab-separated values (TSV) format.

## Installation and usage

### 1. Prerequisites

-   [*SURVIVOR*](https://github.com/fritzsedlazeck/SURVIVOR)

    An installation of the tool *SURVIVOR* (version \>=*1.0.7*) is required to merge SVs (i.e., VCF files) from DNA-seq to generate a consensus or multi-sample vcf file. See [here](https://github.com/fritzsedlazeck/SURVIVOR) for an instruction of installation and usage.

    NOTE: *SV_standard* should be installed on a MacOS or Linux/UNIX operating system
    
