# TCR-BCR-nf
## Pipeline to genotype Tcell and Bcell receptors from bulk or single-cell RNA-seq data

## Description
This pipeline calls software TRUST4 to genotype T cell and B cell receptors.

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- [TRUST4](https://github.com/liulab-dfci/TRUST4)

You can avoid installing all the external software by only installing Docker and using the official TRUST4 docker container (at quay.io/biocontainers/trust4:<tag>, see tags at https://quay.io/repository/biocontainers/trust4?tab=tags). See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information. Note that TRUST4 requires references files IMGT+C.fa and bcrtcr.fa, that can be generated using perl scripts from the trust4 git repository:
```
perl BuildImgtAnnot.pl Homo_sapien > IMGT+C.fa
grep ">" IMGT+C.fa | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq > bcr_tcr_gene_name.txt  
perl BuildDatabaseFa.pl reference.fa annotation.gtf bcr_tcr_gene_name.txt > bcrtcr.fa
```

## Input
  | Type      | Description     |
  |-----------|---------------|
  | input_folder    | Folder with bam files to be processed |

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --IMGTC_fasta    | IMGT+C.fa | fasta file of reference genome from the international ImMunoGeneTics information system [file human_IMGT+C.fa from https://github.com/liulab-dfci/TRUST4] |
| --bcrtcr_fasta    | bcrtcr.fa | fasta file with reference BCR and TCR regions [see https://github.com/liulab-dfci/TRUST4 for generation] |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --barcode   |  None | Run trust4 using a specific barcode (for single cell data only, usually BC for CellRanger output) |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |

## Usage
  ```
  nextflow run iarcbioinfo/TCR-BCR-nf --input_folder bams --IMGTC_fasta TRUST4/IMGT+C.fa --bcrtcr_fasta TRUST4/bcrtcr.fa --barcode BC
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | TRUST_*_barcode_report.tsv    | Report with inferred cell type and TCR and BCR genotypes |
  | TRUST_*_report.tsv | Report with proportion of each TCR and BCR genotype |


## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Nicolas Alcala*    | alcalan@iarc.who.int | Developer to contact for support |

## References
Song, L., Cohen, D., Ouyang, Z. et al. TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nat Methods 18, 627–630 (2021). https://doi.org/10.1038/s41592-021-01142-2
