# MtbRf: A Pan-Lineage Reference Genome of _Mycobacterium tuberculosis_

Data and code for creating and analyzing MtbRf, a novel pan-lineage reference genome of _Mycobacterium tuberculosis_.

## About This Repository
This repository contains the data and source code associated with our research paper titled "MtbRf: A Pan-Lineage Reference Genome of _Mycobacterium tuberculosis_". It includes the reference genome sequence (`mtbrf.fa`) and the corresponding annotations (`mtbrf.gff`) for MtbRf, which were generated using the method described in the paper and the NCBI Prokaryotic Genome Annotation Pipeline (PGAP), respectively.

The sequence and annotation files represent the foundation of our work and will be continuously updated as our research progresses. Therefore, this repository should be considered a preliminary version.

The creation of MtbRf aims to broaden molecular genomic studies of Mycobacterium tuberculosis (MTB) to include lineage-specific analysis, consequently filling the research gap between MTB lineages. Additionally, MtbRf provides yet unaccounted information associated with drug resistance, hence, it is a critical resource for future research in this field.

## Repository Contents
- `mtbrf.fa`: This FASTA file contains the reference genome sequence of MtbRf.
- `mtbrf.gff`: This GFF file contains the annotations for the MtbRf reference genome.

## Performance
### Median mapping rate by *Mycobacterium tuberculosis* lineages
| **Lineage** | **H37Rv (GCA_000195955.2)** | **MtbRf (GCA_963525475.1)** |
|-------------|---------------------------|-------------|
| _L1_        | 99.19%                    | 99.74%      |
| _L2_        | 99.42%                    | 99.77%      |
| _L3_        | 99.36%                    | 99.71%      |
| _L4_        | 99.48%                    | 99.81%      |
| All         | 99.33%                    | 99.73%      |

## Future Updates
This repository is currently a preliminary version and will be updated continuously with improvements to the MtbRf genome and additional analysis code.

We welcome your interest and contributions to this ongoing project.

## Reference
The MtbRf genome and the associated data and code in this repository are part of a research project that is currently under review. Once the paper is accepted and published, we will update this section with full citation details.
