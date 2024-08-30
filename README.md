# MtbRf: A Pan-Lineage Reference Genome of _Mycobacterium tuberculosis_

Data and code for creating and analyzing MtbRf, a novel pan-lineage reference genome of _Mycobacterium tuberculosis_.

## About This Repository

This repository contains the data and source code associated with the research paper titled "Pan-Lineage _Mycobacterium tuberculosis_ Reference Genome for Enhanced Molecular Diagnosis". It includes the reference genome sequence (`mtbrf.fa`) and the corresponding annotations (`mtbrf.gff`) for MtbRf, which were generated using the method described in the paper and the NCBI Prokaryotic Genome Annotation Pipeline (PGAP), respectively.

The creation of MtbRf aims to broaden molecular genomic studies of _Mycobacterium tuberculosis_ (MTB) to include lineage-specific analysis, consequently filling the research gap between MTB lineages. Additionally, MtbRf provides yet unaccounted information associated with drug resistance, making it a critical resource for future research in this field.

## Navigation

- `mtbrf.fa`: FASTA file for the sequence of MtbRf.
- `mtbrf.gff`: GFF formatted annotations for MtbRf.
- `scripts`: Directory containing scripts for analysis pipelines.
- `data`: Directory for related data.

### Associated Links

The MtbRf reference genome is accessible with accession number `GCA_963525475` in public databases:
- EBI: https://www.ebi.ac.uk/ena/browser/view/GCA_963525475
- NCBI: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_963525475

The WGS data for this project is accessible with accession number `PRJEB66375`:
- EBI: https://www.ebi.ac.uk/ena/browser/view/PRJEB66375
- NCBI: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB66375

## Performance

### Median mapping rate by _Mycobacterium tuberculosis_ lineages

| **Lineage** | **MtbRf (GCA_963525475.1)** | **H37Rv-1 (GCA_026185275.1)** | **H37Rv (GCA_000195955.2)** |
|-------------|-----------------------------|-------------------------------|-----------------------------|
| _L1_        | 99.74%                      | 99.20%                        | 99.19%                      |
| _L2_        | 99.77%                      | 99.43%                        | 99.42%                      |
| _L3_        | 99.71%                      | 99.37%                        | 99.36%                      |
| _L4_        | 99.81%                      | 99.49%                        | 99.48%                      |
| All         | 99.73%                      | 99.34%                        | 99.33%                      |

## Citation

Kunhyung Bahk, Joohon Sung, Mitsuko Seki, Kyungjong Kim, Jina Kim, Hongjo Choi, Jake Whang, Satoshi Mitarai, Pan-lineage *Mycobacterium tuberculosis* reference genome for enhanced molecular diagnosis, *DNA Research*, Volume 31, Issue 4, August 2024, dsae023, https://doi.org/10.1093/dnares/dsae023
