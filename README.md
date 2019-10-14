# datascience_tfm

Currently work have been presented as the final project to DataScience master degree.


##### Table of contents
* [Manual pages](#manual)
* [Some considerations before start](#before)
* [Installation and requirements](#install)
* [Analisys of the relation between variables](#relation)
* [Clustering](#cluster)
* [Final prioritization](#prioritization)

---
<a name="manual"></a>
### Manual pages

Documentation for the project utilities is abailable by cloning this repository.

---
<a name="install"></a>
### Installation and requirements
The project requires:

  * **Linux system** (>= Ubuntu 18.04 LTS system recommended; no other distributions were tested)
  * **Python** (>=3.7 recommended). Python languaje can be easily installed throug the [anaconda environment](https://www.anaconda.com/distribution/).

The remaining dependencies can be installed using the included `setup.sh` script. 
Basic instructions:
```bash
git clone https://github.com/nacasfer/datascience_tfm
cd datascience_tfm
bash setup.sh
```
The installer may also be used to check for updates to this and co-dependent packages.

---
<a name="before"></a>
## Some considerations before start

Currently project starts with a set of genomic variants. These variants where collected by our group and no distribution is available. Random example using [1000genomes] (https://)  phase 3 database is provided for script testing.
Aligment, varint calling and variant VCF annotation steps were perform with following modules using UCSC [HG37] fasta reference and a proper BED file:

  1. Mapping step was perform using [BWA] (http://) software
  2. Variant Call step was perform folllowing best practices of [GATK](http://) through `--HaploTypeCaller` module.
  3. Variant annotation step was perform by joint output from two annotation workflows: [Annovar] (https://) and [IonReporter](https://) web service.
  4. Final datasets were obtained by joint of manualy-anotated tsv patient files for the dissease variant candidate. 


---
<a name="relation"></a>
## Analysis of the relation between variables


---
<a name="clustering"></a>
## Clustering


---

<a name="prioritization"></a>
## Final prioritization

