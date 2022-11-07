# datascience_tfm

Currently work have been presented as the first approach to a future propperly DataScience model.


##### Table of contents
* [Manual pages](#manual)
* [Some considerations before start](#before)
* [Installation and requirements](#install)
* [Filling NaN's: Analisys of the relation between variables](#nan)
* [Clustering](#cluster)
* [Final considerations](#final)

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

Currently project starts with a set of genomic variants. These variants where collected by our group and no distribution is available. Random example using [1000genomes phase 3](https://www.internationalgenome.org/category/phase-3/) database is provided for script testing.
Aligment, varint calling and variant VCF annotation steps were perform with following modules using UCSC [HG37] fasta reference and a proper BED file:

  1. Mapping step was perform using [BWA](http://bio-bwa.sourceforge.net/) software
  2. Variant Call step was perform folllowing best practices of [GATK](https://software.broadinstitute.org/gatk/best-practices/) through `--HaploTypeCaller` module.
  3. Variant annotation step was perform by joint output from two annotation workflows: [Annovar](http://annovar.openbioinformatics.org/en/latest/) and [IonReporter](https://ionreporter.thermofisher.com/ir/) web service.
  4. Final datasets were obtained by joint of manualy-anotated tsv patient files for the dissease variant candidate. 


---
<a name="nan"></a>
## Filling NaN's: Analysis of the relation between variables
Filling NaN values is strongly necessary for corret predictions. In this case, NaN values of one feature are filled throug linear, logaritmic and logistic regresion with other features. These models are availables in the "MLmodel" folder, downloadble by cloning this repository. Models are automaticaly appied by using `predict.py` script.

---
<a name="clustering"></a>
## Clustering
Several models are avaiable at `predict.py` script. For not using all just coment de non-selected ones in the code. Otherwise, all the predictions will be storage on a new folder with the name of the model and the name of the sample.

---

<a name="final"></a>
## Final considerations
Currently model is ready to use throug the instruction below. Just change the path of the folder to make the predictions, or crate one as instructions. 

```bash
python3 predict.py
```
If you wish to train models with your own data, just create a main dataset using the  `TrainDataFormat.py` script and run the  `main.py` script as follows. Remember set the correct path.
```bash
python3 main.py
```
For validation, set the correct path to not-trained data and run the `Comp_w_Test.py` as follows.
```bash
python3 Comp_w_Test.py
```
