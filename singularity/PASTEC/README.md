# PASTEC-singularity

## Description

`Singularity` container for the transposable elements classification tool `PASTEC`(from the package [`REPET`](https://urgi.versailles.inra.fr/Tools/REPET)).


## Dependancies

- `singularity` (tested with v.3.9.5)

## Installation

### Development installation

- Clone this repo:
```
git clone https://github.com/TommasoBarberis/TE_Aalb.git
```

- Move to the `PASTEC` directory:
```
cd TE_Aalb/singularity/PASTEC
```

- Build the `.sif` image from the definition file:
```
sudo singularity build PASTEC.sif PASTEC.def
```

### Only user installation

- Pull directly the `.sif` image from the [sylabs.io](sylabs.io) registry:
```
singularity pull library://tommasobarberis98/pastec/pastec:1.0.0
```


## Utilisation

```
PROJECT_DIR=path/to/project/dir
```
This directory has to contain:
- the `fasta` file with consensi sequences: 

    From the `PASTEC` documentation, each sequence must have 60 bp (or less) by line, so if you have a "one-line" `FASTA` file, you can transform it with the following command:
    ```
    sed '/^>/!s/.\{60\}/&\n/g' oneline.fasta > multiline.fasta
    ```

    About the sequence headers, it is highly advised to write them like this : ">XX_i" with XX standing for letters and i standing for numbers. 

- the `PASTEClassifier.cfg` file (provided with this git repo) to update some parameters such as the path to a known database of transposable elements (see more here: [PASTEClassfier-tuto](https://urgi.versailles.inra.fr/Tools/PASTEClassifier/PASTEClassifier-tuto)). Options that you are suggested to update:
    - `project_name`: whatever you want.
    - `TE_nucl_bank`: `/mnt/nucl_bank.fa` such as `RepBase20.05_REPET.embl/repbase20.05_ntSeq_cleaned_TE.fa` (require subscription to [girinst](https://www.girinst.org/repbase/)), but you are free to use any other database. Then you can set:
        - `TE_BLRtx`: yes :arrow_right: homology with known TEs using tblastx.
        - `TE_BLRn`: yes :arrow_right: detection of helitron extremities.
    - `TE_prot_bank`: `/mnt/prot_bank.fa` such as `RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa` (require subscription to [girinst](https://www.girinst.org/repbase/)), but you are free to use any other database. Then you can set:
        - `TE_BLRx`: yes :arrow_right: homology with known TEs using blastx.
    - `HG_nucl_bank`: `/mnt/host_genes.fa` for the cDNA database of the host genes. Then you can set:
        - `HG_BLRn`: yes :arrow_right: homology with host genes.
    - `rDNA_bank`: `/mnt/rdna.fa` for the rDNA database. Then you can set:
        - `rDNA_BLRn`: yes :arrow_right: homology with rDNA.
    - `TE_HMM_profiles`: is already set to the `ProfilesBankForREPET_Pfam32.0.hmm` database, but your are free to use an another database.
    - Adjustable parameters in `[classif_consensus]` section (see [PASTEClassfier-tuto](https://urgi.versailles.inra.fr/Tools/PASTEClassifier/PASTEClassifier-tuto)). 

- the various databases specified in the `PASTEClassifier.cfg` file.


### Interactive mode

- Start a shell through the container:
```
singularity shell --bind $PROJECT_DIR:/mnt pastec_1.0.0.sif
```
- Once in the container, run:
```
python2.7 /opt/PASTEC_linux-x64-2.0/bin/PASTEClassifier.py -i /mnt/consensi.fasta -C /mnt/PASTEClassifier.cfg -p $N_CPU
```

### Batch mode

Create a file `pastec_cmd.sh` that will contain thePASTEC command in your project directory (`$PROJECT_DIR`):
```
#! /bin/bash 

N_CPU=1 #number of CPUs that you want use to run PASTEC

# in the following command the /mnt directory replace your project directory (through the --bind flag in the singularity command)
python2.7 /opt/PASTEC_linux-x64-2.0/bin/PASTEClassifier.py -i /mnt/consensi.fasta -C /mnt/PASTEClassifier.cfg -p $N_CPU
```

Then you can run `PASTEC` using the container as follow:
```
singularity exec --bind $PROJECT_DIR:/mnt PASTEC.sif /mnt/pastec_cmd.sh
```

### Other `PASTEC` options

- classification rules file name (e.g. PASTEClassifierRules.yaml) [optional]
 > `-D DECISIONRULESFILENAME`, `--decisionRules=DECISIONRULESFILENAME` 
- step (0/1/2): default: 0 for all steps, step 1 for detect features, step 2 for classification run tool in parallel
> `-S STEP`, `--step=STEP`
- clean temporary files [optional] [default: False]
> `-c`, `--clean`
- verbosity [optional] [default: 3, from 1 to 4]
> `-v VERBOSITY`, `--verbosity=VERBOSITY`