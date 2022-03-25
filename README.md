# Transposable elements of _Aedes albopictus_

## Managing the conda environment
To install the conda env:
```
conda env create -f TE_Aalb.yml
```

To activate it:
```
conda activate TE_Aalb
```

Create a new `.yml` file:
```
conda env export --from-history > TE_Aalb.yml
```

To update it:
```
conda env update --file TE_Aalb.yml --prune
```

## Organization of the repository

- `misc/`: old scripts and other stuff;
- `notebooks/`: jupyter notebooks and related scripts;
- `scripts/`: scripts used to produce results. 
