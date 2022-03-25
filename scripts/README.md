# Folder organization

- `library_construction/`: scripts for the TE annotation:
    - `octfta_to_fasta.sh`: script that allows to parse OneCodeToFindThemAll output and recover a fasta file containing sequences from TEs full copies.
        - `octfta_to_bed.py`: called by `octfta_to_fasta.sh` to convert the `.tsv` in a `.bed`.
    - `parse_cdhit.py`: parse cd-hit-est results. It recovers FASTA sequences used to call cd-hit-est for each cluster. Then it separates singleton clusters in 'singleton.fa' and it calls a consensus for all other clusters using Refiner (from RepeatModeler2). Consensi sequences can founded in 'consensi.fa' and alignements in 'consensi.stk'. The 'stats.txt' reports general statitics on clusters.
- `genotype/`: scripts for the TE insertions calling.