## Capsule Characterization for *Neisseria meningitidis*

**Dependencies**
```
Python 3.0 or greater (with urllib3 and Biopython)
ncbi-blast+ 2.3.0 or greater
```
How to use:

#### 1) Generate local blast database from the Neisseria capsule genes and insertion element collection:

```python
python3 build_neisseria_dbs.py
```

It is recommended that this script be run in the same directory that "custom_allele_sets" lives in to properly build the collection.

This will produce a directory called "neisseria_capsule_DB," containing all of the blast databases used in the main script.


#### 2) Use characterization script to obtain serogroup information for a collection of FASTA genome assembly files

Usage:

```python

characterize_neisseria_capsule.py [-h] -d INDIR -o OUT [-t THREADS]

Script for predicting serogroup of Neisseria genomes

optional arguments:
  -h, --help            show this help message and exit
  -d INDIR, --indir INDIR
                        Input Dir
  -o OUT, --out OUT     Output Directory
  -t THREADS, --threads THREADS
                        Number of Threads to use (default=1)

```

Example usage:

```python
python3 characterize_neisseria_capsule.py -d my_fasta_file_directory -t 10 -o results
```

This example will take a directory of fasta files (my_fasta_file_directory), and will process them using 10 CPU threads. The results will be found in the output folder "results."

#### 3) Interpreting Output:

Three output directories will be created:

* **GFF:** This directory contains General Feature Format (GFF) files, containing annotations of the genes and insertion elements found in each FASTA assembly. 
* **JSON:** This directory contains intermediary files consisting of a "raw" and "final" JSON for each input assembly. The "raw" json files contain the raw blast results without filtering for identifying the best match for each region.
The final JSON file contains the filtered blast results and are used to generate the serogroup predictions.
* **Serogroup:** This directory contains two files: A tab delimited serogroup file consisting of the assembly name, predicted serogroup, capsule genes identified, and notes containing the serogroup backbone and mutations identified that may prevent capsule expression. The other file is a JSON form of this same data.


 