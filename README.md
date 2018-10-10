# OpenSWATH Snakemake Workflows

## Installation
- Either install all dependencies or use the preconfigured Docker container:

````
    docker pull openswath/develop:latest
````

- Clone the repository to your local installation.

## Usage
- Start up an instance of the Docker container:
````
    docker run --name osw_wf --rm -v workflows:/data -i -t openswath/develop
````

- Copy your DIA files to ``dia_data``.
- Copy your DDA files to ``dda_data``.
- Edit the parameters in ``params`` if necessary.
- Execute the full workflow:
````
    snakemake --snakefile Snakemake.diau -j4
    snakemake --snakefile Snakemake.library -j4
    snakemake --snakefile Snakemake.openswath -j4
````
