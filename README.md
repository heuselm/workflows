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
    docker run --name osw --rm -v ${PWD}:/data -i -t openswath/develop:latest
````

- Copy your DIA files to ``data_dia``.
- Copy your DDA files to ``data_dda``.
- Edit the parameters in ``params`` if necessary.
- Edit the parameters in ``Snakefile.library`` if necessary.
- Edit the parameters in ``Snakefile.openswath`` if necessary.
- Execute the full workflow:
````
    snakemake --snakefile Snakefile.diau -j4
    snakemake --snakefile Snakefile.library -j4
    snakemake --snakefile Snakefile.openswath -j4
````
