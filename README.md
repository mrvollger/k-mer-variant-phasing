# Joint k-mer and variant phasing

This is a pipeline designed to phase variants and PacBio Hifi data using a combination of k-mers and variants. 

However, by running it in different modes you can also use it to phase using only variants and not k-mers (hiphase), or to just generate variant calls without any phasing of reads or the vcf (DeepVariant).

## Usage

First set up a configuration file. See `config/config.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
```
snakemake --profile profile/compute --configfile config/your_config.yaml
```

## Running without parental data

You can run the pipeline the same way, it will only use `hiphase` instead of a combination of `hiphase` and k-mers to phase.
See `config/no_parental.yaml` for an example input file.

## Running only DeepVariant

First set up a configuration file. See `config/deepvariant.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
``` 
snakemake deepvariant --profile profile/compute --configfile config/your_config.yaml 
```

## Running only alignment (pbmm2)

First set up a configuration file. See `config/align.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
``` 
snakemake pbmm2 --profile profile/compute --configfile config/your_config.yaml 
```


### Dependencies 
Requirements for executing the pipeline are:
```
mamba
snakemake>=7.32.0
python<=3.11
```
**Note:**
Snakemake is currently not compatible with Python >=3.12 due to a change in f-strings


Additional dependencies are managed with conda. Be sure to include the following in your `.bashrc` if you want to use the pre-computed conda env. 
```
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
```