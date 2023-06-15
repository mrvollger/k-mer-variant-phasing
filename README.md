# Joint k-mer and variant phasing

This is a pipeline designed to phase variants and PacBio Hifi data using a combination of k-mers and variants. 

However, by running it in different modes you can also use it to phase using only variants and not k-mers (hiphase), or to just generate variant calls without any phasing of reads or the vcf (DeepVariant).

## Usage

First set of a configuration file. See `config/config.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
```
snakemake --profile profile/compute --configfile config/your_config.yaml
```

## Running without parental data

You can run the pipeline the same way, it will just only use `hiphase` instead of a combination of `hiphase` and k-mers to phase.
See `config/no_parental.yaml` for an example input file.

## Running only DeepVariant

First set of a configuration file. See `config/deepvariant.yaml` for a commented example. 
Then run snakemake with the following command pointing to your configuration file.
``` 
snakemake --profile profile/compute --configfile config/your_config.yaml deepvariant
```

### Dependencies 

Dependencies are managed with conda. Be sure to include the following in your `.bashrc` if you want to use the pre-computed conda env. 
```
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
```