# Joint k-mer and variant phasing

# Usage
First set of a configuration file. See `config/config.yaml` for an commented example. 
Then run snakemake with the following command pointing to your configuration file.
```
snakemake --profile profile/compute --configfile config/your_config.yaml
```

## Dependencies 
Dependencies are managed with conda. Be sure to include the following in your `.bashrc` if you want to use the pre-computed conda env. 
```
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
```