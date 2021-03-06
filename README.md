# Kunkel mutagenesis for Bagel project, Siegel Lab UCD 

## Authors

+ Yin He, Transcriptic
+ Ben Miles, Transcriptic
+ Alex Carlin, University of California, Davis 

## Use

To run the protocol, first install the Autoprotocol package and the Transcriptic command line utility, 
and a preview of the protocol. You must be in the repo working directory for this to work. 

```bash
pip install autoprotocol transcriptic
transcriptic preview Kunkel 
```

If you create a new version of the protocol:

```bash 
transcriptic build-release Kunkel 
transcriptic upload-release Kunkel <my_project_name>
```

(see `transcriptic --help` for more). Then, follow the link to publish the protocol on Transcriptic's website. 

Then, you can launch runs: Launch a Run > Organization > Kunkel. The run requires you enter 

1. an aliquot of ssDNA 
2. antibiotic resistance 
3. CSV of the mutants you wish to order (see `make_csv.py` in this repo) 

