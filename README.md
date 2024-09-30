# Retrieve Human Protein Domains

Retrieving human protein domains from [UniProt](https://www.uniprot.org) based on Ensembl transcript ID.

Given a list of Ensembl transcript IDs (i.e. ENST IDs), we:
1. Retrieve the corresponding protein ID, gene name and ID, UniProt ID, and UniProt URL.
1. Retrieve the protein domains from UniProt.
1. Generate an excel file with the IDs and protein domains.


# Configuration File
`./src/config/config.toml`

# Input
A text file containing a list of Ensembl transcript IDs (text file name defined in `./src/config/config.toml`)

For example, a text file containing:
```python 
ENST00000288135
ENST00000302278
ENST00000559488
```

# Output
An excel file containing the corresponding protein domains (excel file name defined in `./src/config/config.toml`).

# Execution Flow
1. Set the configuration parameters in `./src/config/config.toml`
1. Run `./src/main.py`

# Requirements
1. Python >= 3.11
1. pandas
1. XlsxWriter
1. requests
1. From my [`Utils/`](https://github.com/yoramzarai/Utils) repository, the modules `toml_utils.py`, `uniprot_utils.py`, `rest_api_utils.py`, and `ensembl_rest_utils`

