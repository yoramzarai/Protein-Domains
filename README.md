# Retrieve Protein Domains

Given a list of Ensembl transcript IDs (e.g. ENST IDs), we:
1. Retrieve the corresponding protein ID, gene name and ID, UniProt ID, and UniProt URL.
1. Retrieve the protein domains from UniProt.
1. Generate an excel file with the IDs and protein domains.


# Configuration File
`./config/config.toml`

# Input
A text file containing a list of Ensembl transcript IDs (file name defined in `./config/config.toml`)

For example, a `transcript_list.txt` file containing:
```python 
ENST00000288135
ENST00000302278
ENST00000559488
```

# Output
An excel file containing the corresponding protein domains (file name defined in `./config/config.toml`).

# Execution Flow
1. Set the required parameters in `./config/config.toml`
1. Run `main.py`

# Requirements
1. Python >= 3.11
1. pandas
1. XlsxWriter
1. requests
1. From my `Cancer_mut/Utils/` folder, the modules `toml_utils.py`, `uniprot_utils.py`, `rest_api_utils.py`, `ensembl_rest_utils`

