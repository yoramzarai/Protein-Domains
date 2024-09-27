# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments,too-many-locals
"""
Utils for main.py.
"""
import pathlib
from dataclasses import dataclass
import pandas as pd
import toml_utils as tmut
import uniprot_utils as uput  # in my Utils/ folder
import ensembl_rest_utils as erut  # in my Utils/ folder


# configuration Toml file
Cnfg_Toml_file: pathlib.Path = pathlib.Path('./src/config/config.toml')

# Uniprot URL of a UniProt ID DUMMYID (DUMMYID will be replaced by a valid UniProt UD)
Uniprot_url_template: str = "https://www.uniprot.org/uniprotkb/DUMMYID/entry"


@dataclass
class Labels:
    """Labels"""
    Transcript_ID: str = 'Transcript_ID'
    Protein_ID: str = 'Protein_ID'
    Gene_name: str = 'Gene_name'
    Gene_ID: str = 'Gene_ID'
    UniProt_ID: str = 'UniProt_ID'
    UniProt_URL: str = 'UniProt_URL'
    Domains: str = 'Domains'

def load_config() -> dict:
    """Loads the Toml configuration file."""
    cnfg_data = tmut.myToml().load(Cnfg_Toml_file)
    # verify configuration validity
    check_configuration(cnfg_data)
    return cnfg_data

def print_config(cnfg_data: dict) -> None:
    """Pretty print of config data."""
    tmut.print_nested_dicts(cnfg_data)

def check_configuration(cnfg_data: dict) -> None:
    """Check configuration validity."""
    assert cnfg_data['Assembly']['version'] in ["GRCh37", "GRCh38"], f"Aeembly version {cnfg_data['Assembly']['version']} not supported !!"
    assert cnfg_data['Output']['format'] in ["basic", "compact", "expanded"], f"Output format {cnfg_data['Output']['format']} not supported !!"
    if not pathlib.Path(cnfg_data['Transcript']['file']).is_file():
        raise FileNotFoundError(f"Cannot find input transcript file {cnfg_data['Transcript']['file']} !!")
    if not isinstance(cnfg_data['Domains']['uniprot_features'], list):
        raise TypeError(f"Domain:features in {Cnfg_Toml_file} must contain a list of UniProt domains !!")

# this function was taken from myutils.py
def dfs_to_excel_file(dfs: list[pd.DataFrame], excel_file_name: str, sheet_names: list[str],
                      add_index: bool = False,
                      na_rep: str = 'NaN',
                      float_format: str | None = None,
                      extra_width: int = 0,
                      header_format: dict | None = None) -> None:
    """
    Write DataFrames to excel, while auto adjusting column widths based on the data.
    """
    if len(dfs) != len(sheet_names):
        raise ValueError("dfs_to_excel_file: the numbers of dfs and sheet names must match !!")

    with  pd.ExcelWriter(excel_file_name) as writer:
        for sheet_name, df in zip(sheet_names, dfs):
            df.to_excel(writer, sheet_name=sheet_name, index=add_index, na_rep=na_rep, float_format=float_format)

            # Auto-adjust columns' width
            worksheet = writer.sheets[sheet_name]
            for column in df.columns:
                column_width = max(df[df[column].notna()][column].astype(str).map(len).max(), len(column)) + extra_width  # this first removes NA values from the column
                col_idx = df.columns.get_loc(column)
                worksheet.set_column(col_idx, col_idx, column_width)

            if header_format is not None:
                h_format = writer.book.add_format(header_format)
                for col_num, value in enumerate(df.columns.values):
                    worksheet.write(0, col_num, value, h_format)


def get_transcripts(cnfg_data: dict) -> list[str]:
    """Reads the transcripts input file and returns the transcript IDs."""
    try:
        with open(cnfg_data['Transcript']['file'], 'rt', encoding='UTF-8') as fp:
            return [x for x in [line.rstrip() for line in fp] if 'ENST' in x]
    except FileNotFoundError:
        print(f"Can not find input transcripts file {cnfg_data['Transcript']['transcript_file']}. Please check configuration file, under ['Transcript']['transcript_file'] !!")
        raise

def get_transcripts_IDs(cnfg_data: dict, transcripts: list[str]) -> dict[str,dict[str,str]]:
    """
    Retreives different IDs of a transcript.
    """
    rapi = erut.REST_API(cnfg_data['Assembly']['version'])
    info: dict = {}
    for transcript in transcripts:
        # Gene ID and name
        if cnfg_data['IDs']['get_gene_name'] or cnfg_data['IDs']['get_gene_id']:
            ensg_id = rapi.get_transcript_parent(transcript)
            if cnfg_data['IDs']['get_gene_name'] :
                ensg_name = '' if ensg_id == '' else rapi.ENSG_id2symbol(ensg_id)

        uniprot_id = uput.ensembl_id2uniprot_id(transcript)
        info[transcript] = {
            Labels.Protein_ID: rapi.transcript_id2protein_id(transcript) if cnfg_data['IDs']['get_protein_id'] else '',
            Labels.Gene_ID: ensg_id if cnfg_data['IDs']['get_gene_id'] else '',
            Labels.Gene_name: ensg_name if cnfg_data['IDs']['get_gene_name'] else '',
            Labels.UniProt_ID: uniprot_id,
            Labels.UniProt_URL: Uniprot_url_template.replace('DUMMYID', uniprot_id) if cnfg_data['IDs']['get_uniprot_url'] else ''
        }
    return info

def get_uniprot_domains(cnfg_data: dict, transcripts_ids: dict[str,dict[str,str]]) -> dict[str,dict]:
    """Get UniProt domains given the UniProt ID."""
    info: dict = {}
    features: list = cnfg_data['Domains']['uniprot_features']
    for transcript, transcript_ids in transcripts_ids.items():
        if (df_uniprot := uput.retrieve_protein_data_features_subset(transcript_ids[Labels.UniProt_ID], features)).empty:
            print(f"\n[** No {features} UniProt features were found for {transcript} (UniProt ID={transcript_ids[Labels.UniProt_ID]}) **]\n")
        info[transcript] = {
            'domains_df': df_uniprot,  # in a dataframe format
            'domains_list': list(df_uniprot.T.to_dict().values())  # in a list of domains format
        } | transcript_ids
    return info

def _gen_basic_domain_dataframe(cnfg_data: dict, transcripts_domains: dict[str,dict]) -> pd.DataFrame:
    """Generate a dataframe with all transcripts, where each domain is listed in a separate row."""
    all_dfs: list = []
    for k, v in transcripts_domains.items():
        df: pd.DataFrame = v['domains_df'].copy()
        df.insert(0, Labels.Transcript_ID, k)
        df.insert(1, Labels.UniProt_ID, v[Labels.UniProt_ID])
        i = 2
        if cnfg_data['IDs']['get_gene_id']:
            df.insert(i, Labels.Gene_ID, v[Labels.Gene_ID])
            i += 1
        if cnfg_data['IDs']['get_gene_name']:
            df.insert(i, Labels.Gene_name, v[Labels.Gene_name])
            i += 1
        if cnfg_data['IDs']['get_protein_id']:
            df.insert(i, Labels.Protein_ID, v[Labels.Protein_ID])
            i += 1
        if cnfg_data['IDs']['get_uniprot_url']:
            df.insert(i, Labels.UniProt_URL, v[Labels.UniProt_URL])

        all_dfs.append(df)
    return pd.concat(all_dfs).reset_index(drop=True)


def _gen_compact_domain_dataframe(cnfg_data: dict, transcripts_domains: dict[str,dict]) -> pd.DataFrame:
    """Generate a dataframe with all transcripts, where all domains of a transcript are aggregate into a single row."""
    all_dfs: list = []
    for k, v in transcripts_domains.items():
        df = pd.DataFrame({
            Labels.Transcript_ID: k,
            Labels.UniProt_ID: v[Labels.UniProt_ID]
        }, index=[0])
        i = 2
        if cnfg_data['IDs']['get_gene_id']:
            df.insert(i, Labels.Gene_ID, v[Labels.Gene_ID])
            i += 1
        if cnfg_data['IDs']['get_gene_name']:
            df.insert(i, Labels.Gene_name, v[Labels.Gene_name])
            i += 1
        if cnfg_data['IDs']['get_protein_id']:
            df.insert(i, Labels.Protein_ID, v[Labels.Protein_ID])
            i += 1
        if cnfg_data['IDs']['get_uniprot_url']:
            df.insert(i, Labels.UniProt_URL, v[Labels.UniProt_URL])
        df[Labels.Domains] = "|".join([",".join([f"{k}:{v}" for k, v in x.items()]) for x in v['domains_list']])
        all_dfs.append(df)
    return pd.concat(all_dfs).reset_index(drop=True)

# def _gen_compact_domain_dataframe(cnfg_data: dict, transcripts_domains: dict[str,dict]) -> pd.DataFrame:
#     """Generate a dataframe with all transcripts, where all domains of a transcript are aggregate into a single row."""
#     all_dfs: list = []
#     for k, v in transcripts_domains.items():
#         df = pd.DataFrame({
#             Labels.Transcript_ID: k,
#             Labels.UniProt_ID: v[Labels.UniProt_ID],
#             Labels.Gene_ID: v[Labels.Gene_ID],
#             Labels.Gene_name: v[Labels.Gene_name],
#             Labels.Protein_ID: v[Labels.Protein_ID],
#             Labels.UniProt_URL: v[Labels.UniProt_URL]
#         }, index=[0])
#         df[Labels.Domains] = "|".join([",".join([f"{k}:{v}" for k, v in x.items()]) for x in v['domains_list']])
#         drop_cols = ['' if cnfg_data['IDs']['get_gene_id'] else Labels.Gene_ID ] + ['' if cnfg_data['IDs']['get_gene_name'] else Labels.Gene_name] +\
#         ['' if cnfg_data['IDs']['get_protein_id'] else Labels.Protein_ID] + ['' if cnfg_data['IDs']['get_uniprot_url'] else Labels.UniProt_URL]
#         df.drop(columns=[x for x in drop_cols if x != ''], inplace=True)
#         all_dfs.append(df)
#     return pd.concat(all_dfs).reset_index(drop=True)

def generate_output_table(cnfg_data: dict, transcripts_domains: dict[str,dict]) -> tuple[list[pd.DataFrame],list[str]]:
    """Returns the output dataframe containing IDs and domains."""
    match cnfg_data['Output']['format']:
        case 'basic':
            dfs = [_gen_basic_domain_dataframe(cnfg_data, transcripts_domains)]
            transcript_IDs = [Labels.Domains]  # when all transcripts are in the same sheet, we simply call the sheet Labels.Domains]
        case 'compact':
            dfs = [_gen_compact_domain_dataframe(cnfg_data, transcripts_domains)]
            transcript_IDs = [Labels.Domains]  # when all transcripts are in the same sheet, we simply call the sheet Labels.Domains]
        case 'expanded':
            df = _gen_basic_domain_dataframe(cnfg_data, transcripts_domains)
            transcript_IDs, dfs = zip(*list(df.groupby(by=Labels.Transcript_ID)))  # sheet name is the transcript ID
        case _:
            raise ValueError(f"Output format {cnfg_data['Output']['format']} is no supported. Please check configuration file, under ['Output']['format'] !!")
    return dfs, transcript_IDs


def generate_output_file(cnfg_data: dict, transcripts_domains: dict[str,dict]) -> None:
    """Generates the ouput file containing the IDs and domains"""
    dfs, sheet_names = generate_output_table(cnfg_data, transcripts_domains)
    dfs_to_excel_file(dfs, cnfg_data['Output']['file'], sheet_names=sheet_names, add_index=False, extra_width=2)
