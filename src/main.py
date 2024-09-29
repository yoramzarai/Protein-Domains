# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments
"""
Main script for retrieving human protein domains.

1. Set the parameters in ./config/config.toml
2. Run this script.
"""
import Utils.utils as u

def main() -> None:
    """main function."""
    # load configuration
    cnfg_data = u.load_config()

    if (debug := cnfg_data['Debug']['enable']):
        print("Configuration:")
        u.print_config(cnfg_data)

    # load transcripts and remove version (if exists)
    transcripts = [x.split('.')[0] for x in u.load_transcripts_text(cnfg_data)]
    if debug:
        print(f"\nLoaded {len(transcripts):,} transcripts:\n{transcripts}")

    # get IDs
    transcripts_ids = u.get_transcripts_IDs(cnfg_data, transcripts)
    if debug:
        print("\nIDs:")
        for k, v in transcripts_ids.items():
            print(f"{k}: {v}")

    # get domains
    transcripts_domains = u.get_uniprot_domains(cnfg_data, transcripts_ids)

    # generate output table (only in debug mode)
    if debug:
        dfs, transcript_IDs = u.generate_output_table(cnfg_data, transcripts_domains)
        if len(dfs) == 1:
            print("Aggregate domains table:")
            print(dfs[0].to_string())
        else:
            for df, transcript_ID in zip(dfs, transcript_IDs):
                print(f"Domains table for {transcript_ID}:")
                print(df.to_string())

    # write to output file
    u.generate_output_file(cnfg_data, transcripts_domains)
    if debug:
        print(f"\n{cnfg_data['Output']['file']} generated.")


if __name__ == '__main__':
    main()
