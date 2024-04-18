import Bio.Entrez


def main(identifier, email):
    """
    Download the full GenBank record for a given identifier.
    :param identifier: The GenBank identifier to download.
    :param email: The email address to use for the Entrez request.
    Text is printed to standard output.
    """
    Bio.Entrez.email = email
    handle = Bio.Entrez.efetch(db="nucleotide",
                               id=identifier,
                               rettype="gb",
                               retmode="text")
    print(handle.read())


if __name__ == "__main__":
    # Example usage: python download_genome.py NC_012920.1 varelaa@mskcc.org > human_mito.gb
    import sys
    try:
        main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "")
    except IndexError:
        print("Usage: download_genome.py <identifier> <optional: email>")
        sys.exit(1)
