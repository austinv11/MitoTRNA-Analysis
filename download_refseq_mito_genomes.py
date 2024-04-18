import os
import time
from io import StringIO

import Bio.Entrez
import Bio.SeqIO
import pandas as pd
from tqdm import tqdm


def main(email):
    """
    Find all complete mitochondrial genomes in RefSeq, classified by taxonomy.
    :param email: The email to make the queries under.
    """
    Bio.Entrez.email = email

    if not os.path.exists("refseq_mito_genomes"):
        os.mkdir("refseq_mito_genomes")

    complete = False
    retstart = 0
    count = 0
    all_ids = set()
    while not complete:
        # First, we must query for all complete mitochondrial genomes and exclude plant genomes.
        handle = Bio.Entrez.esearch(db="nucleotide",
                                    term="mitochondrion[Title] AND genome[Title] AND complete[Title] AND RefSeq[Keyword] Metazoa[Organism]",
                                    retmax=10_000, retstart=retstart)

        # Now, we must fetch the records.
        record = Bio.Entrez.read(handle)
        handle.close()
        # If 10,000 records were returned, we need to make another query.
        count += len(record["IdList"])
        complete = count >= int(record["Count"])
        retstart += count
        all_ids.update(record["IdList"])

    print(f"Found {len(all_ids)} complete mitochondrial genomes.")

    mito_info = dict(
        identifier=[],
        species=[],
        n_trnas=[],
        path=[],
        taxonomy=[],
    )
    base_path = "refseq_mito_genomes"
    if not os.path.exists(os.path.join(base_path, 'metadata.csv')):
        completed = []
    else:
        metadata = pd.read_csv(os.path.join(base_path, 'metadata.csv'))
        completed = metadata['identifier'].tolist()
        mito_info = metadata.to_dict(orient='list')

    print(f"Found {len(completed)} completed genomes.")

    for i, identifier in tqdm(list(enumerate(all_ids)), "Fetching Genomes..."):
        if int(identifier) in completed:
            continue

        # We will first fetch the taxonomy information to classify the genome.
        completed = False
        while not completed:
            try:
                handle = Bio.Entrez.efetch(db="nucleotide",
                                           id=identifier,
                                           rettype="gb",
                                           retmode="text")
                record_text = handle.read()
                handle.close()
                completed = True
            except Exception as e:
                print(f"Error fetching {identifier}: {e}")
                print("Sleeping for 5 seconds...")
                time.sleep(5)
        record_info = next(Bio.SeqIO.parse(StringIO(record_text), "genbank"))
        # We will extract the taxonomy information from the record.
        taxonomy = record_info.annotations.get("taxonomy", ["Unknown"])
        species = record_info.annotations.get("organism", taxonomy[-1])

        # Now create the path to save the genome.
        save_path = base_path
        for taxon in taxonomy[:-1]:
            save_path = os.path.join(save_path, taxon)
            if not os.path.exists(save_path):
                os.mkdir(save_path)
        save_path = os.path.join(save_path, f"{species}.gb")
        mito_info["identifier"].append(identifier)
        mito_info["species"].append(species)
        mito_info["path"].append(save_path)
        mito_info["taxonomy"].append(";".join(taxonomy))
        mito_info["n_trnas"].append(len([f for f in record_info.features if f.type == "tRNA"]))

        # Now we will save the genome.
        with open(save_path, "w") as f:
            f.write(record_text)

        completed.append(identifier)

        if i % 100 == 0:
            pd.DataFrame(mito_info).to_csv(os.path.join(base_path, "metadata.csv"), index=False)

    # Now we will save the metadata.
    pd.DataFrame(mito_info).to_csv(os.path.join(base_path, "metadata.csv"), index=False)


if __name__ == "__main__":
    import sys

    main(len(sys.argv) > 1 and sys.argv[1] or "")
