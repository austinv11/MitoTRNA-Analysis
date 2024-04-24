import os
import subprocess
from tempfile import NamedTemporaryFile

from Bio.SeqIO import read
import pandas as pd
from tqdm import tqdm


def find_anticodon(taxonomy_string: str, true_trna: str, sequence: str) -> str:
    with NamedTemporaryFile("w+t", delete=False) as ftemp:
        ftemp.write(f">temp\n{sequence}")
        ftemp.close()

        path = ftemp.name

        trna_flag = ["-O", '-g', 'gc_other_mito']  # Search for "Other organellar tRNAs"
        if "Metatheria" in taxonomy_string:  # Marsupials
            trna_flag = ['-M', 'mammal', "-g", "gc_marsu_mito"]  # Use alternative marsupial mitochondrial tRNA coding sequences
        elif "Mammalia" in taxonomy_string:
            trna_flag = ["-M", "mammal"]  # Search for mammalian mitochondrial tRNAs
        elif "Vertebrata" in taxonomy_string:
            trna_flag = ["-M", "vert"]  # Search for vertebrate mitochondrial tRNAs
        elif "Echinodermata" in taxonomy_string:
            trna_flag += ['-g', 'gc_echinoderm_mito']  # Use alternative echinoderm mitochondrial tRNA coding sequences
        elif 'Vertebrata' not in taxonomy_string:  # Check if invertebrate
            trna_flag += ["-g", 'gc_invert_mito']  # Use alternative invertebrate mitochondrial tRNA coding sequences

        # -q = Quiet mode
        # -D = Disable pseudogene detection
        # --brief = No header, just the raw tRNA information
        # --max = "Maximum sensitivity" mode
        # -X 0 = Bit score cutoff, no cutoff since we know that we have a tRNA
        result = subprocess.run(["tRNAscan-SE", "--max", '-q', '-D', '-X', '0', '--brief', *trna_flag, path], stdout=subprocess.PIPE, text=True)
        if result.returncode != 0:
            return None
        outputs = result.stdout.split("\n")
        if len(outputs) < 2:
            return None
        for output in outputs:
            output = output.strip().split("\t")
            # Columns are:
            # 0: Sequence Name
            # 1: Number of tRNAs found
            # 2: tRNA Sequence begin position
            # 3: tRNA Sequence end position
            # 4: tRNA Amino Acid
            # 5: tRNA Anticodon
            # 6: Intron Begin
            # 7: Intron End
            # 8: Infernal score
            # 9: Notes
            if len(output) < 6:
                continue
            amino_acid = output[4].strip()
            if amino_acid not in true_trna:  # Ensure the matched tRNA is the correct one
                continue
            return output[5].strip()
        return None


def main():
    # Extract all tRNAs from mitochondrial genomes
    metadata = pd.read_csv("refseq_mito_genomes/metadata.csv")
    output_file = "all_trnas.csv"
    if os.path.exists(output_file):
        output = pd.read_csv(output_file).to_dict(orient='list')
    else:
        output = dict(
            identifier=[],
            species=[],
            taxonomy=[],
            name=[],
            sequence=[],
            anticodon_sequence=[],
        )
    for i, row in tqdm(metadata.iterrows(), "Extracting tRNAs...", total=len(metadata)):
        identifier = row["identifier"]
        species = row["species"]
        n_trnas = row["n_trnas"]
        path = row["path"]
        taxonomy = row["taxonomy"]
        if n_trnas == 0:
            continue

        # Extract tRNAs by parsing the GenBank file
        record = read(path, "genbank")
        tRNAs = [feature for feature in record.features if feature.type == "tRNA"]
        for tRNA in tRNAs:
            name = tRNA.qualifiers["product"][0]
            anticodon_sequence = None
            try:
                sequence = str(tRNA.location.extract(record.seq))
            except ValueError:
                continue
            if species in output["species"] and name in output["name"] and sequence in output["sequence"]:
                continue
            # Check if annotated anticodon is present so we can skip the tRNAscan-SE step
            if "anticodon" in tRNA.qualifiers:
                anticodon_position = tRNA.qualifiers["anticodon"][0]
                if 'pos:' in anticodon_position:
                    # Try to parse the position, example: '(pos:complement(2798..2800),aa:Ser,seq:tga)'
                    anticodon_position = anticodon_position.split("pos:")[1].split(",")[0].replace("))", ")")
                    # Extract the anticodon sequence
                    if 'complement(' in anticodon_position:
                        anticodon_position = anticodon_position.split("complement(")[1].split(")")[0]
                        anticodon_end, anticodon_start = anticodon_position.split("..")
                        anticodon_start = int(anticodon_start)
                        anticodon_end = int(anticodon_end)
                    else:
                        anticodon_start, anticodon_end = anticodon_position.split("..")
                        anticodon_start = int(anticodon_start)
                        anticodon_end = int(anticodon_end)
                    anticodon_sequence = str(record.seq[anticodon_start:anticodon_end])
                    continue

            # If the anticodon is not annotated, we need to run tRNAscan-SE
            if anticodon_sequence is None:
                anticodon_sequence = find_anticodon(taxonomy, name, sequence)

            if anticodon_sequence is None:
                # Cannot be annotated
                print("WARNING: Could not find anticodon for", name, "in", species)

            output["identifier"].append(identifier)
            output["species"].append(species)
            output["taxonomy"].append(taxonomy)
            output["name"].append(name)
            output["sequence"].append(sequence)
            output["anticodon_sequence"].append(anticodon_sequence or "MISSING")

        if i % 100 == 0:  # Save progress
            pd.DataFrame(output).to_csv(output_file, index=False)

    # Final output
    pd.DataFrame(output).to_csv(output_file, index=False)


if __name__ == '__main__':
    main()
