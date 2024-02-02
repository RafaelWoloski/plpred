from Bio import SeqIO
from Bio.SeqUtils import ProtParam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def compute_aa_composition(sequencia:str) -> dict:
    """
    Computes the aminoacid composition of a given protein sequence

    Parameters
    ----------
    protein_sequence: str
        Sequence of the protein to be processed

    Returns
    -------
    aa_composition: dict
        Dictionary containing the relative abundance of each aminoacid
    """
    protein_analyzer = ProtParam.ProteinAnalysis(str(sequencia))
    aa_composition = protein_analyzer.get_amino_acids_percent()
    return aa_composition

def generate_aa_composition_df(file_path:str, membrane_label:int) -> pd.DataFrame:
    """
    Read a .fasta file and return the aminoacid composition of each protein, as well as assign it a membrane label as a Pandas Dataframe
    The label is defined as 1 if the protein is a membrane protein, and 0 if it's a cytoplasmatic protein

    Parameters
    ----------
    file_path: str
        The path for the .fasta file

    membrane_label: int
        The membrane label, defining if it's 0 for cytoplasm or 1 for membrane

    Returns
    -------
    df: pd.DataFrame
        Returns a Pandas Dataframe with the Aminoacid Composition of each protein as well as it's label
    """
    df = pd.DataFrame()
    handle = open(file_path)
    parser = SeqIO.parse(handle,'fasta')
    for proteina in parser:
        protein_data = compute_aa_composition(proteina)
        protein_data['membrane'] = membrane_label
        df = df.append(protein_data, ignore_index=True)
    return df

def main():
    argument_parser = argparse.ArgumentParser(description='plpred-preprocessing: pre-process data to use in training')
    argument_parser.add_argument('-m', '--membrane_proteins', required=True, help="Path to the file containing membrane proteins (.fasta)")
    argument_parser.add_argument('-c', '--cytoplasm_proteins', required=True, help="Path to the file containing cytoplasm proteins (.fasta)")
    argument_parser.add_argument('-o', '--output', required=True, help="Path to the output file (.csv)" )
    arguments = argument_parser.parse_args()

    df_membrane = generate_aa_composition_df(file_path=arguments.membrane_proteins, membrane_label=1)
    df_cytoplasm = generate_aa_composition_df(file_path=arguments.cytoplasm_proteins, membrane_label=0)
    df_processed = pd.concat([df_membrane, df_cytoplasm])
    df_processed.to_csv(arguments.output, index=False)

if __name__ == "__main__":
    main()