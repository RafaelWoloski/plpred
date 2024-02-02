import argparse
from plpred.preprocessing import compute_aa_composition
from Bio import SeqIO
import pandas as pd
import pickle



def run_model(file_path:str, model_path:str) -> pd.DataFrame:
    """
    Run a membrane prediction on a FASTA file.

    Parameters
    ----------
    file_path:str
        Proteins in FASTA format.
    model_path:str
        Trained model in pickle format

    Returns
    -------
    df_predictions:pd.DataFrame
        Pandas DataFrame containing the membrane proteins predictions
    """
    df_aa_composition = pd.DataFrame()

    with open(model_path, 'rb') as handle:
        model = pickle.load(handle)
    handle = open(file_path)
    parser = SeqIO.parse(handle, 'fasta')
    
    df_predictions = pd.DataFrame(columns=['id', 'membrane'])
    
    for record in parser:
        aa_composition = compute_aa_composition(str(record.seq))
        aa_composition['id'] = record.id
        df_aa_composition = df_aa_composition.append(aa_composition, ignore_index=True)

    x = df_aa_composition.drop(['id'], axis=1)
    ids = df_aa_composition['id']
    y_pred = model.predict(x)

    df_predictions['id'] = ids
    df_predictions['membrane'] = y_pred

    return df_predictions



def main():
    argument_parser = argparse.ArgumentParser(description='plpred-predict: Membrane protein prediction tool')
    argument_parser.add_argument('-i', '--input', required=True, help='Input file (.fasta)')
    argument_parser.add_argument('-o', '--output',required=True, help='Output file (.csv)')
    argument_parser.add_argument('-m', '--model', required=True, help='Trained model (.pickle)')
    arguments = argument_parser.parse_args()

    df_predictions = run_model(file_path=arguments.input, model_path=arguments.model)
    df_predictions.to_csv(arguments.output, index=False)

if __name__ == '__main__':
    main()

