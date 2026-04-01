from peptide_druggability_assessment import batch_process_to_csv
import pandas as pd

pep = pd.read_csv('test_peptides.csv')
batch_process_to_csv(pep.sequence, "test_peptides_results.csv")