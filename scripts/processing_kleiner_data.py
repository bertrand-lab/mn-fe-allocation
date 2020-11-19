from pyteomics import pepxml
import pandas as pd
from matplotlib import pyplot

 # read in proteomics data
pep_data_long = pepxml.read('../data/kleiner_data/Run1and2_U1.pep.xml')
pep_data_short = pepxml.read('../data/kleiner_data/Run4and5_U1.pep.xml')

pep_df_long = pepxml.DataFrame(pep_data_long)
pep_df_short = pepxml.DataFrame(pep_data_short)

def get_organism_id(list_of_proteins):
    '''
    function to parse the protein names and collect the unique names of organisms
    '''
    organism_names = []
    # for each protein name, collect the organism
    for protein_name in list_of_proteins:
        # parse the protein name/org name
        org_name_i = protein_name.split('_')[0]
        # add this to the organism names list
        organism_names.append(org_name_i)
    
    if len(set(organism_names)) == 1:
        return(organism_names[0])
    else:
        return('ambiguous')

def get_all_org_ids(peptide_df):
    '''
    go through the pepxml dataframe and collect the organism ids
    '''
    all_peptide_to_org_mappings = []

    for obs_peptide in peptide_df['protein']:
    
        # if there is only one protein hit
        if len(obs_peptide) == 1:
            # parse the protein name to get the organism
            org_i = obs_peptide[0].split('_')[0]
            all_peptide_to_org_mappings.append(org_i)
        else:
            org_i = get_organism_id(obs_peptide)
            all_peptide_to_org_mappings.append(org_i)
            
    return(all_peptide_to_org_mappings)


pep_df_long['org_mappings'] = get_all_org_ids(pep_df_long)
pep_df_short['org_mappings'] = get_all_org_ids(pep_df_short)

# df['DataFrame Column'] = df['DataFrame Column'].astype(float)
pep_df_long['precursor_intensity'] = pep_df_long['precursor_intensity'].astype(float)
pep_df_short['precursor_intensity'] = pep_df_short['precursor_intensity'].astype(float)

# pep_df.dtypes

def format_output(pep_df_out):

    # df.groupby(['A', 'B'], as_index=False)['C'].sum()
    aggregated_intensity_mean = pep_df_out.groupby(['org_mappings'], as_index = False)['precursor_intensity'].mean()
    aggregated_intensity_median = pep_df_out.groupby(['org_mappings'], as_index = False)['precursor_intensity'].median()
    aggregated_intensity_sum = pep_df_out.groupby(['org_mappings'], as_index = False)['precursor_intensity'].sum()
    aggregated_intensity_count = pep_df_out.groupby(['org_mappings'], as_index = False)['precursor_intensity'].count()

    aggregated_intensity_sum.columns = ['org_mappings', 'sum_precursor_intensity']
    aggregated_intensity_median.columns = ['org_mappings', 'median_precursor_intensity']
    aggregated_intensity_mean.columns = ['org_mappings', 'mean_precursor_intensity']
    aggregated_intensity_count.columns = ['org_mappings', 'number_peps']


# read in data of total protein per group

    actual_protein_amount = pd.read_csv('../data/kleiner_data/uneven_protein_amount_Kleiner2017.csv')
    
    combined_df = aggregated_intensity_mean.merge(actual_protein_amount)
    combined_df_2 = combined_df.merge(aggregated_intensity_sum)
    combined_df_3 = combined_df_2.merge(aggregated_intensity_median)
    combined_df_4 = combined_df_3.merge(aggregated_intensity_count)

    return(combined_df_4)


formatted_long = format_output(pep_df_long)
formatted_short = format_output(pep_df_short)
formatted_long.to_csv('../data/kleiner_data/kleiner_data_processed_long.csv')
formatted_short.to_csv('../data/kleiner_data/kleiner_data_processed_short.csv')
