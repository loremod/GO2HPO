import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class DataManager:
    def __init__(self):
        # Attributes
        self.gene_symbol_id = {}  # Gene Symbol-ID associations (dictionary) 
        self.hpo_name_id = {}     # HPO name-ID associations (dictionary) 
        self.go_id_name_association = {}  # GO ID-name associations (dictionary) 
        self.hpo_gene_data = None  # HPO-gene data (binary matrix)
        self.go_gene_data = None   # GO-gene data (binary matrix)
        self.hpo_selected = []     # Selected HPO list
        self.go_selected = []      # Selected GO list



    
    def getCountDataset(self, df, id="hpo_id", to_count="gene_symbol", new_count_name = 'count'):
        df_count = df.groupby(id)[to_count].nunique()
        df_count = df_count.reset_index()
        df_count.columns = [id, new_count_name]
        return df_count

    def filterByValueBoundaries(self, df, L_bound=None, R_bound=None, is_df_counts = False, id="hpo_id", 
                            to_count="gene_symbol", count_name='count'):
        # return the original df in which both are None (misuse of the function)
        if L_bound == None and R_bound == None:
            return df
        
        countDf = None
        if is_df_counts == False:
            countDf = self.getCountDataset(df, id=id, to_count=to_count, new_count_name=count_name)
        else:
            countDf = df

        filterCountDf = countDf
        if L_bound != None:
            filterCountDf = filterCountDf[filterCountDf[count_name] >= L_bound] 
        if R_bound != None:
            filterCountDf = filterCountDf[filterCountDf[count_name] <= R_bound]

        # Keep only the IDs that meet the threshold
        filtered_ids = filterCountDf[id]
        # Filter the dataframe to keep only rows with IDs that meet the threshold
        filteredDf = df[df[id].isin(filtered_ids)]
        return filteredDf

    # generate hpo-gene matrix
    def generate_hpo_gene_matrix(self, file_path, L_bound = 50, R_bound = 100):
        df_genes = pd.read_csv(file_path, sep="\t")
        unique_gene_hpo_pairs = df_genes[['gene_symbol', 'ncbi_gene_id', 'hpo_id', 'hpo_name']].drop_duplicates()
        filteredHpoDf = self.filterByValueBoundaries(unique_gene_hpo_pairs, L_bound, R_bound)
        sparse_pandas = pd.crosstab(filteredHpoDf['gene_symbol'], filteredHpoDf['hpo_id']).astype(pd.SparseDtype("int", fill_value=0))
        self.hpo_gene_data = sparse_pandas

    # Import the HPO2GENES file
    def importHPO2GeneFile(self, file_path):
        df_hpo_genes = pd.read_csv(PHO2GENES_PATH, sep="\t")
        self.generate_hpo_gene_matrix(df_hpo_genes)

    # Getteer for custom dataset
    def get_dataset(self, hpo_list=None, go_list=None):
        """
        Retrieve datasets based on selected HPO and GO terms and combine them into a single matrix.

        :param hpo_list: List of one or more HPO terms (optional).
        :param go_list: List of GO terms (optional).
        :return: Combined dataset as a single matrix.
        """

        # Use hpo_selected if hpo_list is None
        if hpo_list is None:
            hpo_list = self.hpo_selected

        # Use go_selected if go_list is None
        if go_list is None:
            go_list = self.go_selected

        # Filter the HPO-gene data based on the HPO list
        hpo_filtered_data = pd.DataFrame()
        if self.hpo_gene_data is not None and hpo_list:
            hpo_filtered_data = self.hpo_gene_data.loc[:, hpo_list]

        # Filter the GO-gene data based on the GO list
        go_filtered_data = pd.DataFrame()
        if self.go_gene_data is not None and go_list:
            go_filtered_data = self.go_gene_data.loc[:, go_list]

        # Ensure indices represent the same set of genes
        if not hpo_filtered_data.empty and not go_filtered_data.empty:
            common_genes = hpo_filtered_data.index.intersection(go_filtered_data.index)
            hpo_filtered_data = hpo_filtered_data.loc[common_genes]
            go_filtered_data = go_filtered_data.loc[common_genes]

        # Combine the matrices
        combined_data = pd.concat([hpo_filtered_data, go_filtered_data], axis=1)

        return combined_data

    # Function to filter HPO terms
    def filter_hpo(self, method, params=None):
        """
        Filter HPO terms based on the specified method and parameters.
        
        :param method: Filtering method to use.
        :param params: Optional parameters for the filtering method.
        """
        # Implementation for filtering HPO terms
        pass

    # Function to filter GO terms
    def filter_go(self, method, params=None):
        """
        Filter GO terms based on the specified method and discussions with the professor.
        
        :param method: Filtering method to use.
        :param params: Optional parameters for the filtering method.
        """
        # Implementation for filtering GO terms
        pass





if __name__ == "__main__":
    DATA_PATH = "../../preparation/codice/"
    PHO2GENES_PATH = f"{DATA_PATH}phenotype_to_genes.txt"
    df_genes = pd.read_csv(PHO2GENES_PATH, sep="\t")
    print(df_genes.head())
    data_manager = DataManager()