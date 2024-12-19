import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from scipy.stats import chi2_contingency, fisher_exact

class DataManager:
    def __init__(self):
        # Attributes
        self.gene_id_symbol = {}  # Gene Symbol-ID associations (dictionary) 
        self.hpo_name_id = {}     # HPO name-ID associations (dictionary) 
        self.go_id_name_association = {}  # GO ID-name associations (dictionary) 
        self.hpo_gene_data = None  # HPO-gene data (binary matrix)
        self.go_gene_data = None   # GO-gene data (binary matrix)


    def _associateId2Symbol(self, row, id_column, symbol_column):
        self.gene_id_symbol[row[id_column]] = row[symbol_column]

    def updateGeneIdSymbolAssoc(self, data, id_column='ncbi_gene_id', symbol_column='gene_symbol', reset=False):
        if reset == True:
            self.gene_id_symbol = {}
        data.apply(self._associateId2Symbol, axis=1, args=(id_column,symbol_column))

    def getCountDataset(self, df, id="hpo_id", to_count="ncbi_gene_id", new_count_name = 'count'):
        df_count = df.groupby(id)[to_count].nunique()
        df_count = df_count.reset_index()
        df_count.columns = [id, new_count_name]
        return df_count

    def filterByValueBoundaries(self, df, L_bound=None, R_bound=None, is_df_counts = False, id="hpo_id", 
                            to_count="ncbi_gene_id", count_name='count'):
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
    def generate_hpo_gene_matrix(self, df, L_bound = None, R_bound = None):
        unique_gene_hpo_pairs = df[['ncbi_gene_id', 'gene_symbol', 'hpo_id', 'hpo_name']].drop_duplicates()
        self.updateGeneIdSymbolAssoc(data=unique_gene_hpo_pairs[['ncbi_gene_id', 'gene_symbol']], id_column='ncbi_gene_id', symbol_column='gene_symbol', reset=False)
        filteredHpoDf = self.filterByValueBoundaries(unique_gene_hpo_pairs, L_bound, R_bound)
        sparse_pandas = pd.crosstab(filteredHpoDf['ncbi_gene_id'], filteredHpoDf['hpo_id']).astype(pd.SparseDtype("int", fill_value=0))
        self.hpo_gene_data = sparse_pandas

    # Import the HPO2GENES file
    def importHPO2GeneFile(self, file_path, L_bound = None, R_bound = None):
        df_hpo_genes = pd.read_csv(file_path, sep="\t")
        self.generate_hpo_gene_matrix(df_hpo_genes, L_bound = L_bound, R_bound = R_bound)


    def importGO2GeneFile(self, go_ontology_path, gene2go_path, taxids):
        godag = GODag(go_ontology_path)
        gene2go_reader = Gene2GoReader(gene2go_path, godag=godag)
        human_gene2go = gene2go_reader.get_id2gos(taxids=taxids)

        # Create a row per Gene-Go_term association
        rows = []
        for gene, go_terms in human_gene2go.items():
            for go_term in go_terms:
                rows.append({"Gene": gene, "GO": go_term})

        gene2go_assoc = pd.DataFrame(rows)
        self.go_gene_data = gene2go_assoc.pivot_table(index="Gene", columns="GO", aggfunc=lambda x: 1, fill_value=0)
        self.go_gene_data = self.go_gene_data.astype(pd.SparseDtype("int", fill_value=0))
        self.go_gene_data.columns = self.go_gene_data.columns.get_level_values(0)


    # Getteer for custom dataset
    def get_dataset(self, hpo_list=None, go_list=None):
        """
        Retrieve datasets based on selected HPO and GO terms and combine them into a single matrix.

        :param hpo_list: List of one or more HPO terms (optional).
        :param go_list: List of GO terms (optional).
        :return: Combined dataset as a single matrix.
        """

        # Use all of them if hpo_list is None
        if hpo_list is None:
            hpo_list = self.hpo_gene_data.columns

        # Use all of them if go_list is None
        if go_list is None:
            go_list = self.go_gene_data.columns

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


    def get_id_associations(self):
        return self.gene_id_symbol

    # Show data

    def hpo_head(self, columns = []):
        if columns != []:
            selected_df = self.hpo_gene_data.loc[:, columns]
            return selected_df.head()
        return self.hpo_gene_data.head()
    
    def go_head(self, columns):
        if columns != []:
            selected_df = self.go_gene_data.loc[:, columns]
            return selected_df.head()
        return self.go_gene_data.head()


class LogManager:
    def __init__(self, is_active=True):
        self.is_silent = not is_active
    def activate(self):
        self.is_silent = False
    def deactivate(self):
        self.is_silent = True
    def toggle(self):
        self.is_silent = not self.is_silent
    def log(self, *args):
        if not self.is_silent:
            print(*args)


if __name__ == "__main__":
    # Prepare the logger
    logger = LogManager(is_active=True)

    # Initialize the DataManager class and import all the data
    DATA_PATH = "../../preparation/codice/"
    HPO2GENES_PATH = f"{DATA_PATH}phenotype_to_genes.txt"
    HPO2GENES_PROVA_PATH = f"{DATA_PATH}HPO2Genes_head.csv"
    GO_ONTOLOGY_PATH = f"{DATA_PATH}go-basic.obo" 
    GENE2GO_PATH = f"{DATA_PATH}gene2go"  # Path to the downloaded gene2go file

    data_manager = DataManager()

    # - HPO2Gene File
    logger.log("Importing HPO2Genes file...")
    data_manager.importHPO2GeneFile(HPO2GENES_PATH, L_bound = 50, R_bound = 100)
    logger.log("HPO2Genes file imported.")
    print(data_manager.hpo_head())
    # logger.log("Print the Gene ID - Gene Symbol associations")
    # print(data_manager.get_id_associations())


    # - GO2Gene File
    humanTaxID = 9606
    GO_taxonomies = [humanTaxID]
    logger.log("\nImporting GO2Genes file...")
    data_manager.importGO2GeneFile(go_ontology_path=GO_ONTOLOGY_PATH, gene2go_path=GENE2GO_PATH, taxids = GO_taxonomies)
    logger.log("GO2Genes file imported.")
    print(data_manager.go_head())


    exit()
