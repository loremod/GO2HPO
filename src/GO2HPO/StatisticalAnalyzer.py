from joblib import Parallel, delayed
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np
from typing import Union

class StatisticalAnalyzer:
    def __init__(self, hpo_gene_data, go_gene_data):
        self.hpo_gene_data = hpo_gene_data
        self.go_gene_data = go_gene_data

    def compute_p_value(self, go_column, hpo_column, method = "chi2"):
        #get the series from the column names
        hpo_series = self.hpo_gene_data[hpo_column]
        go_series = self.go_gene_data[go_column]

        # Align both Series to have the same index
        hpo_series, go_series = hpo_series.align(go_series, join='inner')

        # Confusion matrix components
        both = ((go_series == 1) & (hpo_series == 1)).sum()
        only_go = ((go_series == 1) & (hpo_series == 0)).sum()
        only_hpo = ((go_series == 0) & (hpo_series == 1)).sum()
        neither = ((go_series == 0) & (hpo_series == 0)).sum()
        # Contingency table with Laplace Smoothing to prevent some cells to have frequency 0
        contingency_table = [[both + 1, only_hpo + 1],
                            [only_go + 1, neither + 1]]
        p_value = 1
        # Chi-Square test (or switch to Fisher's Exact Test if preferred)
        if method.lower() == "chi2":
            _, p_value, _, _ = chi2_contingency(contingency_table)
        elif method.lower() == "fisher":
            _, p_value = fisher_exact(contingency_table) 
        return p_value
    

    def compute_p_values_batch(self, aligned_data, go_columns, hpo_column, method = "chi2"):
        # aligned_data = self.get_dataset(go_list=go_columns, hpo_list=hpo_column)

        hpo_series = aligned_data[hpo_column]
        # go_series = aligned_data[go_column]

        hpo_present = (hpo_series == 1).values
        hpo_absent = (hpo_series == 0).values

        p_values = []

        for go_col in go_columns:
            go_series = aligned_data[go_col]
            go_present = (go_series == 1).values
            go_absent = (go_series == 0).values

            # Confusion matrix components
            both = (go_present & hpo_present).sum()
            only_go = (go_present & hpo_absent).sum()
            only_hpo = (go_absent & hpo_present).sum()
            neither = (go_absent & hpo_absent).sum()

            # Contingency table
            contingency_table = [[both + 1, only_hpo + 1],
                                [only_go + 1, neither + 1]]

            # Compute p-value
            p_value = 1
            if method.lower() == "chi2":
                _, p_value, _, _ = chi2_contingency(contingency_table)
            elif method.lower() == "fisher":
                _, p_value = fisher_exact(contingency_table)
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            p_values.append(p_value)

        return p_values
    
    def compute_confusion_matrices_vectorized(self, go_matrix, hpo_series):
        """
        Compute confusion matrix components for all GO terms vectorized.

        :param go_matrix: DataFrame or 2D NumPy array, binary GO-term presence/absence.
        :param hpo_series: 1D Series or array, binary HPO-term presence/absence.
        :return: Confusion matrix components as arrays: both, only_go, only_hpo, neither.
        """
        # Convert to NumPy arrays for efficient computation
        go_array = go_matrix.to_numpy()
        hpo_array = hpo_series.to_numpy()

        # Ensure hpo_array aligns with go_array's dimensions
        hpo_present = (hpo_array == 1).reshape(-1, 1)  # Convert to column vector
        hpo_absent = (hpo_array == 0).reshape(-1, 1)   # Convert to column vector

        # Compute confusion matrix components
        both = (go_array & hpo_present).sum(axis=0)  # Both GO and HPO present
        only_go = (go_array & hpo_absent).sum(axis=0)  # GO present, HPO absent
        only_hpo = (~go_array & hpo_present).sum(axis=0)  # HPO present, GO absent
        neither = (~go_array & hpo_absent).sum(axis=0)  # Neither GO nor HPO present

        return both, only_go, only_hpo, neither




    def compute_p_values_vectorized(self, aligned_data, go_columns, hpo_column, method="chi2"):
        """
        Compute p-values for multiple GO columns in a vectorized manner.

        :param aligned_data: Pre-aligned DataFrame containing GO and HPO data.
        :param go_columns: List of GO columns.
        :param hpo_column: The HPO column.
        :param method: Statistical test method ("chi2" or "fisher").
        :return: List of p-values.
        """
        go_matrix = aligned_data[go_columns]
        hpo_series = aligned_data[hpo_column]

        # Compute confusion matrix components using vectorized operations
        both, only_go, only_hpo, neither = self.compute_confusion_matrices_vectorized(go_matrix, hpo_series)

        # Contingency table (with Laplace smoothing)
        contingency_table = np.array([both + 1, only_hpo + 1, only_go + 1, neither + 1]).T.reshape(-1, 2, 2)

        # Compute p-values
        p_values = []
        if method.lower() == "chi2":
            p_values = [chi2_contingency(ct)[1] for ct in contingency_table]
        elif method.lower() == "fisher":
            p_values = [fisher_exact(ct)[1] for ct in contingency_table]
        else:
            raise ValueError(f"Unsupported method: {method}")

        return p_values



    def compute_p_values_parallel(self, aligned_data, go_columns, hpo_column, method="chi2"):
        # Align the dataset once
        # aligned_data = self.get_dataset(go_list=go_columns, hpo_list=hpo_column)

        hpo_series = aligned_data[hpo_column]
        hpo_present = (hpo_series == 1).values
        hpo_absent = (hpo_series == 0).values

        # Helper function for parallel execution
        def compute_p_value_for_column(go_col):
            # go_series = aligned_data[go_col]
            go_present = (aligned_data[go_col] == 1).values
            go_absent = (aligned_data[go_col] == 0).values

            # Confusion matrix components
            both = (go_present & hpo_present).sum()
            only_go = (go_present & hpo_absent).sum()
            only_hpo = (go_absent & hpo_present).sum()
            neither = (go_absent & hpo_absent).sum()

            # Contingency table
            contingency_table = [[both + 1, only_hpo + 1],
                                [only_go + 1, neither + 1]]

            # Compute p-value
            if method.lower() == "chi2":
                _, p_value, _, _ = chi2_contingency(contingency_table)
            elif method.lower() == "fisher":
                _, p_value = fisher_exact(contingency_table)
            else:
                raise ValueError(f"Unsupported method: {method}")

            return p_value
        # Use joblib for parallel computation of p-values
        p_values = Parallel(n_jobs=-1, backend="loky")(
            delayed(compute_p_value_for_column)(go_col) for go_col in go_columns
        )

        return p_values



    # returns a dataframe obtained computing the go-term-wise significance with a specific hpo-term
    # returns the p-values of the statistical tests for each go-term
    # it is possible to apply different statistical tests
    # it is possible to apply a correction
    # it can return the significant go terms only, based on the p_value (if no correction is applied) 
    # or on the adjusted_p_value, if correction is applied
    # By default, the correction parameter is set to None, so no correction is applied
    def compute_significance(self, aligned_data, hpo_column:str, go_columns:list = None,
                                method:str = "chi2", only_significant:bool=True,
                                correction: Union[str, list] = None,
                                approach="batch"):
        # if go_columns is none, all of them are considered (default)
        if go_columns is None:
            go_columns = self.go_gene_data.columns

        p_values = []
        if approach == "batch":
            p_values = self.compute_p_values_batch(aligned_data, go_columns, hpo_column, method=method)
        elif approach == "vectorized":
            p_values = self.compute_p_values_vectorized(aligned_data, go_columns, hpo_column, method=method)
        elif approach == "parallel":
            p_values = self.compute_p_values_parallel(aligned_data, go_columns, hpo_column, method=method)
        elif approach == "naive":
            p_values = [self.compute_p_value(go_col, hpo_column, method=method) for go_col in go_columns]
        elif approach == "restriction": # smart or restriction
            p_values = [] # TO DO
        else:
            p_values = [self.compute_p_value(go_col, hpo_column, method=method) for go_col in go_columns]

        results_df = pd.DataFrame({'GO_Term': go_columns, 'P_Value': p_values})

        if correction == None:
            results_df['Significant']  = results_df['P_Value'] < 0.05
        else:
            if isinstance(correction, str):
                correction = [correction]  # Convert string to list containing the string

            for method in correction:  # Iterate over each correction method in the list
                corrected_results = multipletests(results_df['P_Value'], method=method)
                results_df[f'Adjusted_P_Value_{method}'] = corrected_results[1]  # Corrected p-values
                results_df[f'Significant_{method}'] = corrected_results[0]       # True/False for significance

    
        # Filter significant GO terms
        if only_significant == True:
                if isinstance(correction, str):  # Single correction method
                    return results_df[results_df['Significant']].reset_index(drop=True)
                else:  # Multiple correction methods
                    # Filter rows where at least one "Significant_{method}" column is True
                    significant_columns = [f'Significant_{method}' for method in correction]
                    return results_df[results_df[significant_columns].any(axis=1)].reset_index(drop=True)
        else:
            return results_df


    def compute_significance_optimized(self, go_data, hpo_series,
                                method:str = "chi2", only_significant:bool=True,
                                correction: Union[str, list] = None):

        go_columns = go_data.columns

        p_values = []

        # Compute confusion matrix components using vectorized operations
        both, only_go, only_hpo, neither = self.compute_confusion_matrices_vectorized(go_data, hpo_series)
        # Contingency table (with Laplace smoothing)
        contingency_table = np.array([both + 1, only_hpo + 1, only_go + 1, neither + 1]).T.reshape(-1, 2, 2)

        # Compute p-values
        if method.lower() == "chi2":
            p_values = [chi2_contingency(ct)[1] for ct in contingency_table]
        elif method.lower() == "fisher":
            p_values = [fisher_exact(ct)[1] for ct in contingency_table]
        else:
            raise ValueError(f"Unsupported method: {method}")

        results_df = pd.DataFrame({'GO_Term': go_columns, 'P_Value': p_values})

        if correction == None:
            results_df['Significant']  = results_df['P_Value'] < 0.05
        else:
            if isinstance(correction, str):
                correction = [correction]  # Convert string to list containing the string
            for method in correction:  # Iterate over each correction method in the list
                corrected_results = multipletests(results_df['P_Value'], method=method)
                results_df[f'Adjusted_P_Value_{method}'] = corrected_results[1]  # Corrected p-values
                results_df[f'Significant_{method}'] = corrected_results[0]       # True/False for significance

    
        # Filter significant GO terms
        if only_significant == True:
                if isinstance(correction, str):  # Single correction method
                    return results_df[results_df['Significant']].reset_index(drop=True)
                else:  # Multiple correction methods
                    # Filter rows where at least one "Significant_{method}" column is True
                    significant_columns = [f'Significant_{method}' for method in correction]
                    return results_df[results_df[significant_columns].any(axis=1)].reset_index(drop=True)
        else:
            return results_df