import pandas as pd
import json
import os
import numpy as np
from scipy.sparse import load_npz
from scipy.sparse import save_npz

class DataExporter:

    def export_attributes(self, _object, export_dir: str):
        """
        Exports all attributes of the object.
        
        Parameters:
        - export_dir: Directory to save the exported files.
        """
        if not os.path.exists(export_dir):
            os.makedirs(export_dir)
        
        for attr_name, attr_value in _object.__dict__.items():
            file_path = os.path.join(export_dir, attr_name)
            
            if isinstance(attr_value, dict):
                # Export dictionary as JSON
                with open(f"{file_path}.json", "w") as json_file:
                    json.dump(attr_value, json_file)
            
            elif isinstance(attr_value, list):
                # Export list as JSON
                with open(f"{file_path}.json", "w") as json_file:
                    json.dump(attr_value, json_file)
            
            elif isinstance(attr_value, pd.DataFrame):
                # Export Pandas DataFrame (sparse)
                if isinstance(attr_value.dtypes.iloc[0], pd.SparseDtype):
                    # Save sparse DataFrame as sparse matrix and index/columns separately
                    sparse_matrix = attr_value.sparse.to_coo()
                    save_npz(f"{file_path}.npz", sparse_matrix)
                    attr_value.index.to_series().to_csv(f"{file_path}_index.csv", index=False)
                    attr_value.columns.to_series().to_csv(f"{file_path}_columns.csv", index=False)
                else:
                    # If not sparse, save as CSV
                    attr_value.to_csv(f"{file_path}.csv", index=True)
            
            else:
                print(f"Unsupported attribute type: {type(attr_value)} for attribute {attr_name}")



class DataImporter:
    def import_attributes(self, _object, import_dir: str):
        """
        Imports all attributes from a directory.
        
        Parameters:
        - import_dir: Directory to load the exported files from.
        """
        for file_name in os.listdir(import_dir):
            print(file_name)

        for file_name in os.listdir(import_dir):
            file_path = os.path.join(import_dir, file_name)
            attr_name, ext = os.path.splitext(file_name)
            
            if ext == ".json":
                with open(file_path, "r") as json_file:
                    setattr(_object, attr_name, json.load(json_file))
            
            elif ext == ".csv":
                # Load as non-sparse DataFrame
                if "_index" not in file_name and "_columns" not in file_name:
                    setattr(_object, attr_name, pd.read_csv(file_path, index_col=0))
                else:
                    continue
            
            elif ext == ".npz":
                # Load sparse matrix
                sparse_matrix = load_npz(file_path)
                index = pd.read_csv(os.path.join(import_dir, f"{attr_name}_index.csv")).squeeze("columns")
                columns = pd.read_csv(os.path.join(import_dir, f"{attr_name}_columns.csv")).squeeze("columns")
                sparse_df = pd.DataFrame.sparse.from_spmatrix(sparse_matrix, index=index, columns=columns)
                setattr(_object, attr_name, sparse_df)
            
            else:
                raise TypeError(f"Unsupported file type: {ext} for file {file_name}")
