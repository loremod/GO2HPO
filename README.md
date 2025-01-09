# GO2HPO

GO2HPO is a Python package designed to create association rules that classify genes to related phenotypes based on their annotations.

---

## Features
- Extract and process annotation data related to genes and phenotypes.
- Analyze data to create association rules. (TO DO)
- Allows to be easily configurable and extendable. (TO DO)

---

## Installation

To install GO2HPO, follow these steps:

### **1. Clone the Repository**
Clone this repository to your local machine:
```bash
git clone https://github.com/loremod/GO2HPO.git
cd GO2HPO
```

### **2. Set Up a Virtual Environment**
Create and activate a virtual environment:
```bash
# On Windows
python -m venv venv
venv\Scripts\activate

# On Linux/Mac
python3 -m venv venv
source venv/bin/activate
```

### **3. Install the Package**
Install the package in editable mode (TEMPORARILY):
```bash
pip install -e .
```

---

## Usage

### **Importing the Modules**
After installation, you can import and use the package in your Python scripts or interactive sessions:

```python
from GO2HPO import DataManager, LogManager

# Example usage
data_manager = DataManager()
#....
#....
dataset = data_manager.get_dataset(hpo_list=["HP:0004322"], go_list=["GO:0008150"])
print(dataset)
```

---

## Repository Structure
```
GO2HPO/
├── src/
│   ├── GO2HPO/
│   │   ├── __init__.py
│   │   ├── dataManager.py
│   │   ├── LogManager.py
│   │   ├── dataImportExport.py
│   │   ├── StatisticalAnalyzer.py
├── setup.py
├── requirements.txt
├── README.md
├── LICENSE
```

---

## Dependencies

- `numpy`
- `pandas`
- `joblib`
- Others, like the ones for the generitic algorithm (TO DO)

---


**Lorenzo Modica**

Email: [lorenzomod@hotmail.it](mailto:lorenzomod@hotmail.it)

