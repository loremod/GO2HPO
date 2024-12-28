from setuptools import setup, find_packages

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="GO2HPO",  
    version="0.1.0",
    packages=find_packages(where="src"),  # Automatically finds packages in the "src" folder
    package_dir={"": "src"},  # Tells Python that the source code is in "src"
    install_requires=install_requires,  # I added the requirements automatically, but maybe manually would be better to add only the needed ones
    description="Creates association rules that classify genes to related phenotypes given their annotations",  # Replace with a description
    author="Lorenzo Modica", 
    author_email="lorenzomod@hotmail.it",  
    url="https://github.com/loremod/GO2HPO", 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
