[![](https://img.shields.io/badge/Google%20Colab-Open-orange?logo=google-colab)](https://colab.research.google.com/github/bioexcel/biobb_wf_godmd/blob/main/biobb_wf_godmd/notebooks/biobb_wf_godmd.ipynb)

# Protein conformational transitions calculations tutorial using BioExcel Building Blocks (biobb) and GOdMD

***

This tutorial aims to illustrate the process of computing a **conformational transition** between two known **structural conformations** of a protein, step by step, using the **BioExcel Building Blocks library (biobb)**. 

Examples shown are the calculation of the conformational transition for the **Adenylate Kinase** protein, from the **closed state** (PDB Code [1AKE](https://www.rcsb.org/structure/1AKE)) to the **open state** (PDB Code [4AKE](https://www.rcsb.org/structure/4AKE)). **Adenylate Kinases** are **phosphotransferase enzymes** that catalyze the interconversion of the various **adenosine phosphates** (ATP, ADP, and AMP), and are known to undergo large **conformational changes** during their **catalytic cycle**.

The code wrapped is the ***GOdMD*** method, developed in the **[Molecular Modeling and Bioinformatics](https://mmb.irbbarcelona.org/www/) group** (IRB Barcelona). **GOdMD** determines pathways for **conformational transitions** in macromolecules using **discrete molecular dynamics** and **biasing techniques** based on a combination of **essential dynamics** and **Maxwell-Demon sampling techniques**. A web implementation of the method can be found here: [https://mmb.irbbarcelona.org/GOdMD/index.php](https://mmb.irbbarcelona.org/GOdMD/index.php)

**Exploration of conformational transition pathways from coarse-grained simulations.**<br>
*Sfriso P, Hospital A, Emperador A, Orozco M.*<br>
*Bioinformatics, 129(16):1980-6.*<br>
*Available at: [https://doi.org/10.1093/bioinformatics/btt324](https://doi.org/10.1093/bioinformatics/btt324)*

***

## Settings

### Biobb modules used

 * [biobb_godmd](https://github.com/bioexcel/biobb_godmd): Tools to compute protein conformational transitions using GOdMD.
 * [biobb_io](https://github.com/bioexcel/biobb_io): Tools to fetch biomolecular data from public databases.
 * [biobb_structure_utils](https://github.com/bioexcel/biobb_structure_utils): Tools to modify or extract information from a PDB structure.
 * [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.


### Auxiliary libraries used

* [emboss](https://www.ebi.ac.uk/Tools/emboss/): Software that automatically copes with data in a variety of formats and even allows transparent retrieval of sequence data from the web.
* [jupyter](https://jupyter.org/): Free software, open standards, and web services for interactive computing across all programming languages.
* [plotly](https://plot.ly/python/offline/): Python interactive graphing library integrated in Jupyter notebooks.
* [nglview](https://nglviewer.org/#nglview): Jupyter/IPython widget to interactively view molecular structures and trajectories in notebooks.
* [simpletraj](https://github.com/arose/simpletraj): Lightweight coordinate-only trajectory reader based on code from GROMACS, MDAnalysis and VMD.

### Conda Installation and Launch

```console
git clone https://github.com/bioexcel/biobb_wf_godmd.git
cd biobb_wf_godmd
conda env create -f conda_env/environment.yml
conda activate biobb_wf_godmd
jupyter-notebook biobb_wf_godmd/notebooks/biobb_wf_godmd.ipynb
```

***

## Tutorial

Click here to [open tutorial in Google Colab](https://colab.research.google.com/github/bioexcel/biobb_wf_godmd/blob/main/biobb_wf_godmd/notebooks/biobb_wf_godmd.ipynb)

***

## Version
2025.1 Release

## Copyright & Licensing
This software has been developed in the [MMB group](http://mmb.irbbarcelona.org) at the [BSC](http://www.bsc.es/) & [IRB](https://www.irbbarcelona.org/) for the [European BioExcel](http://bioexcel.eu/), funded by the European Commission (EU Horizon Europe [101093290](https://cordis.europa.eu/project/id/101093290), EU H2020 [823830](http://cordis.europa.eu/projects/823830), EU H2020 [675728](http://cordis.europa.eu/projects/675728)).

* (c) 2015-2026 [Barcelona Supercomputing Center](https://www.bsc.es/)
* (c) 2015-2026 [Institute for Research in Biomedicine](https://www.irbbarcelona.org/)

Licensed under the
[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0), see the file LICENSE for details.

![](https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png "Bioexcel")
