# clibas
### Combinatorial Library Analysis Suite


To get it going:

1. Download the scripts
2. Locate the environment file in clibas/environment 
3. Use it to setup a fresh conda environment. From conda prompt, run:
	conda env create -f clibas_env.yml

4. Activate the environment by running
	conda activate clibas_env
	
5. Build your data analysis pipelines. For a minimal example of executing umap + hdbscan analysis, see clibas/notebooks/umap_hdbscan_demo.html
6. For now, analysis pipeline (.py scripts) should be saved and executed in clibas/src directory (where the config file also lives)
