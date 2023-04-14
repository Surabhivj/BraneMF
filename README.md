# BraneMF: Integration of Biological Networks for Functional Analysis of Proteins

# Citation
**Surabhi Jagtap, Abdulkadir Çelikkanat, Aurélie Pirayre, Frédérique Bidard, Laurent Duval, Fragkiskos D Malliaros, BraneMF: integration of biological networks for functional analysis of proteins, Bioinformatics, Volume 38, Issue 24, 15 December 2022, Pages 5383–5389, https://doi.org/10.1093/bioinformatics/btac691**


![BraneMF](https://user-images.githubusercontent.com/47250394/144092544-0ca33e5a-0d08-4a7a-833b-5edca24a7a61.png)

#### Create a virtual environment with conda

**1.** Create a virtual environment with name 'branemf'
```
conda create -n branemf python=3.6
```

**2.** Activate the virtual environment
```
conda activate branemf
```

**3.** Install the required packages
```
pip install -r requirements.txt
```

**4.** Download the precomputed PPMI matrix files

```
Google drive: https://drive.google.com/drive/folders/1X5Gj5udIPKiLEvzeKWuUucnWw5AZOsnY?usp=sharing 
```
To compute PPMI matrices from String networks, run the following command by providing taxanomy ID of the organism (e.g 4932: yeast, 9606: human)
other organisms: https://stringdb-static.org/download/species.v11.5.txt 

```
 python BraneMF_buildPPMI.py --O 4932

```


**5.** Computation of embeddings

```
run file: branemf.m 
```
**To reproduce the results, use precomuted embeddings and following commands**

**6.** Perform clustering and GO enrichment

```
python cluster_enrichment.py --emb ./data/emb/yeast_branemf_w1_alpha_1.emb --k 40 --sim 20 --genes ./data/yeast_string_genes.txt
```

**7.** Perform protein function prediction
```
python Protein_function_prediction.py --emb ./emb/Integrated_emb/BraneMF_emb_sel/w10/Org4932BraneMF_g1_d128_w10_23May2022.emb  --anno ./data/org_4932_BPl1.txt --trials 10
```
**8.** Perform protein Interaction prediction

  **a.** Prepare the training and test sets.
  ```
  python ppi_pre_preprocess_files.py --new ./data/yeast_string_refnet_2021.txt --old ./data/old_ppis.txt --genes ./data/yeast_string_genes.txt
  ```
  **b.** compute the scores
  ```
  python predict_ppi.py --emb ./data/emb/yeast_branemf_w3_alpha_1.emb --train ./data/train.pkl --test ./data/test.pkl
  ```

**9.** Perform Network reconstruction

```
python net_reconst.py --emb ./data/emb/yeast_branemf_w1_alpha_1.emb --refnet ./data/yeast_string_refnet_2021.txt --genes ./data/yeast_string_genes.txt
```

