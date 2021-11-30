# BraneMF: Integration of Biological Networks for Functional Analysis of Proteins
![workflow-1](https://user-images.githubusercontent.com/47250394/144040612-bda99618-1a26-4f69-b44c-bb0b167d1f8f.png)

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

**4.** Download the yeast data and precomputed embeddings

```
Google drive: 
```
**5.** Computation of embeddings

**6.** Perform clustering and GO enrichment

```
python cluster_enrichment.py --emb ./data/emb/yeast_branemf_w1_alpha_1.emb --k 40 --sim 20 --genes ./data/yeast_string_genes.txt
```

**7.** Perform protein function prediction
