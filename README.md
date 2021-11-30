# BraneMF: Integration of Biological Networks for Functional Analysis of Proteins
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

**4.** Download the yeast data and precomputed embeddings

```
Google drive: 
```
**5.** Computation of embeddings

```
run file: branemf.m 
```

**6.** Perform clustering and GO enrichment

```
python cluster_enrichment.py --emb ./data/emb/yeast_branemf_w1_alpha_1.emb --k 40 --sim 20 --genes ./data/yeast_string_genes.txt
```

**7.** Perform protein function prediction

**8.** Perform protein Interaction prediction

**9.** Perform Network reconstruction
