# TCGA-SKCM Mutational Signature Analysis
### So here's what this project is about, I took real melanoma genomic data from TCGA and ran it through a complete mutational signature analysis pipeline. The goal was simple: figure out what biological processes are actually responsible for the mutations piling up in melanoma genomes.
The analysis includes:
1. Downloading somatic mutation data from TCGA/GDC
2. Processing MAF mutation files
3. Generating SBS96 mutational matrices
4. COSMIC signature fitting using SigProfilerAssignment
5. Exposure visualization and similarity analysis
6. Biological interpretation of mutational signatures

Spoiler: it's mostly the sun.

## Why this project was done
Cancer genomes accumulate mutations through multiple biological and environmental processes. Each process leaves a characteristic mutation pattern known as a mutational signature.
Examples:

Mutational Process	            Associated Signature
UV radiation	                          SBS7
Tobacco smoking	                          SBS4
Aging	                                  SBS1
APOBEC activity	                       SBS2/SBS13
Homologous recombination deficiency	  SBS3

## Why Melanoma?
Honestly, melanoma just made sense as a starting point. It has this insanely high mutation burden compared to most cancers, and UV damage leaves such a strong, clean signal in the data — like SBS7 just dominates. So when you actually run the analysis, the results are almost intuitive. You can see the biology in the numbers, which doesn't always happen.
And the TCGA-SKCM dataset being publicly available made it way more practical. I could actually follow through the whole pipeline without hitting walls.
Melanoma was selected because:
1. melanoma has high mutation burden
2. UV signatures are highly prominent
3. Biological interpretation is straightforward
4. Excellent dataset for learning mutational signatures

The goal of this project was to computationally identify these mutational processes from real TCGA melanoma data.

## What Even Are Mutational Signatures?

So here's the thing — DNA doesn't just randomly break. Different things damage it in very specific ways, and those ways leave patterns. Consistent, reproducible patterns.

UV radiation? It causes C→T mutations, almost always at dipyrimidine sites. Tobacco smoke leaves C→A transversions. And APOBEC enzymes — which are actually part of your own immune system, which I find kind of wild — mutate cytosines in very particular sequence contexts. Your body's defense mechanism quietly leaving marks on your genome. Weird, right?

Now here's where it gets interesting. These patterns are consistent enough across thousands of tumors that researchers started cataloguing them. That catalogue is called COSMIC signatures — and there are now 79+ of them. Each one is basically a fingerprint. A biological fingerprint of whatever process caused the damage.

The way we actually capture this is through the SBS96 framework. Every single mutation gets categorized by its trinucleotide context — meaning, what base mutated, what it changed to, and what bases were sitting on either side of it. Do the math: 6 substitution types × 16 possible flanking combinations = 96 categories total. It sounds complicated but once you see it visually it clicks immediately. Honestly it's kind of elegant for something built to track cancer mutations.

So a mutational signature is essentially just a pattern across those 96 categories. And when you see the same pattern showing up repeatedly across different patients? That's a biological process leaving its mark — consistently enough that you can name it, study it, and eventually maybe target it.

## SBS96 Framework:
Single base substitutions (SBS) are categorized using trinucleotide context.

Example:

A[C>T]G

Meaning:

original base = C

mutated base = T

upstream base = A

downstream base = G

There are:
A. 6 substitution classes
B. 4 upstream possibilities
C. 4 downstream possibilities

Total: 6 × 4 × 4 = 96 contexts
This is known as the SBS96 framework.

## Computational Theory

### Mutational Matrix
Rows: mutation contexts (96)
Columns: tumor samples
Values: mutation counts

Example:

Mutation Type	       Sample1

A[C>T]G	              120

T[C>A]A	              45


### COSMIC Signature Fitting
The observed mutation profile is reconstructed using known COSMIC signatures.

Mathematical model: 
                   M ≈ Σ(wᵢSᵢ)

Where:

M = observed mutation profile

Sᵢ = known COSMIC signature

wᵢ = exposure weight

Method used: Non-negative Least Squares (NNLS)

Reason: Mutation counts cannot be negative

## Dataset

Source: TCGA GDC Portal

Project: TCGA-SKCM (Skin Cutaneous Melanoma)

Data type: Masked Somatic Mutations (MAF format)

Strategy: Whole Exome Sequencing

Samples analyzed: 97 melanoma tumors, 18,025 SNVs

## Tools Used

Tool                                      What it does

GDC Client                                Downloads TCGA data

SigProfilerMatrixGenerator                Builds the SBS96 mutation matrix

SigProfilerAssignment                     Fits COSMIC signatures to your data

scikit-learn                              Cosine similarity between samples

matplotlib                                Visualization

## Full Workflow

### 1. Download the GDC Client

bash

wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip

unzip gdc-client_v1.6.1_Ubuntu_x64.zip

chmod +x gdc-client


### 2. Download Mutation Data

Grab your manifest from the GDC portal (filter by TCGA-SKCM → Masked Somatic Mutation → MAF), then:

bash

./gdc-client download -m gdc_manifest.txt


### 3. Organize the MAF Files

bash

mkdir maf_files

find . -name "*.maf.gz" -exec cp {} maf_files/ \;

gunzip maf_files/*.maf.gz


### 4. Set Up the Environment

bash

conda create -n mutsig python=3.10 -y

conda activate mutsig

pip install SigProfilerMatrixGenerator SigProfilerAssignment SigProfilerExtractor scikit-learn


### 5. Install the Reference Genome

python

from SigProfilerMatrixGenerator import install as genInstall

genInstall.install('GRCh38')

This step will take time to download Reference genome while you can have Coffee or drink.


### 6. Generate the SBS96 Matrix

python

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

matGen.SigProfilerMatrixGeneratorFunc(

    "SKCM",
    
    "GRCh38",
    
    "/path/to/maf_files",
    
    plot=True
)

Output : This generates SBS96, DBS, and indel matrices. The one we care about most is SKCM.SBS96.all.


### 7. COSMIC Signature Fitting

python

import pandas as pd

from SigProfilerAssignment import Analyzer as Analyze

matrix = pd.read_csv(

    "/path/to/output/SBS/SKCM.SBS96.all",
    
    sep="\t", index_col=0
)

Analyze.cosmic_fit(
    
    samples=matrix,
    
    output="/path/to/SKCMAssignment",
    
    input_type="matrix",
    
    context_type="96",
    
    genome_build="GRCh38",
    
    cosmic_version=3.4
)

Under the hood this uses Non-negative Least Squares (NNLS) — because mutation counts can't go negative. It decomposes your observed profile as a weighted sum of known COSMIC signatures.


8. Exposure Heatmap

python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "/path/to/SKCMAssignment/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    sep="\t", index_col=0
)

plt.figure(figsize=(15, 8))
plt.imshow(df.T, aspect='auto')
plt.colorbar(label="Exposure")
plt.xlabel("Samples")
plt.ylabel("Signatures")
plt.title("Mutational Signature Exposure Heatmap")
plt.tight_layout()
plt.show()


9. Cosine Similarity Between Samples

python
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd

df = pd.read_csv(
    "/path/to/output/SBS/SKCM.SBS96.all",
    sep="\t", index_col=0
)

sim = cosine_similarity(df.T)
print(sim)

Values above 0.85 mean the samples basically share the same mutational history. Below 0.3 and they're quite different — probably driven by different secondary processes on top of UV.
