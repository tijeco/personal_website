---
title: Chi-squared analysis on multiple sequence alignments with python
# subtitle: Welcome 👋 We know that first impressions are important, so we've populated your new site with some initial content to help you get familiar with everything in no time.

# Summary for listings and search engines
summary: Welcome 👋 We know that first impressions are important, so we've populated your new site with some initial content to help you get familiar with everything in no time.

# Link this post with a project
projects: []

# Date published
date: "2021-02-21T00:00:00Z"

# Date updated
lastmod: "2021-02-21T00:00:00Z"

# Is this an unpublished draft?
draft: false

# Show this page in the Featured widget?
featured: false

# Featured image
# Place an image named `featured.jpg/png` in this page's folder and customize its options here.
image:
  caption: 'Image credit: [**Unsplash**](https://unsplash.com/photos/CpkOjOcXdUY)'
  focal_point: ""
  placement: 2
  preview_only: false

authors:
- admin


tags:
- python
- protein
- bioinformatics


categories:
- tutorial
- python
- bioinformatics
---
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tijeco/personal_website/blob/master/content/post/chi_squared_analysis_on_multiple_sequence_alignment_with_python/chi_squared_analysis_on_multiple_sequence_alignment_with_python.ipynb)


# Overview

For a general chi-square analysis for an alignment it is calulated
$$chi2 = sum[i from 1 to k] (O_i - E_i)^2 / E_i$$

where k is the size of the alphabet (e.g. 4 for DNA, 20 for amino acids) and the values 1 to k correspond uniquely to one of the nucleotides or amino acids.
O_i is the nucleotide or amino acid frequency in the sequence tested.
E_i is the nucleotide or amino acid frequency expected from the ‘master’ distribution (e.g. the overall frequencies - depends on what one is using).

Whether the nucleotide (or amino acid) composition deviates significantly for the ‘master’ distribution is done by testing the chi2 value using the chi2-distribution with k-1 degrees of freedom (df=3 for DNA or df=19 for amino acids).

# Python functions

## Loading multiple sequence alignment

Using the biopython AlignIO function, we can load in a multiple sequence alignment file of a variety of formats.



```python
alignment_file = "path/to/alignment.fasta"
alignment = AlignIO.read(open(alignment_file), "fasta")
```

## Create composition matrix


```python
import pandas as pd
def compositionMatrix(aln):
	compDict = {}
	fixedCharacters = ["-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
	for record in aln:
		header = record.id
		seq = record.seq
		currentSeqMat = [0]*21
		for i in range(len(seq)):
			aa = seq[i]
			try:
				characterPos = fixedCharacters.index(aa)
				currentSeqMat[characterPos]+= 1
			except:
				print("ERROR:", header, "contains character ("+aa+") not in the list:",fixedCharacters)
		compDict[header] = currentSeqMat
	compDF = pd.DataFrame.from_dict(compDict, orient='index',
                       columns=fixedCharacters)
	return compDF

```

## run chi-squared analysis 



```python
from scipy import stats
import numpy as np
import pandas as pd 

def chi2test(compDF):
	seqTotals = compDF.sum(axis=1)
	gaps = compDF["-"]
	gapsPerSeq = gaps/seqTotals

	nonGap = compDF.loc[:, 'A':'Y']
	nonGapTotals = nonGap.sum().to_frame()
	nonGapSeqTotals = nonGap.sum(axis=1).to_frame()
	numCharacters = nonGapTotals.sum()
	expectedFreq = nonGapTotals / numCharacters

	expectedCountArray = np.dot(nonGapSeqTotals,expectedFreq.transpose())
	expectedCountDF = pd.DataFrame(expectedCountArray,columns =nonGap.columns, index =nonGap.index.values )
	chi2DF = ((expectedCountDF - nonGap)**2)/expectedCountDF
	chi2Sum = chi2DF.sum(axis=1)

	pValueDf = 1 - stats.chi2.cdf(chi2Sum, 19)
	outDF = pd.DataFrame({"Gap/Ambiguity":gapsPerSeq,"p-value":pValueDf})
	outDF.index.name='header'
	
	return outDF

```
