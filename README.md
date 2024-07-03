# PanRNA [![DOI](https://zenodo.org/badge/181719067.svg)](https://zenodo.org/doi/10.5281/zenodo.12626836)
PanRNA is a python written pipeline to map and count reads of an RNA-seq experiment. We divide the pipeline in two parts. The first one, we call Basic Run, is the the most common task on RNA-seq which only map the reads to a reference genome and return counts files for each sample. These can be used for differential expression analysis. The second part, the full PanRNA run, includes an extra task meant to avoid loosing information due to a lack in the reference genome annotation or due to the presence of genes in the samples that are not in the reference genome used.

For more details see the [manual](https://github.com/Mimop/panrna/blob/master/PanRNA_manual.pdf).


![PanRAN](https://github.com/Mimop/panrna/blob/master/panrna.jpg)

# Releases
panrna_v2.1.0 is the original pipeline that was tested with python 2.7
panrna_v2.1.11 should work correctly with python=3.x versions in most systems

# Citing 
If PanRNA is useful for your work please use :

_PanRNA: An RNAseq pipeline for semi-reference guided analysis_ Miguel Morard & Javier Alonso del Real.

https://doi.org/10.5281/zenodo.12626837
https://github.com/Mimop/panrna

*Help* :miguel.morard@iata.csic.es
