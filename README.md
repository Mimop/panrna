# panrna
PanRNA is a python written pipeline to map and count reads of an RNA-seq experiment. We divide the pipeline in two parts. The first one, we call Basic Run, is the the most common task on RNA-seq which only map the reads to a reference genome and return counts files for each sample. These can be used for differential expression analysis. The second part, the full PanRNA run, includes an extra task meant to avoid loosing information due to a lack in the reference genome annotation or due to the presence of genes in the samples that are not in the reference genome used.

For more details see the [manual](https://github.com/Mimop/panrna/blob/master/PanRNA_manual.pdf).


![PanRAN](https://github.com/Mimop/panrna/blob/master/panrna.jpg)
