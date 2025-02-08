# MotAli
Alignment of motifs to detect their similarity and co-occurrence

# Description
De novo motif search is the standard approach to detect nucleotide context specificity of transcription factor (TF) binding sites (BSs) in genomic DNA based on massive sequencing techniques such as ChIP-seq experiments. Although the standard motif model position weight matrix (PWM) is widely aplied, it can not represent all potential sites of direct TF binding. This conclusion is the consequence of hypothesis of the PWM model on the independence from the impacts of different positions of BSs. Thus, de novo motif search should not be restricted by only the PWM model of TF BS, it can apply [BaMM](https://github.com/soedinglab/BaMMmotif2) or [SiteGA](https://github.com/parthian-sterlet/sitega). To proceed with the results of de novo motif search, the lists of enriched motifs, one should apply a tool for motif comparison, e.g. TomTom is an example for the PWM model ([Gupta et al., 2007](https://doi.org/10.1186/gb-2007-8-2-r24)). However, the motif comparison tools currently are developped only for the PWM model. This drawback is very limiting for any motif model besides the PWM. That is why it is required to develop the motif comparison tool which is not depend on the motif model itself, so that the comclusion on the similarity of two motifs is approved or rejected by the degree of overlap between their predictions for a given test set of DNA sequences. The principle of motif comparison based on the similarity of their predictions for the given "full" set of possible DNA sequences were earlier proposed and applied in the tools of [Vorontsov et al. (2013)](https://doi.org/10.1186/1748-7188-8-23) and [Lambert et al. (2017)](https://doi.org/10.1093/bioinformatics/btw489), but these tools were developped for the PWM motif model only.

# Algorithm

# Input data
- FASTA file to recognize two motifs
- two motifs, currently only PWM and SiteGA motif models are supported. Each motif is either (a) the position frequency matrix required to compute PWM, or (b) SiteGA model of motif, the list of locally positioned dinucleotides with their weights, [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545).
- the threshold of the expected recognition rate for both motifs, [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545)
- two tables "Thresholsd vs. Expected Recognition Rate" computed preliminary with the whole genome file of gene promoters, [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545)

# Output data
- Histogram of AUPRC values as function of the mutual orientation of two motifs (they either in the same or opposite strands) and shift between the centers of two motifs
- File with the list of PR curves for all possible orientation and shifts
- Best AUPRC values for (a/b) the mild/stringent criteria on the overlap of the centers of two motifs, and (c) any possible shift beween centers of two motifs
- 
# Source code and command line arguments

# Command line examples
