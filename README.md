# MotAli
Alignment of transcription factor binding sites motifs to detect their similarity and co-occurrence

# Description
De novo motif search is the standard approach to detect nucleotide context specificity of transcription factor (TF) binding sites (BSs) in genomic DNA based on results of massive sequencing techniques such as ChIP-seq experiments. Although the standard motif model position weight matrix (PWM) is widely aplied, it cannot represent all potential sites of direct TF binding. This follows from the hypothesis of the PWM model on the independence from the impacts of different positions of a motif. De novo motif search is not restricted by only the PWM model of TF BS, it can apply the alternative motif models such as [BaMM](https://github.com/soedinglab/BaMMmotif2) or [SiteGA](https://github.com/parthian-sterlet/sitega) extending potential genomic loci of direct binding for tested target TFs ([Tsukanov et al., 2022](https://doi.org/10.3389/fpls.2022.938545)). To proceed with the results of de novo motif search, the list of enriched motifs, one should apply a tool for motif comparison. Unfortunately, available motif comparison tools have been developped only for the PWM model. Hence, de novo motif search application is very limited for any motif model besides the standard model PWM. The standard approach here is to apply the PWM model and to compare any possible enriched motifs with the motifs for binding sites of known transcription factors, e.g. from [JASPAR](https://jaspar.elixir.no/), [Hocomoco](http://hocomoco12.autosome.ru/) or [CisBP](https://cisbp.ccbr.utoronto.ca/) databases of PWM models for DNA binding motifs. The default option to detect the similarity of two motifs of the PWM model is to compare two respective position frequency matrices (PFMs) directly evaluating pairwise similarity of columns of these PFMs, e.g. ([Pietrokovski, 1996](https://doi.org/10.1093/nar/24.19.3836); [Sandelin and Wasserman, 2004](https://doi.org/10.1016/j.jmb.2004.02.048)). Thus, TomTom ([Gupta et al., 2007](https://doi.org/10.1186/gb-2007-8-2-r24)) is a popular tool for the motifs of PWM model. Although a PWM model can be constructed according to predictions of an arbitrary alternative model, and, consequently, the similarity between a PWM and an alternative model can be estimated, this approach is not quite correct, since it certainly destroys partially the sequence specificity of an alternatve model. The similarity of two arbitrary motif model can be assessed indirectly. An arbitrary motif model is a set of weights for certain position-specific nucleotide contexts like frequencies of k-mers. The indirect motif comparison tool is not depend on the exact motif model itself, it depends from profiles of predicted sites of both commpared motif models for a given test set of DNA sequences. Let we have two motifs, each one represents either the standard or alternative model of motif. This principle of indirect motif comparison based on the overlap of their predictions for a given 'full' set of possible DNA sequences were earlier proposed and implemented in the tools of [Vorontsov et al. (2013)](https://doi.org/10.1186/1748-7188-8-23) and [Lambert et al. (2017)](https://doi.org/10.1093/bioinformatics/btw489), but these tools were proposed for the PWM motif model only. 

We suggest the tool to estimate the similarity of two motifs and to predict their co-occurrence. The motif comparison is derived from prediction profiles of these motifs for given FASTA file of DNA sequences. All possible 'close' locations of two motifs with respect to their mutual orientation and the shifts between the centers of motifs are compared with slightly more distant locations of these two motifs according to the statistics of the Area under Precision-Recall curve (AUPRC). Although Precision-Recall plots are less frequently used than more popular Receiver Operating Characteristics (ROC) plots, it was proposed that the PR curve is more accurate than ROC curve for imbalanced datasets ([Saito, Rehmsmeier, 2015](https://doi.org/10.1371/journal.pone.0118432)).

# Algorithm
## Profiles of recognized sites for two input motifs
The tool takes a tested FASTA file, and two motifs. Among these two motifs motifA/motifB have the longer/shorter lengths. If Length(motifA) >= Length(motifB) then the longer/shorter ones are motifA/motifB, othewise vice versa. The input parameter of ERR (Expected Recognition rate) is used transforms the recogntion thresholds of both motifs to the uniform scale. Following previously developed approaches ([MCOT, Levitsky et al., 2019](https://doi.org/10.1093/nar/gkz800)) and ([MetArea, Levitsky et al., 2024](https://doi.org/10.18699/vjgb-24-90)), for each motif ERR is defined as the fraction of all tested positions in the whole-genome set of promoters of protein-coding genes that contain predicted sites. The tables 'Threshold vs. ERR' are preliminary computed for each of input motifs. Each table includes two columns: the first one is the recognition threshold given by a motif recognition model and the second one is the common logarithm of the probability value, -Log<sub>10</sub>(ERR). The lines of tables are sorted from the hihgest to the lowest threshold.
## Alignment of two motifs: DNA strands and shifts of motif centers relative each other
Furhther, the pairs of co-occured motifs with respect to DNA strands are predicted with the parameter of the maximal shift L between the centers of two motifs. There are two types of motif orientation. If two motifs lie in the same DNA strand, than their orientation is 'Direct': pairs in forward and reverse strands are '-> ->' and '<- <-'. Two motifs from the opposite DNA strands refer to Inverted and Everted orientations of two motifs, they are '-> <-' and '<- ->', respectively.  For each motif, its center is defined as the average value of its start and end positions. For example, the motifs of lengths 4 and 5 located at positions 1..4 and 1..5, the centres are located at positions 2.5 and 3. For given parameter of the maximal shift of motif centers, we allow the positive and negative shifts: the pairs of motifs 'motifA .. motifB' / 'motifB .. motifA' respect positive/negative shifts of centers of motifs. 
## Definition of metrics of True Positives and False Positives for PR curve
Counts of predicted pairs of motifs with certain orientation (Direct or Inverted/Everted) and shifts mean the metrics of true positives (TP) predictions. Corresponding pairs mapped in any orientation with any larger shifts from L to 2L are used as the control metrics of false positives (FP). These TP and FP metrics are computed as a function of the ERR values. FP metrics normilized to the TP metrics since TP are defined for each shift and orientation, i.e. FP values transformed FP / 4L values (each motif is predicted with left/right section in any of two strands at maximal distance of L bp). Finally, for each shift and orientation the area under Precision-Recall (PR) curve estimates the alignment of two motifs. Recall is defined as a ratio of the number of predicted pairs (TP) at certain recognition threshold (ERR) to the number of predicted pairs (TP) at the lowest allowable threshold, default ERR = 0.002. Precision is defined as the portion of true positives among the total number of true positives and false positives, Precision = TP/(TP + FP). The value AUPRC<sub>MAX</sub>(motifA,motifB) is defined for given pair of motifs as the maximal AUPRC value among all shifts and both orientations. All shifts are classified into (a) the strigent overlap (motif center of the shorter motif lies inside the longer motif), and (b) all the remaining cases when the centers of two motifs lie within the threshold of shift L. The highest possible value of the area under PR curve measure (AUPRC) of 1 implies perfectcly matching motifs (i.e. two motifs are coincides). The AUPRC values of 0.5 refers to two random motifs. 
## Definition of similarity score for two motifs
For the descibed above metrics AUPRC<sub>MAX</sub>(motifA,motifB) respecting 'heterotypic pair' of motifs 'motifA & motifB', two analogous normalizing values AUPRC<sub>MAX</sub>(motifA,motifA) and AUPRC<sub>MAX</sub>(motifB,motifB) estimate for two homotypic pairs, 'motifA & motifA', and 'motifB & motifB', the overlap of each motif with itself. Finally, the similarity score of two input motifs is computed as follows: 

__Similarity score = AUPRC<sub>MAX</sub>(motifA,motifB) / Max{AUPRC<sub>MAX</sub>(motifA,motifA), AUPRC<sub>MAX</sub>(motifB,motifB)}__

Typically, the similarity score above 0.999 / 0.98 / 0.95 / 0.9 and 0.7 respect the excellent / very high / high / moderate and low levels of similarity between two motifs. A value of 0.5 means that the overlap of two motif  and its absence have equal probabilities (motifs are very dissimilar), and the zero value implies that the overlap of motifs is forbiden with the input threshold of motif recognition. For example, two quite similar yet distinct E-box motifs for the BHLHA15 TF from [Hocomoco](https://doi.org/10.1093/nar/gkad1077),   and [BHA15.H12CORE.0.P.B](https://hocomoco12.autosome.org/motif/BHA15.H12CORE.0.P.B) and [BHA15.H12CORE.1.SM.B](https://hocomoco12.autosome.org/motif/BHA15.H12CORE.1.SM.B) have consensuses CAgcTG and CAtaTG, respectively. Their similarity score is 0.998. The [TomTom tool](https://meme-suite.org/meme/tools/tomtom) reports the p-value of 1.46e-03 for the significance of similarity of these motifs.

# Input data and parameters

Following input data are required:

- Two motifs, only PWM and SiteGA motif models are supported. The PWM motif is represented by the position frequency matrix in [standard format](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm), and the SiteGA model is the [list of locally positioned dinucleotides with their weights](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1), [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545).
- FASTA file to recognize these two motifs, e.g. ChIP-seq peaks, see example format here [top 1000 peaks for mouse BHLHA15 TF](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa)
- Two tables 'Threshold vs. ERR (Expected Recognition Rate)', see example format here [ERR table for the motif of mouse BHLHA15 TF](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.dist). These tables computed preliminary with the FASTA file of promoters of protein-coding genes of whole genome, [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545), example FASTA file for [mouse](https://github.com/parthian-sterlet/mcot-kernel/blob/master/genomes/mm/ups2kb_mm10.seq.tar.gz)

The tool needs two parameters:
- The threshold of the Expected Recognition Rate (ERR) for both motifs, [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545)
- The maximal shift between centers of two motifs, L
  
# Output data
- The similarity scores and the maximal AUPRC values for (a) the stringent criterion on the overlap of the centers of two motifs, and (b) any shifts between the centers of two motifs, including their location with spacers and small overlaps. 
- Histograms of AUPRC values as a function of the mutual orientation of two motifs (they either in the same or opposite strands) and the shift between their centers. Three histograms show distributions for the pair of input motif (heterotypic case, motifA/motifB), and two separate distributions for the first and second motifs (two homotypic cases, motifA/motifA and motifB/motifB). Due to possible even/odd lengths of motifs in a pair, positions in histograms are marked either by integer or non-integer values, ..., -2, -1, 0, +1, +2, ... or ..., -1.5, -0.5, +0.5, +1.5, ..., respectively.
- PR curves for two possible mutual orientations of motifs (motifs in the same and opposite strands) and shifts, -L <= x <= L, x denotes the positive shift means that the longer/shorter motifs are located closer to 5'/3' 
 ends of DNA sequence alignment. Three blocks of data are PR curves for the pair of input motif (heterotypic case), and two separate blocks of curves for the first and second motifs (two homotypic cases).

# Source code

[mco_prc.cpp](https://github.com/parthian-sterlet/motali/blob/main/src/mco_prc.cpp)

# Command line arguments

1. input FASTA file, set of DNA sequences, ChIP-seq peaks, example - [top 1000 peaks (PEAKS039234) for mouse BHLHA15 TF](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) from [GTRD](https://gtrd.biouml.org/#!)
2. type of the first motif model, values 'pwm'/'PWM' and 'sga'/'SGA' imply PWM and SiteGA motif models
3. type of the second motif model, values 'pwm'/'PWM' and 'sga'/'SGA' imply PWM and SiteGA motif models
4. input file for the first motif, either of [PWM](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm) or [SiteGA](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) model
5. input file for the second motif, either of [PWM](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm) or [SiteGA](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) model
6. input file 'Threshold vs. ERR' for the first motif, this file has the same format for PWM and SiteGA models, [ERR table](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.dist). These files for PWM and SiteGA models are generated by [PWM thresholds selection program](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pwm_thr_err/pwm_iz_pwm_thr_dist0.cpp) from [MCOT](https://github.com/parthian-sterlet/mcot-kernel) and [SiteGA thresholds selection program](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) from [SiteGA](https://github.com/parthian-sterlet/sitega) repositories.
7. input file 'Threshold vs. ERR' file for the second motif, this file has the same format for PWM and SiteGA models, [ERR table](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.dist). These files for PWM and SiteGA models are generated by [PWM thresholds selection program](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pwm_thr_err/pwm_iz_pwm_thr_dist0.cpp) from [MCOT](https://github.com/parthian-sterlet/mcot-kernel) and [SiteGA thresholds selection program](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) from [SiteGA](https://github.com/parthian-sterlet/sitega) repositories.
8. double ERR threshold, the default value is 0.002, this value (maximal Expected Recognition Rate) means the lowest threshold used to predicts site and compute the area under PR curve. This ERR value defines the range of recognition thresholds from the table 'Threshold vs. ERR'.
9. integer value, L, maximal shift between the centers of the first and second motifs, the default value is 50 bp
10. output file, histograms of AUPRC values for both orientations of motifs all shifts of their centers, [Histogram](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.hist)
11. output file, PR curves for both orientations of motifs all shifts of their centers, [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.prc)
12. output file, text file of similarity scores and the maximal AUPRC values for the stringent overlap and any shift between the centers of two motifs, [Similarity scores and AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.sta)
13. integer value, 0/1 mean printing and non-printing of histograms of AUPRC values to the output file (see argument #10 above)
14. integer value, 0/1 mean printing and non-printing of PR curves to the output file (see argument #11 above)
15. integer value, 0/1 mean short/detailed description of results in the output file of similarity scores and the maximal AUPRC values (see argument #12 above)

# Command line examples
- One PWM motif vs. another PWM motif

[com_line_pwm_pwm](https://github.com/parthian-sterlet/motali/blob/main/run/com_line_pwm_pwm) command line for two motifs of the PWM model [BHA15.H12CORE.0.P.B.pfm](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B.pfm) and [BHA15.H12CORE.1.SM.B.pfm](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.1.SM.B.pfm) respecting two structurally distinct motif types for [murine BHLHA15 TF](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:22265) from [Hocomoco](http://hocomoco12.autosome.ru/) and [ChIP-seq peaks dataset PEAKS039234](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) for this TF from [GTRD](https://gtrd.biouml.org/#!). Here and below parameters 0.002 and 50 are used for the threshold of ERR and the maximal shift. The output data are [Histogram of AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B_1.SM.B.hist), [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B_1.SM.B.prc) and [Similarity scores and maximal AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B_1.SM.B.sta).

- PWM motif vs. SiteGA motif

[com_line_pwm_sga](https://github.com/parthian-sterlet/motali/blob/main/run/com_line_pwm_sga) command line for two motifs of PWM and SiteGA models. The PWM motif [BHA15.H12CORE.0.P.B.pfm](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm) respects the major structural motif type for [murine BHLHA15 TF](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:22265) from [Hocomoco](http://hocomoco12.autosome.ru/) and [ChIP-seq peaks dataset PEAKS039234](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) for this TF from [GTRD](https://gtrd.biouml.org/#!). The SiteGA motif [PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) is the first ranked enriched motif derived by [SiteGA de novo motif search](https://github.com/parthian-sterlet/sitega) in the 1st iteration of the cross-validation procedure. The output data are [Histogram of AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.hist), [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.prc) and [Similarity scores and maximal AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.sta).

- One SiteGA motif vs. another SiteGA motif

[com_line_sga_sga](https://github.com/parthian-sterlet/motali/blob/main/run/com_line_sga_sga) command line for two motifs of the SiteGA model. These models respect the [murine BHLHA15 TF](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:22265) and [dataset of top 1000 ChIP-seq peaks (PEAKS039234)](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) for this TF from [GTRD](https://gtrd.biouml.org/#!). Two motifs [PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) and [PEAKS039234_BHLHA15_Q9QYC3_MACS2_2_40_cmat1](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_2_40_cmat1) are the first ranked enriched motifs derived by [SiteGA de novo motif search](https://github.com/parthian-sterlet/sitega) in the 1st and 2nd iterations of the cross-validation procedure. The output data are [Histogram of AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/sga_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_i12_cmat_1_1.hist), [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/sga_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_i12_cmat_1_1.prc) and [Similarity scores and maximal AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_i12_cmat_1_1.sta).
