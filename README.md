# MotAli
Alignment of motifs to detect their similarity and co-occurrence

# Description
De novo motif search is the standard approach to detect nucleotide context specificity of transcription factor (TF) binding sites (BSs) in genomic DNA based on massive sequencing techniques such as ChIP-seq experiments. Although the standard motif model position weight matrix (PWM) is widely aplied, it can not represent all potential sites of direct TF binding. This conclusion is the consequence of hypothesis of the PWM model on the independence from the impacts of different positions of BSs. Thus, de novo motif search should not be restricted by only the PWM model of TF BS, it can apply [BaMM](https://github.com/soedinglab/BaMMmotif2) or [SiteGA](https://github.com/parthian-sterlet/sitega). To proceed with the results of de novo motif search, the lists of enriched motifs, one should apply a tool for motif comparison, e.g. TomTom is an example for the PWM model ([Gupta et al., 2007](https://doi.org/10.1186/gb-2007-8-2-r24)). However, the motif comparison tools currently are developped only for the PWM model. This drawback is very limiting for any motif model besides the PWM. That is why it is required to develop the motif comparison tool which is not depend on the motif model itself, so that the comclusion on the similarity of two motifs is approved or rejected by the degree of overlap between their predictions for a given test set of DNA sequences. The principle of motif comparison based on the similarity of their predictions for the given "full" set of possible DNA sequences were earlier proposed and applied in the tools of [Vorontsov et al. (2013)](https://doi.org/10.1186/1748-7188-8-23) and [Lambert et al. (2017)](https://doi.org/10.1093/bioinformatics/btw489), but these tools were developped for the PWM motif model only.

We propose the MotAli tool to compare two motifs representing potential TF BS. The motif comparison is derived from prediction profiles of these motifs for given FASTA file of DNA sequences. All possible close

# Algorithm
Input data of the teted FASTA file, and two motif models are used to predict all co-occured motifs with respect to DNA strands. Co-occ
Motif orientations and shift of motif cetners
Two motifs may lie in the same DNA strand, othewise they lie in the opposite strand. The first case 'Direct' forward and reverse strands represent pairs of motifs '-> ->' and '<- <-', the second means the Inverted ans Everted orientations of two motifs, '-> <-' and '<- ->', respectively. All these cases are formilized by pairs of motifs 'motifA .. motifB' / 'motifB .. motifA' respecting positive/negative spacers. 

# Input data and parameters

MotAli tool requires the following input data:

- FASTA file to recognize two motifs
- two motifs, only PWM and SiteGA motif models are supported. PWM motif is represented by the position frequency matrix in a [standard format](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm), and the SiteGA model is the [list of locally positioned dinucleotides with their weights](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1), [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545).
- two tables "Threshold vs. Expected Recognition Rate" computed preliminary with the FASTA file of promoters of protein-coding genes of whole genome, [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545), example FASTA file for [mouse](https://github.com/parthian-sterlet/mcot-kernel/blob/master/genomes/mm/ups2kb_mm10.seq.tar.gz)

MotAli tool has two parameters:
- the threshold of the Expected Recognition Rate (ERR) for both motifs, [(Tsukanov et al., 2022)](https://doi.org/10.3389/fpls.2022.938545)
- maximal shift between centers of two motifs, L
  
# Output data
- Histogram of AUPRC values as a function of the mutual orientation of two motifs (they either in the same or opposite strands) and the shift between the centers of two motifs
- List of PR curves for two possible mutual orientations of motifs (Direct and Opposite strands) and shifts, 0 <= x <= L.
- Best AUPRC values for (a/b) the mild/stringent criteria on the overlap of the centers of two motifs, and (c) any possible shift beween centers of two motifs

# Source code

[mco_prc.cpp](https://github.com/parthian-sterlet/motali/blob/main/src/mco_prc.cpp)

# Command line arguments

1. input FASTA file, set of DNA sequences, ChIP-seq peaks, example - [top 1000 peaks for mouse BHLHA15 TF](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) from [GTRD](https://gtrd.biouml.org/#!)
2. type of the first motif model, values 'pwm'/'PWM' and 'sga'/'SGA' imply PWM and SiteGA motif models
3. type of the second motif model, values 'pwm'/'PWM' and 'sga'/'SGA' imply PWM and SiteGA motif models
4. input file for the first motif, either of [PWM](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm) or [SiteGA](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) model
5. input file for the second motif, either of [PWM](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm) or [SiteGA](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) model
6. input file 'Threshold vs. ERR' for the first motif, the same format for PWM and SiteGA models, [ERR table](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.dist). These files for PWM and SiteGA models are generated by [PWM thresholds selection program](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pwm_thr_err/pwm_iz_pwm_thr_dist0.cpp) from [MCOT](https://github.com/parthian-sterlet/mcot-kernel) and [SiteGA thresholds selection program](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) from [SiteGA](https://github.com/parthian-sterlet/sitega) repositories.
7. input file 'Threshold vs. ERR' file for the second motif, the same format for PWM and SiteGA models, [ERR table](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.dist). These files for PWM and SiteGA models are generated by [PWM thresholds selection program](https://github.com/parthian-sterlet/mcot-kernel/blob/master/src/pwm_thr_err/pwm_iz_pwm_thr_dist0.cpp) from [MCOT](https://github.com/parthian-sterlet/mcot-kernel) and [SiteGA thresholds selection program](https://github.com/parthian-sterlet/sitega/blob/master/src/sitega_thr_dist_mat.cpp) from [SiteGA](https://github.com/parthian-sterlet/sitega) repositories.
8. double ERR threshold, the default value is 0.01, this value (maximal Expected Recognition Rate) means the most mild threshold used to predict site and compute the area under PR curve. This ERRmax value define the range of recognition thresholds from the table 'Threshold vs. ERR'.
7. integer value, L, maximal shift between the centers of the first and second motifs, the default value is 50 bp
8. output file, histogram of AUPRC values for both orientations of motifs all shifts of their centers, [Histogram](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.hist)
9. output file, PR curves for both orientations of motifs all shifts of their centers, [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.prc)
10. output file, text file for the best AUPRC values, [Common statistics](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.sta)

# Command line examples
- One PWM motif vs. another PWM motif

[com_line_pwm_pwm](https://github.com/parthian-sterlet/motali/blob/main/run/com_line_pwm_pwm) command line for two PWM models [BHA15.H12CORE.0.P.B.pfm](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B.pfm) and [BHA15.H12CORE.1.SM.B.pfm](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.1.SM.B.pfm) respecting two structurally distinct motif types for [murine BHLHA15 TF](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:22265)  from [Hocomoco v12](http://hocomoco12.autosome.ru/) and [ChIP-seq peaks dataset PEAKS039234](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) for this TF. Here and below parameters 0.01 and 30 are used for the threshold of ERR and shift. The output data are [Histogram of AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B_1.SM.B.hist), [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B_1.SM.B.prc) and [Best AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_pwm/BHA15.H12CORE.0.P.B_1.SM.B.sta).

- PWM motif vs. SiteGA motif

[com_line_pwm_sga](https://github.com/parthian-sterlet/motali/blob/main/run/com_line_pwm_sga) command line for PWM and SiteGA models. The PWM motif [BHA15.H12CORE.0.P.B.pfm](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B.pfm) respects the major structural motif type for [murine BHLHA15 TF](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:22265) from [Hocomoco v12](http://hocomoco12.autosome.ru/) and [ChIP-seq peaks dataset PEAKS039234](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) for this TF. The SiteGA motif [PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) is the enriched motif derived by SiteGA de novo motif search in the 1st iteration of cross validation procedure.  The output data are [Histogram of AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.hist), [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.prc) and [Best AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/BHA15.H12CORE.0.P.B_PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_cmat1.sta).

- one SiteGA motif vs. another SiteGA motif

[com_line_sga_sga](https://github.com/parthian-sterlet/motali/blob/main/run/com_line_sga_sga) command line for two SiteGA models. These models respect the [murine BHLHA15 TF](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:22265) and [ChIP-seq peaks dataset PEAKS039234](https://github.com/parthian-sterlet/motali/blob/main/examples/PEAKS039234_BHLHA15_Q9QYC3_MACS2.fa) for this TF. The SiteGA motifs [PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_1_40_cmat1) and [PEAKS039234_BHLHA15_Q9QYC3_MACS2_2_40_cmat1](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_2_40_cmat1) are the enriched motifs derived by SiteGA de novo motif search in the 1st and 2nd iterations of the cross validation procedure.  The output data are [Histogram of AUPRC values]([PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_i12_cmat_1_1.hist](https://github.com/parthian-sterlet/motali/blob/main/examples/sga_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_i12_cmat_1_1.hist)), [PR curves](https://github.com/parthian-sterlet/motali/blob/main/examples/sga_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_i12_cmat_1_1.prc) and [Best AUPRC values](https://github.com/parthian-sterlet/motali/blob/main/examples/pwm_sga/PEAKS039234_BHLHA15_Q9QYC3_MACS2_40_i12_cmat_1_1.sta).
