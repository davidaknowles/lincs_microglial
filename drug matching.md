**Aim 3\. Prioritize drug candidates by matching genetic and chemical signatures** 

**Preliminary results.** 

In collaboration with Dr Rahul Satija, MPI Knowles has recently developed a novel computational framework that leverages single-cell pooled CRISPR perturbation screens as “dictionaries” to explain drug transcriptomic effects and predict mechanisms of action (**Figure 3.1a-b)**. In these screens, genes or combinations of genes can be systematically knocked down to determine the effects on gene expression at single-cell resolution 72–74. One challenge is that since many genes are closely related to one another (e.g., encoding members of the same complex), their perturbation responses can be highly correlated. We therefore use SuSiE to model drug response as a sum of genetic perturbation responses, which provides multiple credible sets of potentially explanatory genes (**Figure 3.1c-d**). We call this approach *RNA fingerprinting*.

**Approach**

We will primarily use the Broad Library of Integrated Network-Based Cellular Signatures (LINCS) due to the breadth of genetic and chemical perturbations and cell lines assayed. LINCS includes over 1M expression profiles using a custom array quantifying just under a thousand genes, which were selected to allow optimal imputation of genome-wide expression. 

A key challenge is that brain cell-types are not well represented in LINCS. However, there is substantial data for proxy cell lines including the myeloid/monocytic THP1 cell line (\>16k profiles), induced pluripotent stem cells (iPSC) neural progenitor cells and NGN-derived neurons (over 18k profiles total), and \>700 profiles from two astrocytoma lines. 

To augment the available data we will use recent advances in using machine learning to predict perturbation response for unseen drug/cell-line combinations. In particular, we will apply the compositional perturbation autoencoder (CPA), which, while primarily developed for single-cell data, was shown to work well with the LINCS data. CPA models observed expression using a latent space where the “basal state” (i.e., without any perturbation) is combined with learned perturbation (e.g., drug) and covariate (e.g., cell-type) vectors. 

Our method matches drug perturbation signatures to genetic knockdowns, accounting for uncertainty due to correlated profiles (blue drug) and multiple targets (green drug). **b**. We assess performance by matching genetic perturbations to each other: here between replicates of a Perturb-seq experiment targeting cleavage and polyadenylation factors75. Performance is still strong even across different cell lines.  **c/d.** The model has high accuracy distinguishing single vs dual gene knockdowns, even when relying on a single-gene knockdown dictionary.* 

Importantly, the perturbation latent vectors are shared across cell-types, so that “counterfactual” predictions can be made: what would the expression of this sample have been if treated with a different drug? CPA was previously applied to the 1000 most represented compounds in LINCS, all of which have at least 168 profiles. Using cross-validation we will assess whether it is possible to predict the effects of drugs with fewer measured profiles (for example, requiring only 100 observed profiles includes 1954 drugs). This will allow us to expand the repertoire of drugs we can consider, especially for the astrocyte proxies and microglia where few drugs were profiled. 

| Cell line | Description | Proxy for | Control profiles | Perturbation profiles | Drugs |
| :---- | :---- | :---: | :---: | :---: | :---: |
| NPC | iPSC neural progenitors | Neurons | 892 | 12634 | 3995 |
| NEU | iPSC NGN2 neurons | Neurons | 259 | 4452 | 2990 |
| THP1 | Monocytic leukemia | Microglia | 1159 | 15820 | 1764 |
| SHSY5Y | Neuroblastoma | Neurons | 108 | 1867 | 224 |
| U251MG | Astrocytoma | Astrocytes | 12 | 360 | 120 |
| YH13 | Astrocytoma | Astrocytes | 13 | 338 | 119 |
| MICROGLIA-PSEN1 | Microglia w PSEN1 mutation | Microglia | 1 | 9 | 9 |

***Table 3.1.** LINCS perturbation data for brain proxy cell types. The number of perturbation profiles is smaller than the number of drugs because multiple doses are tested.* 

We will apply three strategies of increasing sophistication to prioritize drugs. 

1. ***Gene-level prioritization***. For each significant (cell-type, gene) pair from the TWAS (Aim 2\) we will rank the drugs that result in maximal change in the risk-reducing direction in the corresponding proxy cell lines (Table 3.1). 

2. **Genome-wide correlation-based analysis**. We will correlate the vector of per cell-type TWAS effect sizes (Aim 2\) with per-drug perturbation vectors in the proxy cell lines. 

3. **Fingerprinting**. We will reverse our fingerprinting model (see Preliminary Results) to determine the drug (or combination of drugs) mostly likely to give a desired neuroprotective gene signature. The output of the model will be (potentially) multiple credible sets per cell-type. The drugs in each credible set can be interpreted as alternatives to one another (and will tend to have highly correlated effects). A result involving multiple credible sets will suggest combination therapy would be beneficial. 

While beyond the scope of this proposal, we look forward to working with our experimental collaborators to test the effects of the prioritized drugs in model systems. 

Potential pitfalls and alternative approaches. There are two potential pitfalls of this analysis: 1\) that the cell lines available in LINCS do not represent primary brain cell-types sufficiently well, and 2\) that protocol differences between RNA-seq and L1000 will make comparison difficult. For both these concerns, an alternative strategy is to leverage the extensive genetic perturbation data available in LINCS (beyond the genetic perturbations). Rather than directly matching the TWAS effect sizes, we will compare the LINCS genetic perturbation for each TWAS gene to the drug perturbations (in the same cell line). We will assess which drugs are the most frequently highly ranked. This way, expression data is only compared within the same protocol and cell-type. 