# Human Vestibular Schwannoma Single Cell RNA Sequencing Analysis

GEO Series: [GSE216784](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216784)

For the scRNA-seq dataset GEO Series: [GSE216783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216783)

Title：N-Cadherin Dynamically Regulates Schwannoma Migration and Represents a Novel Therapeutic Target in Vestibular Schwannoma

Summary: 
scRNA-seq data from Barrett et al. (2023) was downloaded from NCBI’s Gene Expression Omnibus. scRNA-seq data were obtained as preprocessed count matrices. Initial quality control assessment confirmed that standard quality control had already been applied to the dataset; therefore, all cells from the original data were retained for downstream analyses. To minimize batch effects of input samples, we chose to only include fresh tumor tissue samples (n=11), including SCH1–6, SCH9, SCH13, and SCH14, yielding a total of 63,931 cells. All analyses were performed in Python using scanpy (v1.9.8). For downstream processing, we normalized the count matrix to total library size and applied a log normalization. We then selected the top 3,000 highly variable genes to perform PCA. Batch effects across the nine runs were corrected using the BBKNN algorithm, where the neighborhood graph was constructed based on the top 50 principal components. This integrated graph served as the input for UMAP visualization and Leiden clustering. Final cell type annotations were assigned by mapping each cell to its corresponding label defined in the original study’s metadata file, including Schwann cell subtypes (nmSC and myeSC) and stromal/immune populations. CDH2 expression was plotted as log2 normalized expression on the reduced dimensions UMAP plot. 

Contact Name: Taha Jan (taha.a.jan@vumc.org)

Organization: Vanderbilt University Medical Center
