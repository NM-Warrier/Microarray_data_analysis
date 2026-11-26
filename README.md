# Consensus DEG and Exploratory Machine Learning Pipeline

A general single-script R workflow for directional consensus DEG identification, microarray preprocessing, gene-expression matrix preparation, and exploratory binary or multiclass classification using Random Forest and Elastic Net.

## Requirements

- R
- Required packages: `readxl`, `openxlsx`, `dplyr`, `stringr`, `Biobase`, `affy`, `sva`, `ggplot2`, `randomForest`, `glmnet`, `pROC`, `tidyr`
- Optional package: `pheatmap`

Install packages in R:

```r
install.packages(c(
  "BiocManager", "readxl", "openxlsx", "dplyr", "stringr",
  "ggplot2", "randomForest", "glmnet", "pROC", "tidyr", "pheatmap"
))
BiocManager::install(c("Biobase", "affy", "sva"))
```

## Inputs

| Input | Required format |
| --- | --- |
| DEG results | One Excel workbook per comparison, with one iteration per worksheet. |
| Combined raw microarray object | An RDS containing an `AffyBatch` object from a single compatible array platform. |
| Sample metadata | Excel file containing `SampleID`, `Group`, and `Batch`. Additional outcome columns may be included. |
| Probe annotation | CSV containing `ProbeID` and `GeneSymbol`, with one row per probe ID. |

DEG worksheets need a gene identifier, log2 fold change, and adjusted P value. Supported column aliases include:

| Field | Examples |
| --- | --- |
| Gene identifier | `GeneSymbol`, `gene_symbol`, `Gene`, `SYMBOL`, `gene_id`, `ProbeID` |
| Log2 fold change | `logFC`, `log2FC`, `log2FoldChange`, `estimate` |
| Adjusted P value | `adj.P.Val`, `padj`, `p.adjust`, `FDR` |

Use `estimate` only when it represents log2 fold change. All iteration sheets within a comparison must use the same contrast direction and identifier type. Include iteration sheets only, without summary tables.

For the complete workflow, DEG identifiers must be gene symbols matching the annotated expression matrix. Stable IDs or probe IDs can be used for consensus-only work, but must be mapped consistently before the gene-symbol matrix and ML stages.

Metadata sample IDs must match the raw object's sample names exactly. Supply annotation appropriate to the actual platform; resolve ambiguous mappings before running. This preprocessing workflow is intended for compatible microarrays, not RNA-seq count data or cross-platform pooling.

## Usage

1. Open `general_consensus_and_ML_pipeline.R`.
2. Replace working-directory, input-directory, and input-file placeholders.
3. Update `deg_files` and the significance thresholds.
4. Set `groups_to_keep` to the metadata groups to retain.
5. Review annotation and any independently justified exclusions.
6. Edit `ml_jobs` to specify the consensus feature set, outcome column, and classes for each model.
7. Review repeat counts, training fraction, and model settings.
8. Run the script:

```bash
Rscript general_consensus_and_ML_pipeline.R
```

Each `feature_set` must match the corresponding DEG workbook name after removing `_DEGs_AllIterations.xlsx`. Separate class labels with semicolons in `ml_jobs`. For binary outcomes, list the positive class last. To merge original groups, create an explicit outcome column in the metadata.

The full script runs all stages sequentially. To run consensus identification alone, run only the first section through its completion message.

## Workflow

1. Filter each iteration's DEGs using adjusted P value and absolute log2 fold change.
2. Identify directionally consistent genes and strict consensus genes.
3. Subset the combined raw microarray object and perform joint RMA normalization.
4. Generate PCA plots and flag unusual PC scores for review.
5. Apply ComBat batch correction while preserving the biological group effect.
6. Map probes to gene symbols and select the highest-variance probe per gene.
7. Create strict-consensus expression matrices.
8. Fit class-weighted Random Forest and Elastic Net models across repeated stratified holdouts.
9. Export performance, predictions, confusion matrices, and feature-selection frequencies.

## Consensus definitions

- **Directional consensus:** significant rows for an identifier agree in direction wherever that identifier passes the filters. An identifier may appear in only one iteration.
- **Strict directional consensus:** the identifier passes the filters in every iteration, with a consistent direction.

Nonsignificant rows do not participate in the directional check. Repeated identifiers count once per iteration, while reported means and medians summarize all retained rows. Minimum and mean adjusted P values are descriptive summaries, not combined significance tests.

## Quality control and batch effects

PCA flags are exported for review and do not automatically remove samples. Enter independently justified exclusions in `reviewed_exclusions` and rerun if needed.

Batch correction stops when batch and group effects are confounded. It does not fall back to a model that drops the group effect. With a single batch, correction is skipped. The current multiple-batch workflow requires at least two samples per batch.

## Exploratory ML interpretation

DEG selection, RMA normalization, representative-probe selection, and group-aware batch correction occur before train/test splitting. Therefore, reported performance can be optimistic and must not be presented as independent validation or deployment performance.

For validated prediction, perform learned preprocessing and feature selection within training resamples and evaluate on an independent cohort. Related samples require an appropriately grouped evaluation design.

The script reports descriptive summaries across overlapping repeated holdouts, without treating them as independent replicates or calculating confidence intervals. ROC plots use the first repeat rather than the best-performing repeat. Multiclass ROC uses micro-averaged one-versus-rest indicators.

Random Forest selection frequency counts top-ranked importance features. Elastic Net selection frequency counts nonzero coefficients. These definitions differ and should be considered when interpreting combined stability scores.

## Outputs

| Folder | Contents |
| --- | --- |
| `Directional_Consensus_DEGs` | Filtered DEG rows, directional consensus, and strict consensus tables. |
| `Preprocessing` | Raw subset, normalized and corrected matrices, metadata, PCA plots, and QC flags. |
| `Expression_matrices` | Representative-probe mapping, gene-expression matrix, strict-consensus matrices, and unmatched genes. |
| `Exploratory_ML` | Models, sample splits, predictions, ROC plots, confusion matrices, descriptive performance summaries, and feature-stability outputs. |

R session information is also saved. Use a fresh output directory for each analysis to avoid mixing results from different runs. Saved model objects correspond to individual resamples; the script does not produce a validated deployment model.

## Validation

The script has passed an R syntax check. Full execution and scientific validation require suitable inputs and installed packages.

## References

- [affy package](https://bioconductor.org/packages/affy/)
- [sva and ComBat](https://bioconductor.org/packages/sva/)
- [glmnet documentation](https://glmnet.stanford.edu/)
# Microarray_data_analysis
