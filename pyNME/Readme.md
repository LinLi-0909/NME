## Step1: Prepare expression data

For any normalized gene expression matrix, we select the expression of TF genes and 2000 HVGs, and save the matrix in \*.txt/\*.tsv file.

```R
library(tidyverse)
exp = read.csv("pbmc_exp.csv", header=T) # your expression matrix
genes = exp$X
rownames(exp) = genes
exp = select(exp, -X)

row_variances <- apply(exp_filt, 1, var)
top_2000_rows <- order(row_variances, decreasing = TRUE)[1:2000]
all_genes = rownames(exp_filt)

tf = read.table("hs_tf.txt")
tf_gene = unique(tf$V1)
tf_pos <- match(tf_gene, all_genes)
tf_pos = na.omit(tf_pos)
keep_rows = unique(c(top_2000_rows, tf_pos))
exp_filt <- exp_filt[keep_rows, ]

write.table(exp_filt, "pbmc_filtered_exp.txt", sep = '\t', quote = F) 
```

## Step2: Run GRN

The C++ source code is in `./src/grn.cpp`,  and the compiled software is in `./bin/grn`, which can be run directed on Linux system. To compile the source code manually, please be sure that tou have already installed `mlpack`  package.

```bash
# macOC
brew install mlpack
g++ -o ./bin/grn ./src/grn.cpp -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib
# Anaconda
conda install mlpack
brew install mlpack
g++ -o ./bin/grn ./src/grn.cpp -std=c++17 -I~/anaconda3/envs/mlpack/include -L~/anaconda3/envs/mlpack/lib
```

The sample running command can be

```bash
./src/grn --cm 10 --mi 5 --exp_file pbmc_filtered_exp.txt --species hs --computed_genes 2000 --output_grn ./pbmc_grn.csv
```

You can view the usage of this command by `./src/grn -h`

## Step3: Calculate network entropy and plot

Please follow the `plot_func.ipynb`  in `plot`.