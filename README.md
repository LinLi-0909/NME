# scNME
Single-cell causal network inference based on neighbor cross-mapping entropy
## Overview
## Contents

- [Overview](#overview)
- [NME analysis Guide](./LICENSE)
- [scNME analysis Guide](./LICENSE)
- [License](./LICENSE)
- [Citation](#citation)
- [Contact](#Contact)

## Overview
Gene regulatory networks (GRNs) reveal the complex molecular interactions that govern cell state. However, it is challenging for identifying causal 
relations among genes due to noisy data and molecular nonlinearity. Here, we propose a novel causal criterion, neighbor cross-mapping entropy (NME) 
for inferring GRNs from both steady data and time-series data. NME is designed to quantify “continuous causality” or 
functional dependency from one variable to another based on their function continuity with varying neighbor sizes. 
NME shows superior performance on benchmark datasets, comparing with existing methods. 
By applying to scRNA-seq datasets, NME not only reliably inferred GRNs for cell types but also identified cell states. 
Based on the inferred GRNs and further their activity matrices, NME showed better performance in single-cell clustering and downstream analyses. 
In summary, based on continuous causality, NME provides a powerful tool in inferring causal regulations of GRNs between genes from scRNA-seq data, 
which is further exploited to identify novel cell types/states and predict cell type-specific network modules. <br /> 
![image](https://user-images.githubusercontent.com/63344240/209491331-f360e1a5-786e-48c6-b5ce-39e958373e95.png)
## NME analysis Guide
NME is designed based on matlab.

