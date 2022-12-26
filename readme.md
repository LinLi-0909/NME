## run NME

```bash
cd /path/to/NME/Code
filename=../Data/NME_input_data.csv
# the number of mapping neighbors
p=15
# the number of estimating neighbors
k=5
matlab -nodesktop -nosplash -r "filename='$filename',p=$p,k=$k;calc_NME;quit"
```

## run cNME

```bash
filename=../Data/cNME_input_data.tsv
# save cNME output to
output='../Output/example_cNME_output'
# the nubmer of conditional variables
nz=1
# the number of mapping neighbors
p=5
# the number of estimating neighbors
k=5
matlab -nodesktop -nosplash -r "filename='$filename',nz=$nz,p=$p,k=$k,output='$output';calc_cNME;quit"
```

## run scNME

```bash
filename=../Data/scNME_input_data.csv
# save cNME output to
output='../Output/example_scNME_output'
# the number of mapping neighbors
p=15
# species ='human' or 'mouse'
species='human'
# ncpu for paralell calculating
nworkers=20
# set the pcutoff of Hypothesis testing
p_cutoff=0.05
# save scNME output to
output='../Output/example_scNME_output'
matlab -nodesktop -nosplash -r "filename='$filename',species='$species',nz=$nz,p=$p,k=$k,p_cutoff=$p_cutoff,output='$output',nworkers=$nworkers;calc_scNME;quit"
```

