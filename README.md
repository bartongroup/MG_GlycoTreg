# Glycosyltransferases and glycosidases in Treg activation

Project collaborators: Paul Crocker, Gavuthami Murugesan, Gang Wu

## Description

Analysis of RNA-seq data in mouse in two biological conditions, with regulatory T cells (Tregs) at rest and active.

## Instructions


First, create and activate a conda environment:

```
conda create --name glycotreg --file conda-spec.txt
source activate glycotreg
```

Then, use snakemake to perform all the necessary computations on the cluster:

```
snakemake -c "qsub -V -cwd -o snakelog -e snakelog -pe smp {threads}" --jobs=100
```

Next, edit the `R/setup.R` file and modify the home project directory (in our case it is on a remote system). Finally, knit the main document. 

```
rmarkdown::render("doc/analysis.1.Rmd")
```

The final report with model description and results is in the file `doc/analysis.1.html`.
