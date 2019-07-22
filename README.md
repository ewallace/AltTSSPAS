---
title: "README for AltTSSPAS Repository"
author: Edward Wallace
output: html
---

This repository is for estimates of alternate transcription start site and polyadenylation sites in fungi.

The question: does the relative use of different mRNA end sites change between conditions?

The stats work was started in February 2018. Picked up again in July 2019.

# Contents

## data

un-normalized counts of TSS and PAS usage in various conditions

## src

scripts to analyze data, mostly in .Rmd format

## results

estimation of differentially used termini

# To-do

    [x] add data
    [x] add .Rmd files from 2018
    [x] check original .Rmd files run on data
    [ ] collect functions (negbin fits etc) into single script (started)
    [ ] add alt PAS
    [ ] expand to zero-inflated model?
    [x] assess effect on start codon usage (alt N-termini)
    [ ] assess effect on protein localization 
    [ ] interactive data visualization of clusters
    [ ] .wig raw TSS track visualization