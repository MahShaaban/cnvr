#!/bin/bash

IMG="/data/scratch/DGE/DUDGE/MOPOPGEN/software-env/singularity/penncnv.1.0.5.img"
singularity exec $IMG bash

# call cnv on individuals
for i in offspring mother father; do
    detect_cnv.pl \
        -test \
        -hmm data-raw/PennCNV/lib/hhall.hmm \
        -pfb data-raw/PennCNV/example/example.pfb \
        --confidence \
        data-raw/PennCNV/example/${i}.txt \
        -log data-raw/output/${i}.log \
        -out data-raw/output/${i}.cnv
done

cat data-raw/output/{offspring,father,mother}.cnv > data-raw/output/family.cnv

scan_region.pl \
    data-raw/output/family.cnv \
    data-raw/refs/refGene.txt \
    --refgene --reflink data-raw/refs/refLink.txt \
    > data-raw/output/family.gene.cnv \
    2> data-raw/output/family.gene.log

# copy files to package
mkdir -p inst/extdata
cp data-raw/PennCNV/lib/hhall.hmm inst/extdata
cp data-raw/PennCNV/example/example.pfb inst/extdata
cp data-raw/PennCNV/example/*.txt inst/extdata
cp data-raw/output/*.cnv inst/extdata
cp data-raw/output/*.log inst/extdata

cat inst/extdata/offspring.txt | cut -f 1,3,4,5,6,7 > inst/extdata/offspring.nopos.txt
