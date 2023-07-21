#!/bin/bash
rm singularity/*.simg
rm ref/intervals/*.bed
rm ref/intervals/*.intervals
rm ref/intervals/itv.17.ok
rm ref/refFlat.txt
rm ref/af-only-gnomad.hg38.vcf.gz*
tar zcvf logs.tgz logs
rm -rf logs
