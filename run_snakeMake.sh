#!/usr/bin/sh
snakemake -s ./RNAseq_snakeFile -j 500 \
--cluster-config ./cluster.json \
--latency-wait 30 \
--cluster 'sbatch --mem={cluster.mem} -t {cluster.time} -c {threads}'
