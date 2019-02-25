# project_NB_SE
Neuroblastoma SE landscape project
snakemake -k -w 50 --jobs 50 --use-conda --cluster-config cluster/bsub.yaml --cluster "bsub -W {cluster.walltime} -n {cluster.cores} -M {cluster.memory} -J {cluster.name} -o {cluster.output} -e {cluster.error} {cluster.queue}" 
