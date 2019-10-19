#!/usr/bin/env bash
#properties = {rule}
#$ -lwalltime={properties["cluster"]["time"]}
#$ -lselect={cluster.nodes}:ncpus={cluster.cpus}:mem={cluster.mem}
#$ -N {cluster.name}
#$ -e {cluster.logs}
#$ -o {cluster.logs}
{exec_job}
