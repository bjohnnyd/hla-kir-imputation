#!/usr/bin/env bash
#properties = {properties}
#$ -lwalltime={cluster.time}
#$ -lselect={cluster.nodes}:ncpus={cluster.cpus}:mem={cluster.mem}
#$ -N {cluster.name}
#$ -e {cluster.logs}
#$ -o {cluster.logs}
{exec_job}
