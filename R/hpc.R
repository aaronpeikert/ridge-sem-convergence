library("future.batchtools")
# ssh keys login
#ssh-copy-id -i /home/aaron/.ssh/id_rsa peikert@tardis.mpib-berlin.mpg.de
tardis <- parallelly::makeClusterPSOCK("tardis.mpib-berlin.mpg.de",
                                       port='random', user="peikert",
                                       rscript=c("/opt/software/R/4.1.0/bin/Rscript"), 
                                       homogeneous = TRUE)
ncpus <- 8
plan(list(tweak(cluster, workers=tardis),
          tweak(batchtools_slurm,
                workers = 100,
                template = "/home/mpib/peikert/sem-reg/.batchtools.slurm.singularity.tmpl",
                resources=list(ncpus=ncpus,
                               memory='700m',
                               walltime=172800,
                               partition=c('long'))),
          sequential))
#plan(transparent)
session <- function()with(sessionInfo(), list(rversion = paste0(R.version$major, ".", R.version$minor),
                                              os = running))
# list(local = session(),
#     login_node = value(future(session())),
#     worker = value(future(value(future(session())))))


