######################################################
### Parallel setup for INApest
###
### The current parallel INApest functions use the base R `parallel`
### package with PSOCK clusters and parallel::parLapply(). `parallel`
### ships with R, so doParallel/foreach and abind are no longer required
### for the core parallel functions.
###
### INApestParallelInAScene additionally requires the INA package.
######################################################

### Optional dependency for INApestParallelInAScene
if (!requireNamespace("INA", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
  devtools::install_github("GarrettLab/INA")
}
