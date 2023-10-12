# SleeperPests

Repository contains 2 types of pests spread and management simulation function:
1) serialised version;
2) parallelised version.

   * INApest -- serialised function which models the spread of the pests without tracking local population dynamics;
   * INApestMeta -- this function tracks the spread which is driven by the local population dynamics; it is possible for the user to define timestep function to track local population dynamics;

   * INApestParallel -- parallelised function which models the spread of the pests without tracking local population dynamics;
   * INApestMetaParallel -- parallelised function which tracks the spread driven by the local population dynamics; note that the timestep function is hardcoded and not controlled by the user.
  
     The initial version of INApest() made calls to the function INAscene() from the INA package (available on GitHub but not on CRAN at the moment; see https://github.com/GarrettLab/INA). The setup for the required functions is in ParallelSetup.r.

     The most recent versions of INApest() manually implement the steps that in the previous versions were done by calling INAscene().

Within INApestMeta -- there are 2 versions: parallel and serialised. There is an additional version of INApestMeta tracking of populations across multiple landuses whithin a location. This can be used for large scale modelling where locations span multiple properties, or multiple landuse types and where management of pest populations differs between landuse types.
     

## Examples
  
The example implementation of INApest() is Chilean needle grass. 
Tomato red spider mite is an example implementation of INApestMeta().
   
