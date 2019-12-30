# MPI_Mandelbrot

This program is designed to be run on a cluster of computers although could be run on a single multi-threaded core.
If running on a cluster, ensure that all computers have access to a shared file system (e.g via NFS) where this code should live.
All of the clients should also have the appropriate MPI package installed and the master should be able to login to each client via public key SSH (or at least without password).
The file `hostfile` should contain a list of all hosts that will be running in the cluster including `localhost` (the master).

Once these requirements are met the code can be compiled with `make`.
The code can then be run by

```
mprun -np <n_cores> --hostfile hostfile ./mandelbrot.exe
```

where `<n_cores>` should be replaced by the number of available cores on the cluster.

`mandelbrot-serial.exe` can be run as normal in order to compare multi-threaded perforamnce to serial performance.

## Direction

Over time I intend to add the features I previously developed in [Multithreaded_Mandelbrot](https://github.com/tomchaplin/Multithreaded_Mandelbrot) to this MPI implementation, namely:

* Animation
* Zooming
* A variety of colouring algorithms
