# Hydro_fault_zone
2D hydraulic simulation in a poro-elastic layer (fault zone)

Mathematical formulation is given in Aochi et al. (2014)

Aochi, H., B. Poisson, R. Toussaint, X. Rachez and J. Schmittbuhl, Self-induced seismicity due to fluid circulation along faults, Geophys. J. Int., 196 (3), 1544-1563, doi:10.1093/gji/ggt356, 2014.

Compile
> gfortran ydro_distrib.f

Note: All the parameters are written in program.

RUN
> make hydro # output directory defined by "dir"

> ./a.out

OUTPUT (under ./hydro)
> ls hydro

> hoge1.dat hoge2.dat initial.dat snapshotXXXX.dat

hoge1.dat: given parameters

hoge2.dat: time marker of the simulation progress

initial.dat: snapshot for initial condition; i, j, kappa, phi, pore pressure, zone width, injection rate

snapshotXXXX.dat: snapshot at XXXX (hour); i j, pore pressure, phi, kappa, zone width



## Application and further development 

Douglas, J. and H. Aochi, Using estimated risk to develop exploitation strategies for Enhanced Geothermal Systems, Pageoph, 171, 1847-1858, doi: 10.1007/s00024-013-0765-8, 2014.

Aochi, H., T. Le Guenan, and A. Burnol, Developing subsurface energy exploitation strategies by considering seismic risk, Petrol. Geosci., 23, 298-305, doi:10.1144/petgeo2016-065, 2017.

Aochi, H. and J. Maury, Fault reactivation simulation due to fluidcirculation using multi-scale heterogeneity model, submitted to Pageoph, 2026.
