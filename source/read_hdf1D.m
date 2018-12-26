function [X,Pressure,Density,V,time]=read_hdf1D(filename)

% Read the HDF5 data
Density=h5read(filename,'/hydrodynamic/density');
Pressure=h5read(filename,'/hydrodynamic/pressure');
X=h5read(filename,'/geometry/edges');
V=h5read(filename,'/hydrodynamic/velocity');
time=h5read(filename,'/time');


