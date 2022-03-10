// To install netCDF module: sudo apt-get install libnetcdf-dev libnetcdff-dev
// chpl --fast e/GT_chpl/netcdfWrite.chpl -I/usr/include -L/usr/lib/x86_64-linux-gnu
// Must have C_NetCDF library in the system, and use -I flag to provide path for header file
// and -L to provide path for library
// reference: https://westgrid.github.io/trainingMaterials/materials/cFromChapel20200318.pdf

use NetCDF.C_NetCDF;

proc cdfError(e:c_int) {
    if e != NC_NOERR {
        writeln("Error: ", nc_strerror(e): string);
        exit(2);
    }
}

proc write2DArray(T:[?D],nx:int,ny:int,fileName:string,dataName:c_string = "data"){
    // var dataName = name:c_string;
    var loc_fileName = fileName.c_str();
    var ncid, xDimID, yDimID, varID: c_int;
    var dimIDs: [0..1] c_int;

    var W:[D] c_float;
    forall i in D do{
        W[i] = T[i]: c_float;
    }

    cdfError(nc_create(loc_fileName, NC_NETCDF4, ncid));
    cdfError(nc_def_dim(ncid, "x", nx: size_t, xDimID));
    cdfError(nc_def_dim(ncid, "y", ny: size_t, yDimID));

    dimIDs = [xDimID, yDimID];

    cdfError(nc_def_var(ncid,dataName,NC_FLOAT,2,dimIDs[0],varID));
    cdfError(nc_def_var_deflate(ncid, varID, NC_SHUFFLE, deflate=1, deflate_level=9));
    cdfError(nc_enddef(ncid));
    cdfError(nc_put_var_float(ncid,varID,W[D.low]));
    cdfError(nc_close(ncid));

    return 0;
}

proc WriteTest(){
    var X: [1..10,1..10] real = 1;
    writeln("Exiting with code: " + write2DArray(X,10,10,"TestFile","Data"):string);
}

// WriteTest();

