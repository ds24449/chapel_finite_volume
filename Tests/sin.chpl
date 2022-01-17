/*
    A simple test script checking derivative of sinx 
    chpl --module-dir src/ Tests/Spatial/sinAvgError.chpl && ./sinAvgError
*/
use IO; // To make a graph in python
use DataArray;
use FDSolver;
use linspace;

var saveFile = open("Tests/Data/sinAvgError.txt",iomode.cw);
var saveFileWriter = saveFile.writer();

config var start = 100;
config var end = 1000;
config var step = start;

var errors:[1..end/start] real;

for n in start..end by step{
    var sinArray = new owned DataArray(eltType = real,size = {1..n},dimensions = {"X"}); // An StenArray object with dimension (n-1x1)
    var cosArray = new owned DataArray(eltType = real,size = {1..n},dimensions = {"X"});

    var grid:[1..n] real = linspace(0,2*pi,n,false);

    var h = grid[2]-grid[1];

    forall i in sinArray.dom do {
        sinArray.arr[i] = sin(grid[i]);
    }
    
    forall i in cosArray.dom do {
        cosArray.arr[i] = cos(grid[i]);
    }

    var Solver = new owned FDSolver(sinArray);
    Solver.apply_bc(ax = 0, lt_bc = "periodic", rt_bc = "periodic", accuracy = 2);
    var result = Solver.Finite_Difference(scheme="central",order=1,accuracy=2,step=h,axis=0);
    // var result = sinArray;

    var avgError:real = 0.0;
    for i in result.dom{
        avgError += abs(result.arr[i]-cosArray.arr[i]);
    }
    avgError /= n;
    errors[n/start] = avgError;

    if(n/end == 1) then saveFileWriter.writeln(h);
    else saveFileWriter.write(h," ");
}
// writeln(errors);
saveFileWriter.writeln(errors);
saveFileWriter.close();
saveFile.fsync();
saveFile.close();
