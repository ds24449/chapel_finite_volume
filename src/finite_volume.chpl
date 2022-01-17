use IO;
use DataArray;
use FDSolver;
use linspace;

// proc getConserved(rho, vx, vy, P, gamma:real, vol:real ){
// 	"""
//     Calculate the conserved variable from the primitive
// 	rho      is matrix of cell densities
// 	vx       is matrix of cell x-velocity
// 	vy       is matrix of cell y-velocity
// 	P        is matrix of cell pressures
// 	gamma    is ideal gas gamma
// 	vol      is cell volume
// 	Mass     is matrix of mass in cells
// 	Momx     is matrix of x-momentum in cells
// 	Momy     is matrix of y-momentum in cells
// 	Energy   is matrix of energy in cells
// 	"""
// 	var Mass   = rho * vol;
// 	var Momx   = rho * vx * vol;
// 	var Momy   = rho * vy * vol;
// 	var Energy = (P/(gamma-1) + 0.5*rho*(vx**2+vy**2))*vol;
	
// 	return {Mass, Momx, Momy, Energy};
// }


// proc getPrimitive( Mass, Momx, Momy, Energy, gamma:real, vol:real ){
// 	"""
//     Calculate the primitive variable from the conservative
// 	Mass     is matrix of mass in cells
// 	Momx     is matrix of x-momentum in cells
// 	Momy     is matrix of y-momentum in cells
// 	Energy   is matrix of energy in cells
// 	gamma    is ideal gas gamma
// 	vol      is cell volume
// 	rho      is matrix of cell densities
// 	vx       is matrix of cell x-velocity
// 	vy       is matrix of cell y-velocity
// 	P        is matrix of cell pressures
// 	"""
// 	var rho = Mass / vol;
// 	var vx  = Momx / rho / vol;
// 	var vy  = Momy / rho / vol;
// 	var P   = (Energy/vol - 0.5*rho * (vx**2+vy**2)) * (gamma-1);
	
// 	return {rho, vx, vy, P};
// }

proc getGradient( f, dx){
    // """
    // Calculate the gradients of a field
    // f        is a matrix of the field
    // dx       is the cell size
    // f_dx     is a matrix of derivative of f in the x-direction
    // f_dy     is a matrix of derivative of f in the y-direction
    // """
    var Solver = new owned FDSolver(f);
    Solver.apply_bc(ax = 0, lt_bc = "periodic", rt_bc = "periodic", accuracy = 2);
    var f_dx = Solver.Finite_Difference( scheme = "central", order = 1, accuracy = 2, step = dx, axis = 0);
    var f_dy = Solver.Finite_Difference( scheme = "central", order = 1, accuracy = 2, step = dx, axis = 1);
    return {f_dx, f_dy};
}

proc sincosTest(){
	var saveFile = open("Tests/Data/sincos2DAvgError.txt",iomode.cw);
	var saveFileWriter = saveFile.writer();

	var start = 100;
	var end = 1000;
	var step = start;

	var errors: [1..end/start] real;
	for n in start..end by step{
    	var sincos = new owned DataArray(eltType = real,size = {1..n,1..n},dimensions = {"Y","X"}); // An StenArray object with dimension (n-1x1)
    	var trueValue = new owned DataArray(eltType = real,size = {1..n,1..n},dimensions = {"Y","X"});

    	var grid:[1..n] real = linspace(0,2*pi,n,false);
    	var h = grid[2]-grid[1];

    	forall i in 1..n do {
        	for j in 1..n do{
            	sincos.arr[i,j] = sin(grid[j])*cos(grid[i]);//Sin(x)*Cos(y)
        	}
    	}
    
    	forall i in 1..n do {
        	for j in 1..n do{
            	trueValue.arr[i,j] = -sin(grid[i])*sin(grid[j]) + cos(grid[i])*cos(grid[j]); 
        	}
    	}

    	var Solver = new owned FDSolver(sincos);
    	Solver.apply_bc(["X" => "periodic", "Y" => "periodic"], 2);
    	var result_dx = Solver.Finite_Difference(scheme="central",order=1,accuracy=2,step=h,axis=0);
		var result_dy = Solver.Finite_Difference(scheme="central",order=1,accuracy=2,step=h,axis=1); // This is giving wrong answers
		var result = (result_dx + result_dy):trueValue.type;

   		var avgError:real = 0.0;
    	for i in trueValue.dom{
        	avgError += abs(result.arr[i]-trueValue.arr[i]);
    	}
    	
		avgError /= n*n;
    	errors[n/start] = avgError;

    	if(n/end == 1) then saveFileWriter.writeln(h);
    	else saveFileWriter.write(h," ");
	}
	
	saveFileWriter.writeln(errors);
	saveFileWriter.close();
	saveFile.fsync();
	saveFile.close();
}


sincosTest();