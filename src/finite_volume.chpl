use IO;
use DataArray;
use FDSolver;
use linspace;
use Math;

proc getConserved(rho, vx, vy, P, gamma:real, vol:real ){
	/*
    Calculate the conserved variable from the primitive
	Input:
	rho      is matrix of cell densities
	vx       is matrix of cell x-velocity
	vy       is matrix of cell y-velocity
	P        is matrix of cell pressures
	gamma    is ideal gas gamma
	vol      is cell volume

	Output: 
	Mass     is matrix of mass in cells
	Momx     is matrix of x-momentum in cells
	Momy     is matrix of y-momentum in cells
	Energy   is matrix of energy in cells
	*/
	var Mass   = rho * vol;
	var Momx   = rho * vx * vol;
	var Momy   = rho * vy * vol;
	var Energy = (P/(gamma-1) + 0.5*rho*(vx**2+vy**2))*vol;
	
	return (Mass, Momx, Momy, Energy);
}


proc getPrimitive( Mass, Momx, Momy, Energy, gamma:real, vol:real ){
	/*
    Calculate the primitive variable from the conservative
	INPUT: 
	Mass     is matrix of mass in cells
	Momx     is matrix of x-momentum in cells
	Momy     is matrix of y-momentum in cells
	Energy   is matrix of energy in cells
	gamma    is ideal gas gamma
	vol      is cell volume

	OUTPUT:
	rho      is matrix of cell densities
	vx       is matrix of cell x-velocity
	vy       is matrix of cell y-velocity
	P        is matrix of cell pressures
	*/
	var rho = Mass / vol;
	var vx  = Momx / rho / vol;
	var vy  = Momy / rho / vol;
	var P   = (Energy/vol - 0.5*rho * (vx**2+vy**2)) * (gamma-1);
	
	return (rho, vx, vy, P);
}

proc getGradient( f, dx:real){
    /*
    Calculate the gradients of a field
    f        is a matrix of the field
    dx       is the cell size
    f_dx     is a matrix of derivative of f in the x-direction
    f_dy     is a matrix of derivative of f in the y-direction

	Notes: For discontinuities Slope Limiters are used. Discontinuity occurs for Supersonic Speeds
    */
	// var f_DataArray = new DataArray(f,dimensions = {"Y","X"}); // IF not DataArray
    var Solver = new owned FDSolver(f);
    Solver.apply_bc(["X" => "periodic", "Y" => "periodic"], 2);
    var f_dx = Solver.Finite_Difference( scheme = "central", order = 1, accuracy = 2, step = dx, axis = 0);
    var f_dy = Solver.Finite_Difference( scheme = "central", order = 1, accuracy = 2, step = dx, axis = 1);
    return (f_dx.arr,f_dy.arr);
}

proc extrapolateInSpaceToFace(f: DataArray, f_dx: DataArray, f_dy: DataArray, dx: real){
  	/*
  	Calculate the gradients of a field
  	f        is a matrix of the field
  	f_dx     is a matrix of the field x-derivatives
  	f_dy     is a matrix of the field y-derivatives
  	dx       is the cell size
  	f_XL     is a matrix of spatial-extrapolated values on `left' face along x-axis 
  	f_XR     is a matrix of spatial-extrapolated values on `right' face along x-axis 
  	f_YR     is a matrix of spatial-extrapolated values on `left' face along y-axis 
  	f_YR     is a matrix of spatial-extrapolated values on `right' face along y-axis 
  	*/
  	// directions for np.roll() 
  	var R = -1;   // right
  	var L = 1;    // left
  
  	var f_XL = (f - f_dx):f.type;
  	f_XL.arr *= dx/2; 				// NOTE: Considering this is DataArray
	f_XL = roll(f_XL,R,axis=0);

  	var f_XR = (f + f_dx):f.type;
	f_XR.arr *= dx/2;
  
  	var f_YL = (f - f_dy):f.type;
  	f_YL.arr *= dx/2;
  	f_YL = roll(f_YL,R,axis=1);

  	var f_YR = (f + f_dy):f.type;
	f_YR.arr *= dx/2;
  
  	return (f_XL, f_XR, f_YL, f_YR);
}

proc getFlux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R, gamma:real){
  	/*
  	Calculate fluxed between 2 states with local Lax-Friedrichs/Rusanov rule 
	Input:
		rho_L        is a matrix of left-state  density
		rho_R        is a matrix of right-state density
		vx_L         is a matrix of left-state  x-velocity
		vx_R         is a matrix of right-state x-velocity
		vy_L         is a matrix of left-state  y-velocity
		vy_R         is a matrix of right-state y-velocity
		P_L          is a matrix of left-state  pressure
		P_R          is a matrix of right-state pressure

		gamma        is the ideal gas gamma
	
	OUTPUT:
		flux_Mass    is the matrix of mass fluxes
		flux_Momx    is the matrix of x-momentum fluxes
		flux_Momy    is the matrix of y-momentum fluxes
		flux_Energy  is the matrix of energy fluxes
  	*/
  
  	// left and right energies
  	var en_L = P_L/(gamma-1)+0.5*rho_L * (vx_L**2+vy_L**2);
  	var en_R = P_R/(gamma-1)+0.5*rho_R * (vx_R**2+vy_R**2);

  	// compute star (averaged) states
  	var rho_star  = 0.5*(rho_L + rho_R);
  	var momx_star = 0.5*(rho_L * vx_L + rho_R * vx_R);
  	var momy_star = 0.5*(rho_L * vy_L + rho_R * vy_R);
  	var en_star   = 0.5*(en_L + en_R);
  
  	var P_star = (gamma-1)*(en_star-0.5*(momx_star**2+momy_star**2)/rho_star);
  
  	// compute fluxes (local Lax-Friedrichs/Rusanov)
  	var flux_Mass   = momx_star;
  	var flux_Momx   = momx_star**2/rho_star + P_star;
  	var flux_Momy   = momx_star * momy_star/rho_star;
  	var flux_Energy = (en_star+P_star) * momx_star/rho_star;
  
  	// find wavespeeds
  	var C_L = sqrt(gamma*P_L/rho_L) + abs(vx_L);
  	var C_R = sqrt(gamma*P_R/rho_R) + abs(vx_R);
  	var C = max( C_L, C_R );
  
  	// add stabilizing diffusive term
  	flux_Mass   -= C * 0.5 * (rho_L - rho_R);
  	flux_Momx   -= C * 0.5 * (rho_L * vx_L - rho_R * vx_R);
  	flux_Momy   -= C * 0.5 * (rho_L * vy_L - rho_R * vy_R);
  	flux_Energy -= C * 0.5 * ( en_L - en_R );

  	return (flux_Mass, flux_Momx, flux_Momy, flux_Energy);
}

proc applyFluxes(F, flux_F_X, flux_F_Y, dx, dt){
  	/*
  	Apply fluxes to conserved variables
  		F        is a matrix of the conserved variable field
  		flux_F_X is a matrix of the x-dir fluxes
  		flux_F_Y is a matrix of the y-dir fluxes
  		dx       is the cell size
  		dt       is the timestep
  	*/
  	// directions for np.roll() 
  	var R = -1;   // right
  	var L = 1;    // left
  
  	// update solution
  	F += - dt * dx * flux_F_X;
	F +=   dt * dx * roll(flux_F_X,L,axis=0); 
  	F += - dt * dx * flux_F_Y;
	F +=   dt * dx * roll(flux_F_Y,L,axis=1);
  
  	return F;
}



proc main( Mass, Momx, Momy, Energy, gamma:real, vol:real ) {
	// simulation parameters
	var t = 0;
	var gamma = 5/3;
	var courant_fac = 0.4;
	var tEnd = 2;

	// simulation loop
	while(t < tEnd) {
		var rho, vx, vy, P: real;
		// fetch the primitive variables
		(rho, vx, vy, P) = getPrimitive(Mass, Momx, Momy, Energy, gamma, vol);

		// calculate the time-step
		var dt = courant_fac * np_min(dx / sqrt(gamma*p/rho) + sqrt(vx**2 + vy**2));
		/* TODO : create the np_min function */

		// calculate the gradients
		var rho_dx, rho_dy, vx_dx, vx_dy, vy_dx, vy_dy, P_dx, P_dy: real;
		(rho_dx, rho_dy) = getGradient(rho, dx);
		(vx_dx, vx_dy) = getGradient(vx, dx);
		(vy_dx, vy_dy) = getGradient(vy, dx);
		(P_dx, P_dy) = getGradient(P, dx);

		// extrapolate half-step in time
		var rho_prime = rho - 0.5*dt * ( vx * rho_dx + rho * vx_dx + vy * rho_dy + rho * vy_dy);
		var vx_prime = vx - 0.5*dt * ( vx * vx_dx + vy * vx_dy + (1/rho) * P_dx );
		var vy_prime = vy - 0.5*dt * ( vx * vy_dx + vy * vy_dy + (1/rho) * P_dy );
		var P_prime = P - 0.5*dt * ( gamma*P * (vx_dx + vy_dy)  + vx * P_dx + vy * P_dy );

		// extrapolate in space to face centers
		var rho_XL, rho_XR, rho_YL, rho_YR, vx_XL, vx_XR, vx_YL, vx_YR, vy_XL, vy_XR, vy_YL, vy_YR, P_XL, P_XR, P_YL, P_YR:real;
		(rho_XL, rho_XR, rho_YL, rho_YR) = extrapolateInSpaceToFace(rho_prime, rho_dx, rho_dy, dx);
		(vx_XL, vx_XR, vx_YL, vx_YR) = extrapolateInSpaceToFace(vx_prime, vx_dx, vx_dy, dx);
		(vy_XL, vy_XR, vy_YL, vy_YR) = extrapolateInSpaceToFace(vy_prime, vy_dx, vy_dy, dx);
		(P_XL, P_XR, P_YL, P_YR) = extrapolateInSpaceToFace(P_prime, P_dx, P_dy, dx);

		// compute fluxes (local Lax-Friedrichs/Rusanov)
		var flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X, flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y:real;
		(flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X) = getFlux(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, gamma);
		(flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y) = getFlux(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, gamma);
		
		// update solution
		var Mass = applyFluxes(Mass, flux_Mass_X, flux_Mass_Y, dx, dt);
		var Momx = applyFluxes(Momx, flux_Momx_X, flux_Momx_Y, dx, dt);
		var Momy = applyFluxes(Momy, flux_Momy_X, flux_Momy_Y, dx, dt);
		var Energy = applyFluxes(Energy, flux_Energy_X, flux_Energy_Y, dx, dt);
  
		// update time
		t += dt;
	}

}

/* TODO:
		- Assign datatypes in arguments

	STATUS: 
		- Compilation Errors: 0
		- Runtime Errors: 0
*/