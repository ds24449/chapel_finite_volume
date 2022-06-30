use IO;
use DataArray;
use FDSolver;
use linspace;
use Math;
use ntCDF;

class FiniteVolume {
	const TAG = FiniteVolume : string;

	// Variables
	var rho, vx, vy, P;
	var dx, vol: real;

	// Simulation Parameters
	var resolution: int;
	var gamma, courantFactor, t, tOut, tEnd: real;
	var useSlopeLimiting, saveData: bool;

	proc init (rho, vx, vy, P, dx, vol, resolution, gamma, courantFactor, t, tOut, tEnd, useSlopeLimiting = false, saveData = true) {
		this.rho = rho;
		this.vx = vx;
		this.vy = vy;
		this.P = P;

		this.dx = dx;
		this.vol = vol;

		this.resolution = resolution;
		this.gamma = gamma;
		this.courantFactor = courantFactor;
		this.t = t;
		this.tOut = tOut;
		this.tEnd = tEnd;
		this.useSlopeLimiting = useSlopeLimiting;
		this.saveData = saveData;
	}


	proc getConserved(){
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

	proc getPrimitive( Mass, Momx, Momy, Energy){
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

	proc getGradient(f){
		/*
		Calculate the gradients of a field
		f        is a matrix of the field
		dx       is the cell size
		f_dx     is a matrix of derivative of f in the x-direction
		f_dy     is a matrix of derivative of f in the y-direction

		Notes: For discontinuities Slope Limiters are used. Discontinuity occurs for Supersonic Speeds
		*/
		// var f_DataArray = new DataArray(f,dimensions = {"Y","X"}); // IF not DataArray
		// var Solver = new owned FDSolver(f_DataArray);
		// Solver.apply_bc(["X" => "no_bc", "Y" => "no_bc"], 2);
		// var f_dx = Solver.Finite_Difference( scheme = "central", order = 1, accuracy = 2, step = dx, axis = 0);
		// var f_dy = Solver.Finite_Difference( scheme = "central", order = 1, accuracy = 2, step = dx, axis = 1);
		// return (f_dx.arr,f_dy.arr);
		var f_dx = (roll(f, -1, axis = 0) - roll(f, 1, axis = 0)) / (2*dx);
		var f_dy = (roll(f, -1, axis = 1) - roll(f, 1, axis = 1)) / (2*dx);
		return (f_dx, f_dy);
	}

	proc extrapolateInSpaceToFace(f, f_dx, f_dy){
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
		// directions for np_roll() 
		var R = -1;   // right
		var L = 1;    // left
	
		var f_XL = f - (f_dx * dx/2);
		f_XL = roll(f_XL,R,axis=0);

		var f_XR = f + (f_dx * dx/2);
	
		var f_YL = f - (f_dy * dx/2);
		f_YL = roll(f_YL,R,axis=1);

		var f_YR = f + (f_dy * dx/2);
	
		return (f_XL, f_XR, f_YL, f_YR);
	}

	proc getFlux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R){
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
		var C = np_maximum( C_L, C_R );
	
		// add stabilizing diffusive term
		flux_Mass   -= C * 0.5 * (rho_L - rho_R);
		flux_Momx   -= C * 0.5 * (rho_L * vx_L - rho_R * vx_R);
		flux_Momy   -= C * 0.5 * (rho_L * vy_L - rho_R * vy_R);
		flux_Energy -= C * 0.5 * ( en_L - en_R );

		return (flux_Mass, flux_Momx, flux_Momy, flux_Energy);
	}

	proc applyFluxes(F, flux_F_X, flux_F_Y, dt){
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

	proc slopeLimit(f, f_dx, f_dy){
		/*
		Apply slope limiter to slopes
			f        is a matrix of the field
			dx       is the cell size
			f_dx     is a matrix of derivative of f in the x-direction
			f_dy     is a matrix of derivative of f in the y-direction
		*/
		// directions for np.roll()
		var R = -1;   // right
		var L =  1;    // left
		// var N = f_dx.domain.dim(0).high;

		//TODO: Remove like_ones, change minimum and maximum function to work with int and real
		f_dx = np_maximum(like_zeros(f.domain), np_minimum(like_ones(f.domain), ( (f-roll(f,L,axis=0))/dx)/(f_dx + 1.0e-8*(f_dx==0):real))) * f_dx;
		f_dx = np_maximum(like_zeros(f.domain), np_minimum(like_ones(f.domain), (-(f-roll(f,R,axis=0))/dx)/(f_dx + 1.0e-8*(f_dx==0):real))) * f_dx;
		f_dy = np_maximum(like_zeros(f.domain), np_minimum(like_ones(f.domain), ( (f-roll(f,L,axis=1))/dx)/(f_dy + 1.0e-8*(f_dy==0):real))) * f_dy;
		f_dy = np_maximum(like_zeros(f.domain), np_minimum(like_ones(f.domain), (-(f-roll(f,R,axis=1))/dx)/(f_dy + 1.0e-8*(f_dy==0):real))) * f_dy;
		return (f_dx, f_dy);
	}

	proc run(){

		// Get conserved variables
		var (Mass, Momx, Momy, Energy) = getConserved();
		var outputCount = 1;

		// Simulation Main Loop
		while t < tEnd{

			LOG(TAG , t:string);
			// get Primitive variables
			(rho, vx, vy, P) = getPrimitive(Mass, Momx, Momy, Energy);

			var temp_arr = dx / (sqrt(gamma * P/rho) + sqrt(vx**2 + vy**2)); 
			var new_temp = min_of_arr( temp_arr );
			var dt = courantFactor*new_temp;
			var saveThisTurn = false;
			if (t + dt > outputCount*tOut){
				dt = outputCount * tOut - t;
				saveThisTurn = true;
			}
			// calculate gradients
			var (rho_dx, rho_dy) = getGradient(rho);
			var (vx_dx,  vx_dy) = getGradient(vx);
			var (vy_dx,  vy_dy) = getGradient(vy);
			var (P_dx,   P_dy) = getGradient(P);


			// slope limit gradients
			if useSlopeLimiting{
				(rho_dx, rho_dy) = slopeLimit(rho, rho_dx, rho_dy);
				(vx_dx,  vx_dy) = slopeLimit(vx, vx_dx,  vx_dy);
				(vy_dx,  vy_dy) = slopeLimit(vy, vy_dx,  vy_dy);
				(P_dx,   P_dy) = slopeLimit(P, P_dx,   P_dy);
			}
			// extrapolate half-step in time
			var rho_prime = rho - 0.5*dt * (vx * rho_dx + rho * vx_dx + vy * rho_dy + rho * vy_dy);
			var vx_prime = vx - 0.5*dt * (vx * vx_dx + vy * vx_dy + (1/rho) * P_dx);
			var vy_prime = vy - 0.5*dt * (vx * vy_dx + vy * vy_dy + (1/rho) * P_dy);
			var P_prime = P - 0.5*dt * (gamma * P * (vx_dx + vy_dy) + vx * P_dx + vy * P_dy);


			// extrapolate in space to face centers
			var (rho_XL, rho_XR, rho_YL, rho_YR) = extrapolateInSpaceToFace(rho_prime, rho_dx, rho_dy);
			var (vx_XL,  vx_XR,  vx_YL,  vx_YR) = extrapolateInSpaceToFace(vx_prime,  vx_dx,  vx_dy);
			var (vy_XL,  vy_XR,  vy_YL,  vy_YR) = extrapolateInSpaceToFace(vy_prime,  vy_dx,  vy_dy);
			var (P_XL,   P_XR,   P_YL,   P_YR) = extrapolateInSpaceToFace(P_prime,   P_dx,   P_dy);

			// compute fluxes (local Lax-Friedrichs/Rusanov)
			var (flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X) = getFlux(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR);
			var (flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y) = getFlux(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR);

			// update solution
			Mass = applyFluxes(Mass, flux_Mass_X, flux_Mass_Y, dt);
			Momx = applyFluxes(Momx, flux_Momx_X, flux_Momx_Y, dt);
			Momy = applyFluxes(Momy, flux_Momy_X, flux_Momy_Y, dt);
			Energy = applyFluxes(Energy, flux_Energy_X, flux_Energy_Y, dt);

			// update time
			t += dt;

			if ((saveData && saveThisTurn) || (t >= tEnd)){
				LOG(TAG, outputCount:string);
				write2DArray(rho, resolution, resolution, fileName = "rho_" + outputCount:string + ".nc");
				// write2DArray(Mass, N, N, fileName = "Mass_" + outputCount:string + ".nc");
				// write2DArray(Momx, N, N, fileName = "Momx_" + outputCount:string + ".nc");
				// write2DArray(Momy, N, N, fileName = "Momy_" + outputCount:string + ".nc");
				// write2DArray(Energy, N, N, fileName = "Energy_" + outputCount:string + ".nc");

				outputCount += 1;
			}
		}
		return (rho, vx, vy, P);
		LOG(TAG, "Program Ended");
	}
}

proc main() {
	// Simulation parameters
	var N                      = 128; // resolution
	var boxsize                = 1.0;
	var gamma		           = 5.0/3.0; // ideal gas gamma
	var courant_fac            = 0.4;
	var t                      = 0.0;
	var tEnd                   = 2.0;
	var tOut                   = 0.02; // draw frequency
	var useSlopeLimiting       = false;
	var saveData 		       = true; // switch on for plotting as the simulation goes along

	/*FINITE VOLUME*/


	// MESH
	var dx:real = boxsize / N;
	var vol = dx**2;
	var xlin = linspace(0.5*dx, boxsize-0.5*dx, N);
	var Y:[0..N-1,0..N-1] real;
	var X:[0..N-1,0..N-1] real;

	forall i in 0..N-1{
		forall j in 0..N-1{
			X[i, j] = xlin[i];
		}
	}

	forall i in 0..N-1{
		forall j in 0..N-1{
			Y[i, j] = xlin[j];
		}
	}

	// Generate Initial Conditions - opposite moving streams with perturbation
	var w0 = 0.1;
	var sigma = 0.05/sqrt(2.0);
	var rho = 1.0 + (abs(Y-0.5) < 0.25):real;
	var vx = -0.5 + (abs(Y-0.5) < 0.25):real;
	var vy = w0*sin(4*pi*X) * (exp(-(Y-0.25)**2 / (2 * sigma**2)) + exp(-(Y-0.75)**2/(2*sigma**2)));
	var P = 2.5 * like_ones(X.domain);

	//in rho, in vx, in vy, in P, in vol, in dx, in resolution, in gamma, in courantFactor, in t, in tOut, in tEnd 

	var finiteVol = new FiniteVolume(rho, vx, vy, P, dx, vol, N, gamma, courant_fac, t, tOut, tEnd, useSlopeLimiting, saveData);
	finiteVol.run();
}

main();