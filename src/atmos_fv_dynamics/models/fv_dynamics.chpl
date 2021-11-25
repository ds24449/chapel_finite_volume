// We are gonna use LIMA version
module fv_dynamics_mod{
    proc fv_dynamics(in im:int, in jm:int, in km:int, in jfirst:int, in jlast:int, 
     in nq:int, in pnats:int, in ndt:int, in p_map:bool, in consv:real, inout u:[?dom_u] real, 
     inout v:[?dom_v] real, inout delp:real, inout pt:[?dom_pt] real, inout q:[?dom_q] real, 
     inout ps:real, inout pe:[?dom_pe] real, inout pk:[?dom_pk] real, inout pkz:[?dom_pkz] real, 
     inout phis:real, inout omga:[?dom_omga] real,inout peln:[dom_pe] real, in ptop:real, in om:real, 
     in r_vir:real, in cp:real, in rg:real, in cappa:real, in ae:real,inout ua:[?dom_a] real,
     inout va:[dom_a] real, in Time_next: time_type){
        /*
            INPUT PARAMETERS:
            Intent in:
            Integer
                im      : dimension in east-west
                jm      : dimension in North-South
                km      : number of Lagrangian layers
                jfirst  : starting latitude index for MPI
                jlast   : ending latitude index for MPI
                nq      : total # of tracers to be advected
                pnats   : number of non-advected tracers
                ndt     : the large time step in seconds

            Real
                ptop    : constant pressure at model top (pascal)          
                om      : angular velocity of earth's rotation  
                r_vir   : constant for the virtual effect
                cp      : heat capacity of air at constant pressure
                rg
                cappa   : rg/cp
                ae      : radius of the earth (m)
                consv   : energy conserving correlation factor: 0 = correction, 1 = full correction

            Bool 
                p_map   : partial remaping 
            
            time_type
                Time_next
            
            Intent inout
            Real
                ps      : surface pressure (pascal)   [NOTE] dom = (im,jfirst:jlast)
                delp    : pressure thickness (pascal) [NOTE] dom = (im,jfirst:jlast,km)

            Real, ghosted prog arrays
                u   : u-wind (m/s) [NOTE] dom = (im,jfirst-ng_d:jlast+ng_s,km)
                v   : v-wind (m/s) [NOTE] dom = (im,jfirst-ng_d:jlast+ng_d,km)
                pt  : temperature (K) [NOTE] dom = (im,jfirst-ng_d:jlast+ng_d,km)
                q   : [NOTE] dom = (im,jfirst-ng_d:jlast+ng_d,km,nq)
                ! tracers (e.g., specific humidity)
                ! tracer mass / moist_air_mass
                
                pkz: finite-volume mean of pk [NOTE] dom = (im,jfirst:jlast,km)

            [NOTE] peln are not needed as input if consv = 0.
                pe      :pressure (pascal) at layer edges [NOTE] dom = (im,km+1,jfirst:jlast)
                peln    :log pressure (pe) at layer edges [NOTE] dom = (im,km+1,jfirst:jlast)
                ua      :A grid u-wind (m/s) [NOTE] dom = (im,jfirst:jlast,km)
                va      :A grid v-wind (m/s) [NOTE] dom = (im,jfirst:jlast,km)

            OUTPUT (input values are not used):
                pk :pe**cappa (cappa = rg/cp) [NOTE] dom = (im,jfirst:jlast,km+1)

            [NOTE] outputed values can be used by physdrv
            omga :vertical pressure velocity (pa/sec) [NOTE] dom = (im,jfirst:jlast,km)
            This is the rate of change of the Lagrangian interface in pascal/sec unit.
        */


        // LOCAL VARIABLES
        var i, j, k, it, padv, iq : int; 
        var umax : real = 300.0; //estimated upper bound of the maximum u-wind (m/s)
        var dt, bdt, kappa, cp_dyn, frac : real;
        // frac: tracer time split fraction

        // Local auto arrays
        var  wk1: [im+2, jfirst:jlast+2] real;
        var  wk2: [im+2, jfirst:jlast+2] real;
        var  tte: [jfirst:jlast] real;

        // temporary work arrays
        var cx: [im,jfirst-ng_d:jlast+ng_d,km] real;
        var mfx: [im,jfirst:jlast,km] real;
        var cy: [im,jfirst:     jlast+1,   km] real;
        var mfy: [im,jfirst:     jlast+1,   km] real;
        var dpt: [im,jfirst-1:   jlast+1,   km] real;
        var vc: [im,jfirst-2:   jlast+2,   km] real;
        var delpf: [im,jfirst-ng_d:jlast+ng_d,km] real;
        var uc: [im,jfirst-ng_d:jlast+ng_d,km] real;
        var dwz: [im,jfirst-1:   jlast,     km+1] real;
        var pkc: [im,jfirst-1:   jlast+1,   km+1] real; 
        var wz: [im,jfirst-1:   jlast+1,   km+1] real;
        // Save a copy of pe for the computation of vertical velocity

        var pem: [im,km+1,jfirst:jlast] real;

        // [DEV-NOTE] ignoring use_shared_pointer pragma
        
        var js2g0,jn2g0: int;

        // #ifdef SW_DYN
        // ! The following code segment is for shallow water test cases
        //     real yy, tday, aoft, pi, gv, ssec
        //     data ssec /0./

        //     kappa  = 1.
        //     cp_dyn = 1.
        // #else
        //     kappa  = cappa
        //     cp_dyn = cp
        // #endif

        js2g0 = max(2,jfirst)
        jn2g0 = min(jm-1,jlast)
        padv = n1-pnats // LOCAL VARIABLE 
        

    }

    writeln("Hello World?");
}