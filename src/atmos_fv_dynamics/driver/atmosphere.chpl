/* [DEV MESSAGES]: 
    - time_manager_mod   @DONE
        - time_type     : Type
        - get_time      : proc
        - set_time      : proc
        - + operator    : operator
        
    - fv_pack
        - fv_init       : proc
        - cold_start    : variable

    - fv_restart_mod
        - fv_restart    : proc

    - fv_diagnostics
        - fv_time       : var
        - fv_diag_init  : proc

    - hs_forcing_mod
        - hs_forcing_init   :proc
    

*/

// TODO: change Timi_inits calling like set_time and other stuff
class atmosphere_mod{
    private:
        var dt_atmos : real;
        var sec, days, seconds: int;

    proc init(in Time_init: time_type,in Time: time_type,in Time_step: time_type){
        var axes: [0..3] int;
        var ss, ds: int;

        // [DEV] WRITE VERSION AND NAMELIST - NOT USEFULL HERE IG? 

        // Compute physics/atmos time step in seconds
        get_time(Time_step,sec); //FIXME: call member method from class object
        dt_atmos = sec:real;

        // Initialize FV dynamical core
        fv_init( sec );
        fv_restart( days, seconds );

        if(!cold_start){
            // Check consistency in Time
            fv_time = set_time(seconds, days);
            get_time(Time, ss, ds); // FIXME: call member method from class object

            if(seconds != ss || days != ds){
                writeln("Fatal Error: Time Inconsistenct between fv_rst and INPUT/atmos_model.res");
            }else{
                fv_time = Time;
            }
        }

        fv_diag_init(axes, Time);

        if(nlev > 1) hs_forcing_init(axes,Time);
    }
}