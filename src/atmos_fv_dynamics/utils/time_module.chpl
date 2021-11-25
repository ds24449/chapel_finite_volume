import Math.mod;
import Math.floor;

class time_type{
//! IS IT NECESSARY TO IMPLEMENT

// private:
    // const days_per_month: [0..11] int = {31,28,31,30,31,30,31,31,30,31,30,31};
    var seconds_per_day = 86400;
    var days_in_400_year_period = 146097;
    var ticks_per_second = 1;

    var seconds: int;
    var days: int;
    var ticks: int;

// public:
    proc init(){
        this.seconds = 0;
        this.days = 0;
        this.ticks = 0;
    }
    proc init(in seconds:int, in days:int = 0, in ticks:int = 0){
        this.seconds = seconds;
        this.days = days;
        this.ticks = ticks;
    }

    proc get_time(in ret_day:bool = true){
        // return a tuple (days,seconds,ticks)
        // TODO: fix this so we could enable it to not return days 
        var seconds,days,ticks: int;

        seconds = this.seconds;
        ticks = this.ticks;

        if(!ret_day) then seconds = seconds + this.days*seconds_per_day;

        if(ret_day) then return (this.days,seconds,ticks);
        else return (0,seconds,ticks);

        //TODO: PUT CHECKS HERE FOR OVERFLOW OF INT64
        
    }
    proc set_time(in seconds:int = 0,in days:int = 0,in ticks:int = 0){
        var seconds_new, days_new, ticks_new:int;
        // var new_time:time_type;

        seconds_new = seconds + (floor(ticks/ticks_per_second:real)):int;
        ticks_new = mod(ticks,ticks_per_second); 
        days_new = days + (floor(seconds_new/seconds_per_day:real)):int;
        seconds_new = mod(seconds_new,seconds_per_day);

        if(seconds_new < 0 || ticks_new < 0){
            writeln("[Function set_time] Bad Result");
            exit();
        }
        if(days_new < 0){
            writeln("[Function set_time] ?? Time is negative ??");
        }else{
            this.days = days_new;
            this.seconds = seconds_new;
            this.ticks = ticks_new;
        }
    }

    proc increment_time(in seconds:int,in days:int = 0,in ticks:int = 0,in allow_neg_inc_local:bool = false){
        if(!allow_neg_inc_local && (seconds < 0 || days < 0 || ticks < 0)){
            writeln("[Function increment_time] One or More time increments are negative");
        }

        // var t_new: time_type = new time_type();
        this.set_time(this.seconds + seconds, this.days + days, this.ticks + ticks);
        // return t_new;
    }

    proc decrement_time(in seconds:int,in days:int = 0,in ticks:int = 0){
        this.increment_time(-seconds, -days, -ticks, true);
    }

    proc +(in t1:time_type,in t2:time_type){
        return this.increment_time(t1,t2.seconds,t2.days,t2.ticks);
    }

    proc -(in t1:time_type,in t2:time_type){
        return this.increment_time(t1,-t2.seconds,-t2.day,-t2.ticks,true);
    }
}

// test 1: set_time and get_time w/o ticks
// test 2: set_time and get_time w/ ticks
// test 3: + operator check
// test 4: test of increment_time and decrement_time
// test 5: test of negative increments in increment_time and decrement_time
// test 6: test of negative seconds and/or ticks

proc test(){
    writeln("Testing Functionalities");
    writeln("Running Test - 1");

    var time: time_type = new time_type(seconds = 2, days = 1);
    time.ticks_per_second = 10;
    var res1 = time.get_time();
    writeln(res1);
    res1 = time.get_time(false);
    writeln(res1);
    writeln("-------------------------------------------------------");

    writeln("Running Test - 2");
    time.set_time(seconds = 10, days = 2, ticks = 10);
    res1 = time.get_time();
    writeln(res1);
    writeln("-------------------------------------------------------");

    // writeln("Running Test - 3");
    // time.set_time(seconds = 0, days = 2, ticks = 5);
    // var time2: time_type = new time_type();
    // time2.set_time(0,2,6);
    // writeln((time + time2).get_time());
    // writeln("-------------------------------------------------------");

    writeln("Running Test - 4");
    time.set_time(seconds=0, days=2);
    time.increment_time(seconds=0, days=1);
    writeln(time.get_time());
    time.decrement_time(seconds=0, days=2);
    writeln(time.get_time());
    time.set_time(seconds=0, days=2, ticks=5);
    time.increment_time(seconds=400, days=1, ticks=14);
    writeln(time.get_time());
    time.decrement_time(seconds=400, days=1, ticks=14);
    writeln(time.get_time());
    writeln("-------------------------------------------------------");

    writeln("Running Test - 5");
    time.set_time(seconds=0, days=2);
    time.increment_time(seconds=0, days=-1,allow_neg_inc_local = true);
    writeln(time.get_time());
    time.decrement_time(seconds=0, days=-1);
    writeln(time.get_time());
    time.set_time(seconds=0, days=2, ticks=5);
    time.increment_time(seconds=-400, days=-1, ticks=-14,true);
    writeln(time.get_time());
    time.decrement_time(seconds=-400, days=-1, ticks=-14);
    writeln(time.get_time());
    writeln("-------------------------------------------------------");

    writeln("Running Test - 6");
    time.set_time(seconds=-86399, days=2, ticks=-1);
    writeln(time.get_time());
    time.set_time(-86399,2,-9);
    writeln(time.get_time());
    time.set_time(86400,2,9);
    writeln(time.get_time());

}