iter linspace(type dtype, start, stop, num, in endpoint:bool=true) {
    assert(num > 0, "number of points must be > 0");
    if num==1 then endpoint = false;
    const ninterval = (if endpoint then num-1 else num):real;
    const dt = (stop-start)/ninterval;

    for i in 0..#num do yield (start+i*dt):dtype;
}
iter linspace(start, stop, num, endpoint=true) {
    for x in linspace(real(64), start, stop, num, endpoint) do yield x;
}


proc roll(a:[?D],in shift, axis = 0){
    /*
    Roll array elements along a given axis.
    Elements that roll beyond the last position are re-introduced at
    the first.
    Parameters
    ----------
    a : array_like
        Input array.
    shift : int or tuple of ints
        The number of places by which elements are shifted.  If a tuple,
        then `axis` must be a tuple of the same size, and each of the
        given axes is shifted by the corresponding number.  If an int
        while `axis` is a tuple of ints, then the same value is used for
        all given axes.
    axis : int or tuple of ints, optional
        Axis or axes along which elements are shifted.  By default, the
        array is flattened before shifting, after which the original
        shape is restored.
    Returns
    -------
    res : ndarray
        Output array, with the same shape as `a`.
    */
    var a1:domain(a.rank);
    var a2:domain(a.rank);
    
    if(a.rank == 1){
        
        if(shift<0){
            shift *= -1;
            a1 = {D.low+shift+1..D.high};
            a2 = {D.low..D.low+shift};
        }else{
            a1 = {D.high-shift..D.high};
            a2 = {D.low..D.high-shift-1};
        }

        var rez:[D] a.eltType;
        rez[D.low..a1.size] = a[a1];
        rez[a1.size+1..D.high] = a[a2];
        return rez;
    }
    else if(a.rank == 2){
        // TODO: Code for 2D arrays, if possible ND arrays
        if(axis == 0){
            if(shift<0){
                shift *= -1;
                a1 = {a.dim(0).low+shift+1..a.dim(0).high, a.dim(1)};
                a2 = {a.dim(0).low..a.dim(0).low+shift, a.dim(1)};
            }else{
                a1 = {a.dim(0).high-shift..a.dim(0).high, a.dim(1)};
                a2 = {a.dim(0).low..a.dim(0).high-shift-1, a.dim(1)};
            }
        }
        else{
            if(shift<0){
                shift *= -1;
                a1 = {a.dim(0), a.dim(1).low+shift+1..a.dim(1).high};
                a2 = {a.dim(0), a.dim(1).low..a.dim(1).low+shift};
            }else{
                a1 = {a.dim(0), a.dim(1).high-shift..a.dim(1).high};
                a2 = {a.dim(0), a.dim(1).low..a.dim(1).high-shift-1};
            }
        }
        var rez:[D] a.eltType;
        rez[a.dim(0).low..a1.dim(0).size,a.dim(1).low..a1.dim(1).size] = a[a1];
        rez[a1.dim(0).size+1..a.dim(0).high,a1.dim(1).size+1..a.dim(1).high] = a[a2];
        return rez;
    }

}


proc rollTest(){
    var a:[1..10] int = [0,1,2,3,4,5,6,7,8,9];
    assert(roll( a, shift = 1, axis = 0) == [8, 9, 0, 1, 2, 3, 4, 5, 6, 7]);
    assert(roll( a, shift =-1, axis = 0) == [2, 3, 4, 5, 6, 7, 8, 9, 0, 1]);
    writeln("Rolling on 1-D Array Passed");

    var a2 = reshape(a,{1..2,1..5});
    // writeln(a2.dim(1).high);
    writeln(roll( a2, shift = 1, axis = 0));
    // roll( a2, shift =-1, axis = 1);
}