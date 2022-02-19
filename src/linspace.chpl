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
    var a1: D.type;
    var a2: D.type;
    
    if(a.rank == 1){
        
        if(shift<0){
            shift *= -1;
            a1 = {D.low+shift..D.high};
            a2 = {D.low..D.low+shift-1};
        }else{
            a1 = {D.high-shift+1..D.high};
            a2 = {D.low..D.high-shift};
        }

        var rez:[D] a.eltType;
        rez[D.low..a1.size] = a[a1];
        rez[a1.size+1..D.high] = a[a2];
        return rez;
    }
    else if(a.rank == 2){
        // TODO: Code for 2D arrays, if possible ND arrays
        var first_half:domain(2);
        var second_half:domain(2);
        if(axis == 0){
            if(shift<0){
                shift *= -1;
                a1 = {a.dim(0).low+shift..a.dim(0).high, a.dim(1)};
                a2 = {a.dim(0).low..a.dim(0).low+shift-1, a.dim(1)};
            }else{
                a1 = {a.dim(0).high-shift+1..a.dim(0).high, a.dim(1)};
                a2 = {a.dim(0).low..a.dim(0).high-shift, a.dim(1)};
            }
            second_half = {a1.dim(0).size..a.dim(0).high,a.dim(1)};
            first_half = {a.dim(0).low..a1.dim(0).size-1,a.dim(1)};
        }
        else{
            if(shift<0){
                shift *= -1;
                a1 = {a.dim(0), a.dim(1).low+shift..a.dim(1).high};
                a2 = {a.dim(0), a.dim(1).low..a.dim(1).low+shift-1};
            }else{
                a1 = {a.dim(0), a.dim(1).high-shift+1..a.dim(1).high};
                a2 = {a.dim(0), a.dim(1).low..a.dim(1).high-shift};
            }
            second_half = {a.dim(0),a1.dim(1).size..a.dim(1).high};
            first_half = {a.dim(0),a.dim(1).low..a1.dim(1).size-1};
        }
        // writeln("[DEBUG - I]" + first_half:string + " " +a1:string + second_half:string + " " +a2:string);
        var rez:[D] a.eltType;
        rez[first_half] = a[a1];
        rez[second_half] = a[a2];
        return rez;
    }

}

proc like_ones(in D,type eltType = real){
    var res:[D] eltType = 1;
    return res;
}

proc min_of_arr(A:[?D]) {
/*
	Calculate the smallest element in the array
    	Parameters
    	----------
    	A : array_like
        	Input array.

*/
	var (minVal, minLoc) = minloc reduce zip(A, A.domain);
	return minVal;
}

proc max_of_arr(A:[?D]) {
/*
	Calculate the largest element in the array
    	Parameters
    	----------
    	A : array_like
        	Input array.
*/
	var (maxVal, maxLoc) = maxloc reduce zip(A, A.domain);
	return maxVal;
}

proc np_min(A:[?d1], B:[?d2]) {
/*
	Calculates the minimum element from the 2 arrays
    	Parameters
    	----------
    	A : array_like
        	Input array.
        B : array_like
            Input array-2
*/
	var min_a = min_of_arr(A);
	var min_b = min_of_arr(B);
	var res = if min_a < min_b then min_a else min_b;
	return res;
}

proc np_max(A:[?d1], B:[?d2]) {
/*
	Calculates the maximum element from the 2 arrays
    	Parameters
    	----------
    	A : array_like
        	Input array-1
        B : array_like
            Input array-2
*/
	var max_a = max_of_arr(A);
	var max_b = max_of_arr(B);
	var res = if max_a > max_b then max_a else max_b;
	return res;
}

// ------------------------------- * TESTS * -----------------------------------

proc rollTest(){
    var a:[1..10] int = [0,1,2,3,4,5,6,7,8,9];
    assert(roll( a, shift = 1, axis = 0) == [9, 0, 1, 2, 3, 4, 5, 6, 7, 8]);
    assert(roll( a, shift =-1, axis = 0) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]);
    assert(roll( a, shift = 2, axis = 0) == [8, 9, 0, 1, 2, 3, 4, 5, 6, 7]);
    assert(roll( a, shift =-2, axis = 0) == [2, 3, 4, 5, 6, 7, 8, 9, 0, 1]);
    writeln(roll( a, shift = 1, axis = 0));
    writeln(roll( a, shift =-1, axis = 0));
    writeln(roll( a, shift = 2, axis = 0));
    writeln(roll( a, shift =-2, axis = 0));
    // writeln("Rolling on 1-D Array Passed");

    var a2 = reshape(a,{1..2,1..5});
    // assert(roll( a2, shift = 1, axis = 0) == [[5, 6, 7, 8, 9], [0, 1, 2, 3, 4]]);
    // assert(roll( a2, shift =-1, axis = 0) == [[5, 6, 7, 8, 9], [0, 1, 2, 3, 4]]);
    // assert(roll( a2, shift = 1, axis = 1) == [[4, 0, 1, 2, 3], [9, 5, 6, 7, 8]]);
    // assert(roll( a2, shift =-1, axis = 1) == [[1, 2, 3, 4, 0], [6, 7, 8, 9, 5]]);
    writeln(roll( a2, shift = 1, axis = 0));
    writeln(roll( a2, shift =-1, axis = 0));
    writeln(roll( a2, shift = 1, axis = 1));
    writeln(roll( a2, shift =-1, axis = 1));
}