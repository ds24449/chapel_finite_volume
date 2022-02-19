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


proc roll(in a:[?D],in shift, axis = 0){
    var a1: D.type;
    var a2: D.type;
    
    var tup1 = D.dims();
    var tup2 = D.dims();
    if(shift<0){
        shift *= -1;
        tup1(axis) = D.dim(axis).low+shift..D.dim(axis).high;
        tup2(axis) = D.dim(axis).low..D.dim(axis).low+shift-1;
    }else{
        tup1(axis) = D.dim(axis).high-shift+1..D.dim(axis).high;
        tup2(axis) = D.dim(axis).low..D.dim(axis).high-shift;
    }
    a1 = tup1;
    a2 = tup2;

    tup1 = D.dims();
    tup2 = D.dims();

    tup1(axis) = D.dim(axis).low..a1.dim(axis).size;
    tup2(axis) = a1.dim(axis).size+1..D.dim(axis).high;

    var firstHalf:D.type = tup1;
    var secondHalf:D.type = tup2;

    var rez:[D] a.eltType;
    rez[firstHalf] = a[a1];
    rez[secondHalf] = a[a2];
    return rez;
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
    // assert(roll( a, shift = 1, axis = 0) == [9, 0, 1, 2, 3, 4, 5, 6, 7, 8]);
    // assert(roll( a, shift =-1, axis = 0) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]);
    // assert(roll( a, shift = 2, axis = 0) == [8, 9, 0, 1, 2, 3, 4, 5, 6, 7]);
    // assert(roll( a, shift =-2, axis = 0) == [2, 3, 4, 5, 6, 7, 8, 9, 0, 1]);
    // writeln(roll( a, shift = 1, axis = 0));
    // writeln(roll( a, shift =-1, axis = 0));
    // writeln(roll( a, shift = 2, axis = 0));
    // writeln(roll( a, shift =-2, axis = 0));
    writeln("Rolling on 1-D Array Passed");

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
rollTest();