module DataArray {
    pragma "no doc"
    enum DType {Int64, Real64, Bool, String, Undefined};

    pragma "no doc"
    proc toDType(type t) {
        select t {
            when int do
                return DType.Int64;
            when real do 
                return DType.Real64;
            when bool do 
                return DType.Bool;
            when string do 
                return DType.String;
            otherwise
                return DType.Undefined;
        }
    }

    class AbstractDataArray {
        pragma "no doc"
        var d_type: DType;

        pragma "no doc"
        const rank: int;

        pragma "no doc"
        proc init(type d_type, rank: int) {
            this.d_type = toDType(d_type);
            this.rank = rank;
        }

        pragma "no doc"
        proc _add(lhs): owned AbstractDataArray {
            halt("Pure virtual method");
        }

        proc add(rhs: borrowed AbstractDataArray): owned AbstractDataArray {
            halt("Pure virtual method");
        }

        pragma "no doc"
        proc _subtract(lhs): owned AbstractDataArray {
            halt("Pure virtual method");
        }

        proc subtract(rhs: borrowed AbstractDataArray): owned AbstractDataArray {
            halt("Pure virtual method");
        }

        pragma "no doc"
        proc _eq(lhs): bool {
            halt("Pure virtual method");
        }

        proc eq(rhs: borrowed AbstractDataArray): bool {
            halt("Pure virtual method");
        }

        pragma "no doc"
        proc convertTo(type to_convert): owned AbstractDataArray {
            halt("Pure virtual method");
        }
    }

    /*
    A DataArray is a container over an array with additional support for dimension labels for each axis of the array.
    */
    class DataArray: AbstractDataArray {
        /* The type of elements contained in the DataArray. */
        type eltType;

        /* A int paramter indicating the rank of the DataArray. */
        param rank: int;

        /* A bool parameter indicating whether any of the domain’s dimensions will be characterized by a strided range */
        param stridable: bool;

        /* The domain of the Chapel array owned by the DataArray.
           Refer to `the Chapel docs <https://chapel-lang.org/docs/primers/domains.html>`_ for more information on domains. 
        */
        var dom: domain(rank, stridable = stridable);
        var arr: [dom] eltType;

        /* An index space indicating the labels for each of the axis of the array owned by the DataArray. */
        var dimensions: domain(string);

        /* Initializes a DataArray with the given size and value as the default value of the given type. */
        proc init(type eltType, size: domain, dimensions: domain(string)) where isDefaultInitializable(eltType) {
            super.init(eltType, size.rank);
            this.eltType = eltType;
            this.rank = size.rank;
            this.stridable = size.stridable;
            this.dom = size;

            var arr: [size] eltType;
            this.arr = arr;

            this.dimensions = dimensions;
        }

        /* Initializes a DataArray with the given size and value as the given value. */
        proc init(size: domain, dimensions: domain(string), in default_value) where isDefaultInitializable(default_value.type) {
            super.init(default_value.type, size.rank);
            this.eltType = default_value.type;
            this.rank = size.rank;
            this.stridable = size.stridable;
            this.dom = size;
            
            var arr: [size] eltType = default_value;
            this.arr = arr;

            this.dimensions = dimensions;
        }

        /* Initializes a DataArray with the given array. */
        proc init(in arr, dimensions: domain(string)) {
            super.init(arr.eltType, arr.domain.rank);
            this.eltType = arr.eltType;
            this.rank = arr.domain.rank;
            this.stridable = arr.domain.stridable;
            this.dom = arr.domain;

            this.arr = arr;

            this.dimensions = dimensions;
        }

        pragma "no doc"
        override proc _add(lhs: borrowed DataArray): owned AbstractDataArray where this.rank == lhs.rank && isOperable(lhs, this) {
            var rhs: borrowed DataArray = this;
            var arr = lhs.arr + rhs.arr;
            return new owned DataArray(arr, lhs.dimensions);
        }

        /* Utility method to add two ``DataArray``. 
        
          **Note**: This can be chained with other operators. 
        */
        override proc add(rhs: borrowed AbstractDataArray): owned AbstractDataArray {
            return rhs._add(this);
        }

        pragma "no doc"
        override proc _subtract(lhs: borrowed DataArray): owned AbstractDataArray where this.rank == lhs.rank && isOperable(lhs, this) {
            var rhs: borrowed DataArray = this;
            var arr = lhs.arr - rhs.arr;
            return new owned DataArray(arr, lhs.dimensions);
        }

        /* Utility method to subtract two ``DataArray``. 
        
          **Note**: This can be chained with other operators. 
        */
        override proc subtract(rhs: borrowed AbstractDataArray): owned AbstractDataArray {
            return rhs._subtract(this);
        }

        pragma "no doc"
        override proc _eq(lhs: borrowed DataArray): bool where isOperable(lhs, this) {
            var rhs: borrowed DataArray = this;

            if lhs.arr._value == rhs.arr._value then
                return true;

            if lhs.rank != rhs.rank then
                return false;

            if lhs.dom.size: uint != rhs.dom.size: uint then
                return false;

            if isRectangularDom(lhs.dom) && isRectangularDom(rhs.dom) {
                for d in 0..#rhs.rank do
                    if rhs.dom.dim(d).size: uint != lhs.dom.dim(d).size: uint then
                        return false;
            }

            var ret = true;
            forall (l, r) in zip(lhs.arr, rhs.arr) with (&& reduce ret) do
                ret &&= (l == r);
            return ret;
        }

        /* Utility method to check if two ``DataArray`` are equal. */
        override proc eq(rhs: borrowed AbstractDataArray): bool {
            return rhs._eq(this);
        }
    }

    pragma "no doc"
    proc isOperable(lhs: borrowed DataArray, rhs: borrowed DataArray) param {
        if !isPrimitive(lhs.eltType) && !isPrimitive(rhs.eltType) {
            return lhs.arr[lhs.dom.alignedLow].dims(rhs.arr[rhs.dom.alignedLow]);
        }
        return isCoercible(lhs.eltType, rhs.eltType);
    }

    operator +(lhs: borrowed AbstractDataArray, rhs: borrowed AbstractDataArray): owned AbstractDataArray {
        return lhs.add(rhs);
    }

    operator -(lhs: borrowed AbstractDataArray, rhs: borrowed AbstractDataArray): owned AbstractDataArray {
        return lhs.subtract(rhs);
    }

    operator ==(lhs: borrowed AbstractDataArray, rhs: borrowed AbstractDataArray): bool {
        return lhs.eq(rhs);
    } 

    operator !=(lhs: borrowed AbstractDataArray, rhs: borrowed AbstractDataArray): bool {
        return !(lhs == rhs);
    }    
}