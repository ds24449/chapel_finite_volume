module FDSolver{
    use DataArray;
    use List;

    class FDSolver{
        const weights = [
            [(0.0,0.0,0.0,-1.0/2,0.0,1.0/2,0.0,0.0,0.0),(0.0,0.0,1.0/12,-2.0/3,0.0,2.0/3,-1.0/12,0.0,0.0),(0.0,-1.0/60,3.0/20,-3.0/4,0.0,3.0/4,-3.0/20,1.0/60,0.0),(1.0/280,-4.0/105,1.0/5,-4.0/5,0.0,4.0/5,-1.0/5,4.0/105,-1.0/280)],
            [(0.0,0.0,0.0,1.0,-2.0,1.0,0.0,0.0,0.0),(0.0,0.0,-1.0/12,4.0/3,-5.0/2,4.0/3,-1.0/12,0.0,0.0),(0.0,1.0/90,-3.0/20,3.0/2,-49.0/18,3.0/2,-3.0/20,1.0/90,0.0),(-1.0/560,8.0/315,-1.0/5,8.0/5,-205.0/72,8.0/5,-1.0/5,8.0/315,-1.0/560)],
            [(0.0,0.0,-1.0/2,1.0,0.0,-1.0,1.0/2,0.0,0.0),(0.0,1.0/8,-1.0,13.0/8,0.0,-13.0/8,1.0,-1.0/8,0.0),(-7.0/240,3.0/10,-169.0/120,61.0/30,0.0,-61.0/30,169.0/120,-3.0/10,7.0/240),(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)]
        ];
        const forward_wts = [
            [(-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),(-3.0/2,2.0,-1.0/2,0.0,0.0,0.0,0.0,0.0,0.0),(-11.0/6,3.0,-3.0/2,1.0/3,0.0,0.0,0.0,0.0,0.0),(-25.0/12,4.0,-3.0,4.0/3,-1/4.0,0.0,0.0,0.0,0.0),(-137.0/60,5.0,-5.0,10.0/3,-5.0/4,1.0/5,0.0,0.0,0.0),(-49.0/20,6.0,-15.0/2,20.0/3,-15.0/4,6.0/5,-1.0/6,0.0,0.0)],
            [(1.0,-2.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0),(2.0,-5.0,4.0,-1.0,0.0,0.0,0.0,0.0,0.0),(35.0/12,-26.0/3,19.0/2,-14.0/3,11.0/12,0.0,0.0,0.0,0.0),(15.0/4,-77.0/6,107.0/6,-13.0,61.0/12,-5.0/6,0.0,0.0,0.0),(203.0/45,-87.0/5,117.0/4,-254.0/9,33.0/2,-27.0/5,137.0/180,0.0,0.0),(469.0/90,-223.0/10,879.0/20,-949.0/18,41.0,-201.0/10,1019.0/180,-7.0/10,0.0)],
            [(-1.0,3.0,-3.0,1.0,0.0,0.0,0.0,0.0,0.0),(-5.0/2,9.0,-12.0,7.0,-3.0/2,0.0,0.0,0.0,0.0),(-17.0/4,71.0/4,-59.0/2,49.0/2,-41.0/4,7.0/4,0.0,0.0,0.0),(-49.0/8,29.0,-461.0/8,62.0,-307.0/8,13.0,-15.0/8,0.0,0.0),(-967.0/120,638.0/15,-3929.0/40,389.0/3,-2545.0/24,268.0/5,-1849.0/120,29.0/15,0.0),(-801.0/80,349.0/6,-18353.0/120,2391.0/10,-1457.0/6,4891.0/30,-561.0/8,527.0/30,-469.0/240)]
        ];
        const backward_wts = [
            [(0.0,0.0,0.0,-1.0,1.0),(0.0,0.0,1.0/2,-2.0,3.0/2)],
            [(0.0,0.0,1.0,-2.0,1.0),(0.0,-1.0,4.0,-5.0,2.0)],
            [(0.0,-1.0,3.0,-3.0,1.0),(3.0/2,-7.0,12.0,-9.0,5.0/2)]
        ];

        var lcl_cpy: DataArray;
        var lcl_dom: domain;

        proc init(const original: DataArray){
            this.lcl_cpy = new DataArray(original.arr, original.dimensions);
            this.lcl_dom = original.dom;
        }

        proc apply_bc(const ax:int, lt_bc:string, rt_bc:string, accuracy = 1){
            var old_dom = this.lcl_dom;

            if(lt_bc == "periodic"){
                for i in 1..accuracy{
                    if(this.lcl_cpy.rank == 1){
                        this.lcl_cpy.arr[this.lcl_dom.low - i] = this.lcl_cpy.arr[this.lcl_dom.high + 1 - i];
                        this.lcl_cpy.arr[this.lcl_dom.high + i] = this.lcl_cpy.arr[i];
                    }
                    else if(this.lcl_cpy.rank == 2){
                        var n = old_dom.high[ax];
                        var left_b:domain(2);
                        var left_orig:domain(2);

                        var right_b:domain(2);
                        var right_orig:domain(2);

                        if(ax == 0){
                            left_b = {1-i..1-i,1..n};
                            left_orig = {i..i,1..n};

                            right_b = {n+i..n+i,1..n};
                            right_orig = {n-i+1..n-i+1,1..n};
                        }else{
                            left_b = {1..n,1-i..1-i}; // These values are to be replaced by right_orig
                            left_orig = {1..n,i..i}; // These values will replace right_b

                            right_b = {1..n,n+i..n+i};
                            right_orig = {1..n,n-i+1..n-i+1};
                        }
                        this.lcl_cpy.arr[left_b] = this.lcl_cpy.arr[right_orig];
                        this.lcl_cpy.arr[right_b] = this.lcl_cpy.arr[left_orig];
                    }
                }
            }
            
        }

        proc apply_bc(const dict,accuracy = 2){ //TODO: Take an Associative Array {co-ordinate: bounds type} 
            // X co-ordinate {'x' => "periodic"};
            this.lcl_cpy.dom = this.lcl_cpy.dom.expand(accuracy);
            for d in dict.domain{
                var axis:int;
                if(dict.domain.idxType == string){
                    for (i,j) in zip(0..this.lcl_cpy.dimensions.rank,this.lcl_cpy.dimensions){
                        if(d == j){
                            axis = i;
                            break;
                        }
                    }
                }else{
                    axis = d;
                }

                if(dict[d].type == 2*string){
                    var left_bc = dict[d][0];
                    var right_bc = dict[d][1];
                    apply_bc(axis,left_bc,right_bc,accuracy);
                }
                else{
                    var bc = dict[d];
                    apply_bc(axis,bc,bc,accuracy);
                }
            }   
        }

        proc derivative(const weight,const extent,const axis:int = 0){ 
            // if(weight.size != extent.size) then writeln("Weight and Extent length Mis-Match in derivative function");

            var data:DataArray = new owned DataArray(this.lcl_cpy.arr[this.lcl_dom],this.lcl_cpy.dimensions); // change the domain here

            if(this.lcl_dom.rank == 1){
                forall i in this.lcl_dom{
                    data.arr[i] = 0;
                    for (k,j) in zip(weight,extent){
                        data.arr[i] += k*this.lcl_cpy.arr[i+j];
                    }
                }
            }
            else{
                forall i in this.lcl_dom{
                    var sum:real = 0.0;

                    for (k,j) in zip(weight,extent){
                        var temp = i;
                        temp[axis] += j;
                        sum += (k*this.lcl_cpy.arr[temp]);
                    }
                    data.arr[i] = sum;
                }
            }
            //data.arr.updateFluff();
            return data;
        } 

        proc Finite_Difference(scheme = "central",order:int(32) = 1,accuracy:int(32) = 2,step:real(64),axis:int = 0,debugFlag:bool=false){
            var wts:list(real(64));
            select scheme{
                when "forward" do{
                    var extnt_temp = 0..(accuracy+order-1);
                    for j in extnt_temp{
                        wts.append(forward_wts[order-1][accuracy-1][j]);
                    }
                    var temp = this.derivative(wts,extnt_temp,axis=axis);
                    temp.arr /= (step**(order));
                    return temp;
                }
                when "backward" do{
                    var extnt_temp = -(accuracy+order-1)..0;
                    for j in extnt_temp{
                        wts.append(backward_wts[order-1][accuracy-1][4+j]);
                    }
                    var temp = this.derivative(wts,extnt_temp,axis=axis);
                    temp.arr /= (step**(order));
                    return temp;
                }
                when "central" do{
                    var extnt_temp = -accuracy/2..accuracy/2;
                    for j in extnt_temp{
                        wts.append(weights[order-1][(accuracy-1)/2][4+j]);
                    }
                    var temp = this.derivative(wts,extnt_temp,axis=axis);
                    temp.arr /= (step**(order));
                    return temp;
                }
                otherwise{
                    writeln("Error: Wrong Scheme Name");
                    return new DataArray(this.lcl_cpy.arr,this.lcl_cpy.dimensions);
                }
            }
        }


    }
}