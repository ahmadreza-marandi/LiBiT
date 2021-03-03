# Static_Pricing_problem
Codes related to the paper: A. Marandi, V Lurkin (2020), Static Pricing Problems under Mixed Multinomial Logit Demand (https://arxiv.org/abs/2005.07482).

In case of using these codes, you are obliged to properly cite the paper.  

## Table of contents
* [Technologies](#technologies)
* [Data](#Data)
* [Functions and their features](#Functions-and-their-features)
* [Examples](#Examples)

## Technologies
For this project, we use multiple Julia packages. The code has been updated and tested on Julia 1.5.3. The packages and the tested versions are listed below:
* ForwardDiff v0.10.15 (https://github.com/JuliaDiff/ForwardDiff.jl)
* LinearAlgebra (https://github.com/JuliaLang/julia/tree/master/stdlib/LinearAlgebra)
* JuMP v0.21.5 (https://jump.dev/)
* CPLEX v0.6.6 (https://github.com/jump-dev/CPLEX.jl): make sure to properly install IBM ILOG CPLEX Optimization Studio beforehand
* MAT v0.8.1 (https://github.com/JuliaIO/MAT.jl): used to import and export some parameters
* CPUTime v1.0.0 (https://github.com/schmrlng/CPUTime.jl)
* JLD v0.10.0 (https://github.com/JuliaIO/JLD.jl): used for some importing and exporting
* Distributions v0.23.12 (https://github.com/JuliaStats/Distributions.jl)
* Cuba v2.2.0 (https://github.com/giordano/Cuba.jl)
* FiniteDiff v2.8.0 (https://github.com/JuliaDiff/FiniteDiff.jl)
* NLopt v0.6.2 (https://github.com/JuliaOpt/NLopt.jl)
* SCIP v0.9.6 (https://github.com/scipopt/SCIP.jl): make sure to instal SCIP beforehand.

## Data
This package contains many data related to pricing problems related to parking choice. In the Data.jl file, we store all the used data in the paper and in this section we discuss how to use them. All the codes return the following outputs:
* Beta_parameter: a matrix whose rows are related to the parking choices (FSP, PSP, and PUP) and columns are related to customer classes. This matrix represents $`\beta^p_{in}`$.
* q_parameter: a matrix whose rows are related to the parking choices (FSP, PSP, and PUP) and columns are related to customer classes. This matrix represents $`q_{in}`$.
* NUM_POINTS: This is the diemension of the problem. 
* N: This is the number of customer classes.
* UB_p: This is the vector containing the upper bounds on the prices
* LB_p: This is the vector containing the lower bounds on the prices.

All the data are store in a function format. So, in order to generate the parameters, you may need to call the function with suitable arguments. 
* ``` Logit_10() ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for MNL model with N=10;
* ``` Logit_50() ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for MNL model with N=50;
* ``` Mixed_Logit_10(β) ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for Mixed Logit model with N=10, given the value β.
* ``` Mixed_Logit_50(β) ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for Mixed Logit model with N=50, given the value β.
* ``` Mixed_Logit_n10_random(R_AT) ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p) for discrete Mixed Logit model with R_AT number of points with N=10 customer classes. For this function R_AT can be 10, 100, 1000, 10000, 100000, 1000000. 


## Functions and their features
This package contains many functions. In this section we discuss each functions, what they do, and how they can be used:
* ``` discrete_Mixed_logit_function(p) ``` returns the objective value of the Discrete Mixed Logit model at point p. 
* ``` Distribution_function(p) ``` returns the objective value of the MNL model at point p.
* ``` homogenious_MNL(time_limit,data_function,R_AT) ``` returns a 2-tuple. Input: time_limit is the limit on the time to run the code, data_function is one of the functions in Data.jl about discrete choice models (for instance Mixed_Logit_n10_random), and R_AT is the number of point in the discrete choice model.  Output: The first element is the objective value and the second element is the obtained solution. This function is the implementation of Algorithm 2 in [Li et al. (2019)](https://pubsonline.informs.org/doi/abs/10.1287/msom.2017.0675). 
* ``` SCIP_degenerate(function_data, time_limit) ``` returns a 4-tuple. Input: function_data is the data related to a MNL model (for instance Logit_10 in Data.jl), and time_limit is the time limit for SCIP to solve the pricing problem. The first element of the output is the time taken by SCIP, the second element is the termination status of SCIP, the third element is the optimality gap, and the last element is the obtained solution by SCIP.
* ``` Mixed_logit_function(p::Vector) ``` returns the objective value of the pricing problem with a continuous logit model. To use this function, the following considerations are needed:
  * this function uses the function ``` Mixed_Logit_distribution(p,β,i) ```. The data of the logit model should be put here. In the first line of this function, the data is acquared (for instance Mixed_Logit_50(β)). Then, the information about the distribution of β should be given. Currently we use the normal distrbution Normal(μ, Σ) using the code ``` pdf(MvNormal(μ, Σ), β) ```.
  * to take the derivative, we need a bound. These are the vectors ``` a ``` (lower bounds on β) and ``` b ``` (upper bound on β) inside the function ``` Mixed_logit_function_i(p, i) ```.
* ``` Mixed_logit_function_discreteDis(p::Vector) ``` returns the value of the discretize relaxation of the continuous mixed logit at point p. For this fucntion the following considerations are needed:
   * all the information about the data is given in the function ``` Mixed_logit_function_i_discreteDis(p::Vector, i) ```. In this function, we currently use the data for ``` Mixed_logit_50 ``` function where the box determined with  ``` a ``` (lower bounds on β) and ``` b ``` (upper bound on β) are discritized by $`R^2`$ points. You need to specify ``` R ``` in this function (default value is 10);
* ```  cuts(Scenario)  ``` retunrs a 2-tuple. This functions generate the cut to be used in Voronoi diagram. Input: set of scenario in the form of ``` Array{Array{Float64,2},1} ``` (each component is a scenario); Output: first output is a matrix A and the second output is the vector b together generate cuts in the form Aξ⫺b.
* ``` solve_local_opt(type_of_problem,time_limit) ``` returns a 4-tuple. This fucntion is the implimentation of LiBiT in our paper. Input: type_of_problem is string with values either MNL, Mixed, or Discrete, time_limit is the time limit in seconds for LiBiT. Output: element 1 is the lower bound on the optimal value, element 2 is the best obtain solution, element 3 is the upper bound on the optimal value, element 4 is the list containing the information extracted from LiBiT after each iteration.
   * MNL: uses the fucntions ``` Distribution_function ``` and ``` Distribution_function_i ```
   * Mixed: uses the functions ``` Mixed_logit_function ``` and ``` Mixed_logit_function_i ```
   * Discrete: uses the fucntions ``` Mixed_logit_function_discreteDis ``` and ``` Mixed_logit_function_i_discreteDis ```.
 
In LiBit, we use different methods. The rest of functions are the ones used each iteration of LiBiT:
* ``` trust_region_iteration(A,b,Tolerr,obj_fun::Function ,initial_point,time_limit) ```: In LiBiT we use the trust region approach. This function performs one iteration of the trust region approach. In other words, it gets matrix A and vector b (to define the feasible region Ap⩾b), and the objective function ``` obj_fun ```. Given the initial solution ``` initial_point ``` and within time limt time_limit, it finds the best solution using trust region method with tolorence error of Tolerr. The function returns the solution obtain from one iteration of trust region method (2nd element) and the corresponding objective value (1st element). 
* ``` randomPoint(A,b,p_0,LB_p,UB_p) ```: in LiBiT, we may use random feasible point selection. So, this function returns a (random) solution inside feasible region Ap⩾b, and given vectors of lower bound LB_p and upper bound UB_p (for p). The fucntion trys to find either a point close to the p_0 on the other half of the feasible region or a random point. 
* ``` overestimator_general(A,b,P_0::Array{Array{Float64,2},1},obj_fun::Function,LB_τ,UB_τ,LB_p,UB_p,N,time_limit) ``` is the function related to problem (12). 
   * Ax⫺b, p⩾0 is the feasible region
   * P_0 is the ser of points 
   * obj_fun is the objective function
   * LB_τ,UB_τ,LB_p,UB_p are the corresponding lower and upper bounds of p and τ
   * N is the number of customer classes
   * time_limit is the limit on time in seconds

   The output will be the upper bound on the optimal value. 
   
