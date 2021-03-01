# Static_Pricing_problem
Codes related to the paper: A. Marandi, V Lurkin (2020), Static Pricing Problems under Mixed Multinomial Logit Demand

In case of using these codes, you are obliged to properly cite the paper.  

## Table of contents
* [Technologies](#technologies)
* [Data](#Data)
* [Functions and their features](#Functions-and-their-features)

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

## Data
This package contains many data related to pricing problems related to parking choice. In the Data.jl file, we store all the used data in the paper and in this section we discuss how to use them.
All the data are store in a function format. So, in order to generate the parameters, you may need to call the function with suitable arguments. 
* ``` Logit_10() ```: 


