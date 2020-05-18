# Julia 1.0.5
#This document contains codes written for the paper
# A. Marandi, V Lurkin (2020), Static Pricing Problems under Mixed Multinomial Logit Demand
#written by Ahmadreza Marandi
# All rights are reserved.
#
#
using ForwardDiff
using LinearAlgebra
using JuMP
using CPLEX
using MAT
using CPUTime
using Cubature
using JLD
using Distributions
using Cuba
using FiniteDiff
using NLopt

struct MyProblem
	model
	p
	y
end

function eye(n)
	Matrix{Float64}(I,n,n);
end
function Distribution_function(p::Vector)
	#objective function related to MNL
	f=sum( p[i]*( exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/(sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS)))  for n=1:N, i=2:NUM_POINTS)		
	return f;
end
function homogenious_MNL(time_limit)
	# Algorithm 2 of the paper Li et al. (2019)
	Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p=Mixed_Logit_n10_r100();
	#considering one customer behaviour
	b=Array{Float64,2}(-Beta_parameter[:,1,:]);
	a=Array{Float64,2}(q_parameter[:,1,:]); 
	A=exp.(a);
	m=R;
	n=NUM_POINTS;
	w=1/R*ones(R,1);
	hatp_1=1 ./ b[:,1];
	hatp_1[1]=0;
	diff=Inf;
	Tolerr=1e-4;
	time_elapsed=0;
	maxf=Inf;
	while diff>Tolerr
		display(" ===$maxf===$hatp_1 =====time=$time_elapsed")
		CPUtic()
		d=zeros(NUM_POINTS,1);
		q_0=zeros(m,1);
		q=zeros(n,m);
		for k=1:m
			q_0[k]=1/(1 + sum( A[j,k] * exp( -b[j,k] * hatp_1[j] ) for j=2:n));
			for i=2:n
				q[i,k] = q_0[k] * A[i,k] * exp( -b[i,k] * hatp_1[i] );
			end
		end
		q_cum= sum(w[k] * q[:,k] for k=1:m);

		for i=2:n
			d[i]=(1 / ( sum( b[i,k] * w[k] * q[i,k] / q_cum[i] for k=1:m ) ) )+
					sum( (w[k] * b[i,k] * q[i,k] / (sum( w[ℓ] * b[i,ℓ] * q[i,ℓ] for ℓ=1:m ))) * (sum( hatp_1[j] * q[j,k] for j=2:n)) for k=1:m ) - hatp_1[i];
		end
		function π_obj(α)
			p=hatp_1+α*d;
			q_0=zeros(m,1);
			q=zeros(n,m);
			for k=1:m
				q_0[k]=1/(1 + sum( A[j,k] * exp( -b[j,k] * p[j] ) for j=2:n));
				for i=2:n
					q[i,k] = q_0[k] * A[i,k] * exp( -b[i,k] * p[i] );
				end
			end
			return sum( w[k] * sum(p[i] * q[i,k] for i=1:n) for  k=1:m);
		end
		opt = Opt(:GN_DIRECT_L, 1);
		opt.lower_bounds = [0];
		opt.upper_bounds = [1];
		# opt.xtol_rel = 1e-;
		opt.maxtime= time_limit-time_elapsed;
		opt.max_objective = π_obj;
		(maxf,α_p,ret) = optimize(opt, [1])
		maxf=π_obj(α_p[1]);
		hatp_1=hatp_1+α_p[1]*d;
		diff=norm(α_p[1]*d);
		time_elapsed = time_elapsed+ CPUtoc();
	end
end
function Mixed_Logit_distribution(p::Vector,β::Vector,i)
	#component of the continous mixed logit inside the integral
	#p is the cost Vector
	#β is the random variable used in the continous mixed logit
	#i is the index of the alternative service

	Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p=Mixed_Logit_10(β);# given β, reading the data.
	NUM_POINTS=size(p,1);
	μ=[-0.788;-32.3];# mean of the random varialbe
	Σ=[1.06 12/8;12/8 14.12];# covariance matrix of the random varialbe
	f=zeros(N,1);
	# for each customer, calculate the component function
	for n=1:N
		f[n]=(exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/ sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS));
	end
	
	return sum(f[n] for n=1:N)*pdf(MvNormal(μ, Σ), β)
end
function Mixed_logit_function_i(p::Vector, i)
	#component of the continous mixed logit related to the integral
	#p is the cost Vector
	#i is the index of the alternative service
	a=[-3.6;-68.52];# lowerbounds on the integral
	b=[1.94;3.92];#upperbounds on the integral

	f=cuhre((x,val)->val[1]=Mixed_Logit_distribution(p,a+x.*(b-a) ,i),2,1)# calculating the integral
	return 1/(prod(b-a)*sum(f[1]));#sum is to make the value scalar
end
function Mixed_logit_function(p::Vector)
	#calculating the objective function with respect to the continuous mixed logit
	# p is the vector of prices
	return sum(p[i]/Mixed_logit_function_i(p, i) for  i=1:NUM_POINTS)
end
function Mixed_logit_function_nlopt(p::Vector, grad::Vector)
	# the objective function used in NLopt package for continuous mixed logit
	if length(grad) > 0# gradient 
		g = x -> FiniteDiff.finite_difference_gradient(Mixed_logit_function, x);
		grad=g(p);
	end
	val=sum(p[i]/Mixed_logit_function_i(p, i) for i=1:NUM_POINTS);
	return val
end
function cuts(Scenario)
	#Scenario is the active scenarios; each point is in the row format without ","
	# returns A,b as the cut Aξ⫺b corresponding to the first scenario

	m=size(Scenario,1);

	# constructing the cuts
	A=(Scenario[2]-Scenario[1])';
	bhelp=(Scenario[1]+Scenario[2]);
	b=0.5.*(Scenario[2]-Scenario[1])'*bhelp;
	for i=3:m
		A=[A;(Scenario[i]-Scenario[1])'];
		bhelp=(Scenario[1]+Scenario[i]);
		b=[b;0.5.*(Scenario[i]-Scenario[1])'*bhelp];
	end
	return -A,-b;
end
function Distribution_function_NLopt(p::Vector, grad::Vector)
	# the objective function used in NLopt package for discrete mixed logit
	if length(grad) > 0 #gradient calculations
		g = x -> FiniteDiff.finite_difference_gradient(Distribution_function, x);
		grad=g(p);
	end
	return sum( p[i]*( exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/(sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS)))  for n=1:N, i=2:NUM_POINTS)
end

function Distribution_function_i(p::Vector,i)
	# objective fucntion of the discrete mixed logit
	f= 1/sum(1/sum(exp((Beta_parameter[j,n]*p[j]-Beta_parameter[i,n]*p[i])+q_parameter[j,n]-q_parameter[i,n]) for j=1:NUM_POINTS) for n=1:N); 
	return f;
end

function solve_local_opt(type_of_problem,time_limit)
	#returning the optimal value by conducting branching
	# type_of_problem is whether the problem is "Mixed" (for continuous mixed logit problem), or "MNL" (for multinomial logit problem),
	#time_limit is the time limit in seconds
	iteration_history=[-Inf Inf 0];
	break_points=[Float64.(LB_p)];
	push!(break_points,UB_p);
	Tolerr=0.0001;
	best_solution=[];
	bestf=-Inf;
	best_old=-Inf;
	A=[eye(NUM_POINTS);
		-eye(NUM_POINTS)];
	b=[LB_p;
		-UB_p];
	Diff=Inf;	
	fixed_diff=0;
	CompleteLoop=0;
	depth=1;
	# The algorithm do the split, then find new points to break, then find the cuts and do the splitting of the previous sets.  
	Num_breaks_prev=1;
	A_parent=[A];
	b_parent=[b];
	point_parents=[break_points];

	Upperbound_value=Inf;
	time_elapsed=0;
	while abs(bestf-Upperbound_value)>Tolerr &&time_elapsed<time_limit
		display("====================New_branch======================")
		iteration_history=[iteration_history; bestf Upperbound_value time_elapsed];
		CPUtic()
		best_old=copy(bestf);
		#split the feasible reagon
		m_break=size(break_points,1);
		A_leaf=Array{Any}(undef,m_break,Num_breaks_prev);
		b_leaf=Array{Any}(undef,m_break,Num_breaks_prev);
		LB_p=Array{Any}(undef,m_break,Num_breaks_prev);
		UB_p=Array{Any}(undef,m_break,Num_breaks_prev);
		points_leaf=Array{Any}(undef,m_break,Num_breaks_prev);

		Best_solution_leaf=Array{Any}(undef,m_break,Num_breaks_prev);
		Bestf_leaf=zeros(m_break,Num_breaks_prev);
		Best_Upperbound_value=zeros(m_break,Num_breaks_prev);
		status=zeros(m_break,Num_breaks_prev);# to track whether we have feasible region. all the regions are assumed to be infeasible at the beginning
		time_elapsed=time_elapsed+CPUtoq();
		for i=1:m_break
			CPUtic()
			Scenario_help=copy([break_points[i:i];break_points[1:i-1];break_points[i+1:m_break]]);
			A_lea,b_lea=cuts(Scenario_help);
			time_elapsed=time_elapsed+CPUtoq();
			for j=1:Num_breaks_prev
				CPUtic()
				A_leaf[i,j]=[A_parent[j];A_lea];b_leaf[i,j]=[b_parent[j];b_lea];
				#check feasibility of the area (at least one of the points should be in any of the region)
				num_feas_point=0;
				initial_point=zeros(NUM_POINTS,1);
				for mm=1:m_break
					if sum(A_leaf[i,j]*break_points[mm].>=b_leaf[i,j])==size(A_leaf[i,j],1)
						num_feas_point=num_feas_point+1;
						initial_point=break_points[mm];
					end
				end
				if num_feas_point>0
					status[i,j]=1;
				end
				time_elapsed=time_elapsed+CPUtoq();
				display("==Number of points $(m_break)===========Remaining node $((m_break-i)*(Num_breaks_prev)+Num_breaks_prev-j)========$(bestf)==========$(Upperbound_value) ===================$time_elapsed===")
				if status[i,j]==1 && time_elapsed<time_limit
					CPUtic()
					if type_of_problem=="MNL" 
						obj_fun=Distribution_function;
						obj_fun_i=Distribution_function_i;
					elseif type_of_problem=="Mixed" 
						obj_fun=Mixed_logit_function;
						obj_fun_i=Mixed_logit_function_i;
					end
					#finding the lower bound
					Bestf_leaf[i,j],Best_solution_leaf[i,j]=trust_region_iteration(A_leaf[i,j],b_leaf[i,j],Tolerr,obj_fun,initial_point,time_limit-time_elapsed)
					Best_solution_leaf[i,j]=abs.(Best_solution_leaf[i,j]); #for -0.0 to be positive
					time_elapsed=time_elapsed+CPUtoq();
					if time_elapsed<time_limit
						CPUtic()
						# founding the bounds needed for the overstimator
						LB_τ,UB_τ,LB_p[i,j],UB_p[i,j] = bound_identifier(A_leaf[i,j],b_leaf[i,j],obj_fun_i,N,time_limit-time_elapsed);
						#finding the upperbounds
						if LB_p!="solved"
							help=overestimator_general(A_leaf[i,j],b_leaf[i,j],break_points,obj_fun_i,LB_τ,UB_τ,LB_p[i,j],UB_p[i,j],N,time_limit-time_elapsed);
						else
							help="solved"
						end
						#updating the information
						if help=="solved"
							Best_Upperbound_value[i,j]=round.(Bestf_leaf[i,j],digits=4);
							status[i,j]=0;
						else
							Best_Upperbound_value[i,j]=help;
						end
						time_elapsed=time_elapsed+CPUtoq();
						
						if Bestf_leaf[i,j]>=bestf
							bestf=copy(Bestf_leaf[i,j]);
							best_solution=copy(Best_solution_leaf[i,j]);
							display("===========$(bestf)==========$(Upperbound_value) ===================$time_elapsed====")
							iteration_history=[iteration_history; bestf Upperbound_value time_elapsed];
						end
						#terminate branching if the upperbound is at most the same as the best found solution's objective function
						if Best_Upperbound_value[i,j]<=bestf
							status[i,j]=0;# we don't further branch this leaf as we can't get a better solution there
						end
						##

						

						#checking whether the found solution is among the existence points, if not put it in the set 
						CPUtic()
						m_break_help=size(break_points,1);
						point_exist=0;
						for m=1:m_break_help
							if sum(abs.(break_points[m]-Best_solution_leaf[i,j]))<Tolerr #check whether the solution is the same as before or not. if so mark it
								point_exist=1;
								break;
							end
						end
						# randomly choosing a point if the obtained solution is not new
						if point_exist==0 
							push!(break_points,Best_solution_leaf[i,j]);
						elseif help != "solved"
							max_i=1;max_j=1;max_width=0;
							for i_1=1:i
								for j_1=1:j
									if status[i_1,j_1]==1
										if maximum(UB_p[i_1,j_1]-LB_p[i_1,j_1])>max_width
											max_width=maximum(UB_p[i_1,j_1]-LB_p[i_1,j_1]);
											max_i=copy(i_1);
											max_j=copy(j_1);
										end
									end
								end
							end
							if status[max_i,max_j]==1
								p=randomPoint(A_leaf[max_i,max_j],b_leaf[max_i,max_j],Best_solution_leaf[i,j],LB_p[max_i,max_j],UB_p[max_i,max_j]);
								point_exist=0;
								for m=1:m_break_help
									if sum(abs.(break_points[m]-p))<Tolerr #check whether the solution is the same as before or not. if so mark it
										point_exist=1;
										break;
									end
								end
								if point_exist==0
									push!(break_points,round.(abs.(p),digits=4));
								end
							end
						end
						time_elapsed=time_elapsed+CPUtoq();
					end
				end
			end

		end

		#checking the tolerance found between the previous best and this iteration
		CPUtic()
		if best_old==-Inf
			Diff=bestf;
		else
			Diff=abs(best_old-bestf);
			if Diff<Tolerr
				fixed_diff=fixed_diff+1;
			end
			if bestf>best_old
				fixed_diff=0;
			end
		end
		# checking whether there is a new point added to the list, otherwise report that the points loop has been occurred
		m_break_help=size(break_points);

		if m_break_help==m_break #check whether a new point has been added. If not the loop will be terminated
			CompleteLoop==1;
		end
		# reconstruction of the parents: we go to the next branching
		ind_feasible=findall(status.==1);
		Num_breaks_prev=size(ind_feasible,1);
		if Num_breaks_prev>0
			A_parent=Array{Any}(undef,Num_breaks_prev,1);
			b_parent=Array{Any}(undef,Num_breaks_prev,1);
			for i=1:Num_breaks_prev
					A_parent[i]=A_leaf[ind_feasible[i]];
					b_parent[i]=b_leaf[ind_feasible[i]];
			end
			depth=depth+1;
			Upperbound_value=min(Upperbound_value,maximum(Best_Upperbound_value))
		else
			Upperbound_value=bestf;
		end

		time_elapsed=time_elapsed+CPUtoq();
		display("===========$(bestf)==========$(Upperbound_value) ===================$time_elapsed====")
		
	end
	return bestf,best_solution,Upperbound_value,iteration_history;

end

function trust_region_iteration(A,b,Tolerr,obj_fun::Function ,initial_point,time_limit)
	# do the iteration of the trust-region
	# the feasible region of Ap⩾b 
	
	function g(x)#gradient function
		return FiniteDiff.finite_difference_gradient(obj_fun, x);
	end

	NUM_POINTS=size(A,2);
	constraint_size=size(A,1);
    f0=0;
    initial_p=copy(initial_point);
    f0=obj_fun(initial_p);
    radius=0.1;
    time_elapsed=0;
  	f1=Inf;
    while abs(f0-f1)>Tolerr && norm(g(initial_p))>0.01*Tolerr &&time_limit>time_elapsed
		CPUtic()
		p1=trust_region(initial_p,obj_fun,A,b,radius,time_limit-time_elapsed)
		time_elapsed=time_elapsed+CPUtoq();
		f1=obj_fun(p1);
		while f1>f0
			CPUtic()
			initial_p=copy(p1);
			f0=copy(f1);
			p1=trust_region(initial_p,obj_fun,A,b,radius,time_limit-time_elapsed);
			f1=obj_fun(p1);
			radius=1;
			time_elapsed=time_elapsed+CPUtoq();
		end
		radius=radius/2;
    end
    return f0,initial_p;
end
function randomPoint(A,b,p_0,LB_p,UB_p)
	#if we find p_0 and it is not a new solution, we get to a new solution
	# Ax⫺b is the polytope with the lower bound vector LB_p and upperbound UB_p 
	ind=findall(UB_p-LB_p.==maximum(UB_p-LB_p));
	Right_side=1;
	if p_0[ind[1]]<=LB_p[ind[1]]+round((UB_p[ind[1]]-LB_p[ind[1]])/2,digits=4)
		Right_side=0;
	end
	rand_generator=0;
	if p_0[ind[1]]==LB_p[ind[1]]+round((UB_p[ind[1]]-LB_p[ind[1]])/2,digits=4)
		rand_generator=1;
	end
	NUM_POINTS=size(A,2);
	constraint_size=size(A,1);
	c=rand(NUM_POINTS,1)-0.5*ones(NUM_POINTS,1);
	rand_point = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_PERIND=0));
    @variable(rand_point,p[1:NUM_POINTS]);
    if Right_side==1
	    @constraint(rand_point,p[ind[1]]<=LB_p[ind[1]]+round((UB_p[ind[1]]-LB_p[ind[1]])/2,digits=4));
	else
		@constraint(rand_point,p[ind[1]]>=LB_p[ind[1]]+round((UB_p[ind[1]]-LB_p[ind[1]])/2,digits=4));
	end
    @constraint(rand_point,A*p.>=b)
    if rand_generator==0
	    @objective(rand_point,Min,sum((p-p_0)'*(p-p_0)));
	else
		@objective(rand_point,Max,sum( c'*p));
	end
    sol=solve(rand_point);
    if sol == :Optimal
	    initial_p=getvalue(p)[:];
	else
		rand_point = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_PERIND=0));
	    @variable(rand_point,p[1:NUM_POINTS]);
	    @constraint(rand_point,A*p.>=b)
	    @objective(rand_point,Max,sum( c'*p));
	    sol=solve(rand_point);
		initial_p=getvalue(p)[:];
	end
    return round.(initial_p,digits=4)
end

function overestimator_general(A,b,P_0,obj_fun,LB_τ,UB_τ,LB_p,UB_p,N,time_limit)

	#finding the upperbound of the problem
	#the feasible region is Ax⫺b
	#obj_fun is a function reprensting the objective fucntion
	#LB_τ and UB_τ are the lower and upper bounds of τ
	#LB_p and UB_p are the lower and upper bounds of τ
	# N is the number of customers
	#
	#
	if LB_τ=="solved" #because of the rounding errors we get infeasible. So, the area is small
		return "solved"
	end
	if sum((UB_τ-LB_τ))<0.0001 
		return "solved"
	end
	if sum((UB_p-LB_p))<0.0001 
		return "solved"
	end

	LB_f=1 ./UB_τ;
	UB_f=1 ./LB_τ;

	
	A=round.(A,digits=4);
	b=round.(b,digits=4);
	LB_τ=round.(LB_τ,digits=4);
	UB_τ=round.(UB_τ,digits=4);
	LB_p=round.(LB_p,digits=4);
	UB_p=round.(UB_p,digits=4);
	LB_f=round.(LB_f,digits=4);
	UB_f=round.(UB_f,digits=4);

	LP=Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_PERIND=0,CPX_PARAM_TILIM=time_limit));
	P_size=size(A,2);
	@variable(LP, W[1:P_size,1:P_size]>=0);
	@variable(LP, τ[1:P_size]>=0)
	@variable(LP, γ[1:P_size]>=0)
	violance=0;
	@variable(LP, p[1:P_size]>=0)
	@constraint(LP, A*p.>=b)
	for i=2:P_size
		@constraint(LP, τ[i]<=UB_τ[i]+violance);
		@constraint(LP, τ[i]>=LB_τ[i]-violance);
		@constraint(LP, γ[i]>=0);
		@constraint(LP, γ[i]<=1);
		for j=1:P_size
			@constraint(LP,W[i,j] >=LB_τ[i]*p[j]+τ[i]*LB_p[j]-LB_τ[i]*LB_p[j])
			@constraint(LP,W[i,j] >=UB_τ[i]*p[j]+τ[i]*UB_p[j]-UB_τ[i]*UB_p[j])
			@constraint(LP,W[i,j] <=UB_τ[i]*p[j]+τ[i]*LB_p[j]-UB_τ[i]*LB_p[j])
			@constraint(LP,W[i,j] <=τ[i]*UB_p[j]+LB_τ[i]*p[j]-LB_τ[i]*UB_p[j])
			@constraint(LP,W[i,j] -UB_τ[i]*p[j] <=violance)
			@constraint(LP,W[i,j] -LB_τ[i]*p[j] >=-violance)
			@constraint(LP, W[i,j] <=τ[i]*UB_p[j]+violance);
			@constraint(LP, W[i,j] >=τ[i]*LB_p[j]-violance);
		end
		if abs(LB_τ[i])>0.00001
			@constraint(LP, γ[i]<=τ[i]/LB_τ[i]);
		end
		if abs(UB_τ[i])>0.00001
			@constraint(LP, γ[i]>=τ[i]/UB_τ[i]);
		end
		
		for t=1:size(A,1)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-b[t]*τ[i]>=LB_τ[i]*(sum(A[t,j]*p[j] for j=1:P_size)-b[t])-violance)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-b[t]*τ[i]<=UB_τ[i]*(sum(A[t,j]*p[j] for j=1:P_size)-b[t])+violance)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-sum(A[t,j]*p[j]*LB_τ[i] for j=1:P_size)>=(τ[i]- LB_τ[i])*b[t]-violance)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-sum(A[t,j]*p[j]*UB_τ[i] for j=1:P_size)<=(τ[i]- UB_τ[i])*b[t]+violance)
		end
		Num_points=size(P_0,1)
		for k=1:Num_points
			p_0=round.(P_0[k],digits=4);
			function f_derivative(p)
				f=x->obj_fun(x,i);
				g=x->FiniteDiff.finite_difference_gradient(f, x);
				return g(p);
			end
			Matrix_f_Derivative=f_derivative(p_0);
			value_f_i=obj_fun(p_0,i);
			@constraint(LP, τ[i]*(value_f_i)+ sum(-τ[i]*(Matrix_f_Derivative[j]*p_0[j])+(W[i,j]*Matrix_f_Derivative[j]) for j = 1:P_size )<=γ[i]+violance)
			for j=1:P_size
				
				if value_f_i>0
					@constraint(LP, 4*W[i,j]+(p_0[j]-1/value_f_i)^2+2*(p_0[j]-1/value_f_i)*(p[j]-p_0[j])-2*(p_0[j]-1/value_f_i)*(τ[i]-(1/value_f_i)) <= (UB_p[j]+UB_τ[i])^2 ); #linearization of w+(p+τ)^2	⩽ (UB_p[j]+UB_τ[i,n])^2
				end
			end
			if LB_τ[i]>0
				@constraint(LP,γ[i]>=LB_τ[i]*(value_f_i+sum(Matrix_f_Derivative[j]*(p[j]-p_0[j]) for j=1:P_size))+τ[i]*LB_f[i]-LB_τ[i]*LB_f[i])#McCormick
				@constraint(LP,value_f_i+sum(Matrix_f_Derivative[j]*(p[j]-p_0[j]) for j=1:P_size)<=γ[i]/LB_τ[i])
			elseif LB_τ[i]<0
				@constraint(LP,γ[i]<=LB_τ[i]*(value_f_i+sum(Matrix_f_Derivative[j]*(p[j]-p_0[j]) for j=1:P_size))+τ[i]*UB_f[i]-LB_τ[i]*UB_f[i])#McCormick
				@constraint(LP,value_f_i+sum(Matrix_f_Derivative[j]*(p[j]-p_0[j]) for j=1:P_size)>=γ[i]/LB_τ[i])
			end
			if UB_τ[i]>0
				@constraint(LP,γ[i]>=UB_τ[i]*(value_f_i+sum(Matrix_f_Derivative[j]*(p[j]-p_0[j]) for j=1:P_size))+τ[i]*UB_f[i]-UB_τ[i]*UB_f[i])#McCormick
			elseif UB_τ[i]<0
				@constraint(LP,γ[i]<=UB_τ[i]*(value_f_i+sum(Matrix_f_Derivative[j]*(p[j]-p_0[j]) for j=1:P_size))+τ[i]*LB_f[i]-UB_τ[i]*LB_f[i])#McCormick
			end
		end
	end
	@objective(LP,Max,sum(W[i,i] for i=2:P_size))
	sol=solve(LP);
	if sol== :Infeasible || sol == :CPX_STAT_ABORT_DUAL_OBJ_LIM
		return "solved"# this is because the infeasibility comes from rounding errors
	end
	W=getvalue(W);
	return sum(W[i,i] for i=2:P_size)
end


function trust_region(P_k,obj_fun::Function,A,b,radius,time_limit)
	# Ap ⫺ b is the feasible set 
	# P_k is the obtained feasible solution in the previous round
	# radius shows the size of the ball
	P_size=size(A,2);
	Model_k=Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_PERIND=0,CPX_PARAM_TILIM=time_limit));
	@variable(Model_k,p[1:P_size])
	@constraint(Model_k,A*p.>=b)
	function g(x)
		return FiniteDiff.finite_difference_gradient(obj_fun, x);# gradient function
	end
	@objective(Model_k,Max,sum((g(P_k)).*(p)));
	@variable(Model_k, auxABS[1:P_size]);
	@constraint(Model_k, sum(auxABS)<=radius)
	for i=1:P_size
		@constraint(Model_k, auxABS[i]>=p[i]-P_k[i])
		@constraint(Model_k, -auxABS[i]<=p[i]-P_k[i])
	end
	solve(Model_k)
	P_new=getvalue(p);
	return P_new
end

function bound_identifier(A,b,obj_fun::Function,N,time_limit)
	#bound identifier for τ and p
	#Ap⫺b
	#obj_fun is the objective function
	#N is the number of customer
	A=round.(A,digits=4);
	b=round.(b,digits=4);
	P_size=size(A,2);
	LB_τ=zeros(P_size);
	UB_τ=zeros(P_size);
	LB_p=zeros(P_size);
	UB_p=zeros(P_size);
	p_init=zeros(P_size,1);
	for i=1:P_size
		LP=Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_PERIND=0));
		@variable(LP, p[1:P_size]>=0)
		@constraint(LP, A*p.>=b)
		@objective(LP,Max,p[i])
		sol=solve(LP)
		if sol==:Infeasible ||sol==:CPX_STAT_ABORT_DUAL_OBJ_LIM || sol== :CPX_STAT_NUM_BEST 
			return "solved","solved","solved","solved";
		end
		UB_p[i]=getvalue(p)[i];
		LP=Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_PERIND=0));
		@variable(LP, p[1:P_size]>=0)
		@constraint(LP, A*p.>=b)
		@objective(LP,Min,p[i])
		sol=solve(LP)
		p_init=getvalue(p);
		LB_p[i]=getvalue(p)[i];
	end
	#we know that τ_in = 1/f_in. So, to find the min we can minize 1/f_in, as it is convex, and to find the
	#max we can simply minimize f_in and then use the argmin as the argmax for τ. To solve the convex min we use ipopt
	for i=1:P_size
		function ff(x)
			return -1/obj_fun(x,i);
		end						
		f0,p=trust_region_iteration(A,b,0.00001, ff,p_init,time_limit)
		LB_τ[i]=1/obj_fun(p,i);
		function f(x)
			return -obj_fun(x,i);
		end
		f0,p=trust_region_iteration(A,b,0.00001,f,p_init,time_limit)
		UB_τ[i]=1/obj_fun(p,i);
	end
	return LB_τ,UB_τ,LB_p,UB_p;
end


function solve_nlopt()
	#solving the problem using NLopt
	opt = Opt(:GN_ESCH, P_size);
	opt.lower_bounds = LB_p;
	opt.upper_bounds = UB_p;
	opt.maxtime= 61200;
	opt.max_objective = Mixed_logit_function_nlopt
	t=@timed (minf,p,ret) = optimize(opt, [0;0;0])
	return t
end