# Julia 1.5.3
#This document contains data used in the paper
# A. Marandi, V Lurkin (2020), Static Pricing Problems under Mixed Multinomial Logit Demand
#written by Ahmadreza Marandi
# All rights are reserved.
#
#
using MAT

function Logit_10()
	NUM_POINTS=2+1; #alternatives
	N=10; #customers

	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0;
	];
	 Age_veh =[	0
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		0];
	 Low_inc =[	1
	 		1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		1
	 		1
	 		0];
	 Res =[ 1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		1];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=		   Beta_AT * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + Beta_AT * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + Beta_AT * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=Beta_FEE + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=Beta_FEE + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	#				 
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit_10.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
 	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end
function Mixed_Logit_10(β)
	NUM_POINTS=2+1; #alternatives
	N=10; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0;
	];
	 Age_veh =[	0
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		0];
	 Low_inc =[	1
	 		1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		1
	 		1
	 		0];
	 Res =[ 1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		1];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=		   β[1] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + β[1] * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + β[1] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=β[2] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=β[2] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	#				 
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit_10.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
 	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end
function Mixed_Logit_50(β)
	NUM_POINTS=2+1; #alternatives
	N=50; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	 Age_veh =[1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	 Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	 Res =[ 1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=		   β[1] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + β[1] * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + β[1] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=β[2] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=β[2] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	#				 
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit_10.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
 	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end
function Logit_50()
	NUM_POINTS=2+1; #alternatives
	N=50; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	 Age_veh =[1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	 Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	 Res =[ 1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=Beta_AT * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + Beta_AT * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + Beta_AT * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=Beta_FEE + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=Beta_FEE + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end

function Mixed_Logit_n10_random(R_AT)
	NUM_POINTS=2+1; #alternatives
	N=10; #customers
	R=R_AT;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	 ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	source=@__DIR__;
	address=source*"\\random beta\\Beta_N"*string(N)* "R_"*string(R_AT)*".mat";
	input=matread(address)
	Beta_AT=input["Beta_AT"];#AT coefficient (10 customers x R_AT draws)
	 Beta_FEE=input["Beta_FEE"];#FEE coefficient (10 customers x R_FEE draws)
	
	#Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[
	 	1	0
	 	2	1
	 	3	1
	 	4	0
	 	5	0
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	0];
	  Origin= Origin[:,2:end]	;
	  Age_veh =[
	  	1	0
	 	2	0
	 	3	0
	 	4	1
	 	5	0
	 	6	0
	 	7	1
	 	8	0
	 	9	0
	 	10	0];
	  Age_veh= Age_veh[:,2:end];
	  Low_inc =[
	  	1	1
	 	2	1
	 	3	1
	 	4	1
	 	5	0
	 	6	1
	 	7	1
	 	8	1
	 	9	1
	 	10	0];
	  Low_inc= Low_inc[:,2:end];
	  Res =[
	 	1	1
	 	2	1
	 	3	1
	 	4	0
	 	5	1
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	1];
	   Res= Res[:,2:end];

	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r]=0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 	# 	 file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;

end

