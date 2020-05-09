# Load packages
using FileIO, CSV, DataFrames, Plots, DelimitedFiles, Distributions, StatsBase, LinearAlgebra, SparseArrays, JLD, JLD2, Random

# Load functions
include("VLStm.jl")

function SaveOutput(incA,incI,incP,K02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv,df,w,tAs,tRAs,tEs,origin,str)
    # Reshape status matrix
    Status_wide = reshape(Status,length(tA),length(t));

    # Add observed data in unsimulated window
    tA[vec(tAs.<=w)] = tAs[tAs.<=w];
    tR[vec(tRAs.<=w)] = tRAs[tRAs.<=w];
    tEs1 = Vector{Union{Missing,Int64}}(missing,size(df,1));
    tEs1[K02_10] = tEs;
    tE[convertToBool(tEs1.<=w)] = tEs[tEs.<=w];
    infw = ((df.KA.>origin) .& (df.KA.<=origin+w) .& ismissing.(df.MIG_IN)) .& .!(df.INTMIG_IN.==1);
    infw = convertToBool(infw);
    tI[infw] = df.KA[infw].-origin;
    P02_10sim = Tuple(any(Status_wide.==6,dims=2));
    Dorw = ((df.KARX.>origin) .& (df.KARX.<=origin+w) .& ismissing.(df.MIG_IN)) .& .!(df.INTMIG_IN.==1) .& P02_10sim;
    Dorw = convertToBool(Dorw);
    tDor[Dorw] = df.KARX[Dorw].-origin;
    Recw = ((df.KARX.>origin) .& (df.KARX.<=origin+w) .& ismissing.(df.MIG_IN)) .& .!(df.INTMIG_IN.==1) .& .!P02_10sim;
    Recw = convertToBool(Recw)
    tR[Recw] = df.KARX[Recw].-origin;

    # Write event times into dataframe
    dfo = DataFrame([tA,tE,tI,tDor,tP,tR],[:tA,:tE,:tI,:tDor,:tP,:tR]);

    # Save simulation output
    CSV.write("sim_output_para1" * str * ".csv",dfo)
    writedlm("status_mtrx_para1" * str * ".csv",Status_wide,",")
    writedlm("init_stat_para1" * str * ".csv",Status_wide[:,1],",")
    writedlm("PKDL_infctsnss_para1" * str * ".csv",hv,",")

    save("incA_para1" * str * ".jld","incA",incA)
    save("incI_para1" * str * ".jld","incI",incI)
    save("incP_para1" * str * ".jld","incP",incP)

end

# Load data
df = CSV.read("data_final2.csv");
rename!(df,:HHNEWLAT => :latitude);
rename!(df,:HHNEWLNG => :longitude);
# Select which paras to include
para = 1;

# Subset data
subs = in.(df.PARA,[para]);
df = df[subs,:];

# # Plot longitude and latitude coordinates
# plot(df.longitude,df.latitude)

# Load MCMC output for model parameters and missing data
p_str = "p_Para1_fixed_pI_lambda0.csv";
p1_str = "p1_Para1_fixed_pI_lambda0.csv";
tAs_str = "tAs_Para1_fixed_pI_lambda0.csv";
tRAs_str = "tRAs_Para1_fixed_pI_lambda0.csv";
tEs_str = "tEs_Para1_fixed_pI_lambda0.csv";
tRsANONR_str = "tRsANONR_Para1_fixed_pI_lambda0.csv";
tRsAONR_str = "tRsAONR_Para1_fixed_pI_lambda0.csv";
p = readdlm(p_str,',',Float64);
p1 = readdlm(p1_str,',',Float64);
tAs = readdlm(tAs_str,',',Int);
tRAs = readdlm(tRAs_str,',',Int);
tEs = readdlm(tEs_str,',',Int);
tRsANONR = try
    readdlm(tRsANONR_str,',',Int);
catch
    Matrix(undef,0,size(p,1));
end
tRsAONR = try
    readdlm(tRsAONR_str,',',Int);
catch
    Matrix(undef,0,size(p,1));
end

# No. of posterior samples to run simulations for
nsmpls = 1;
# No. simulations for each sample
nsims = 1;

# Start and end years and months for the data
startyr = 2002; #2000; #2001; #1998; #2003;
startmo = 1;
endyr = 2010; #2018;
endmo = 12;
# Window not simulated at start of data (in months)
w = 12;

# Run simulations with PKDL
# MATLAB code: pars5=nbinfit(tRP(~isnan(tRP)&~PothrObs)-tP(~isnan(tRP)&~PothrObs)-1)
r5 = 1.183947478584665;
p5 = 0.065759912826616;
pP = 0.17; # proportion of VL cases who develop PKDL
pI0 = 0.15; # proportion of infected individuals who develop VL
pD = 16/((1-pI0)/pI0*1018); # estimate proportion of asx individuals who develop PKDL as ratio of asx incidence to incidence of PKDL w/o prior KA
h1 = 9/26/(10/15); # infectiousness of macular/papular PKDL cases
h3 = 18/21/(10/15); # infectiousness of nodular PKDL cases
h2 = (h1+h3)/2; # infectiousness of plaque PKDL cases
hs = [h1;h2;h3];
ps = [101/138;31/138;6/138]; # proportions of different PKDL types (macular/papular, plaque, nodular)
hu = sum(ps.*hs); # average infectiousness of unexamined PKDL cases
hs = [h1;h2;h3;hu];
ps = [101/190;31/190;6/190;52/190]; # proportions of different PKDL types including unexamined cases

# Set random seed
Random.seed!(1234)

# Run simulations with different proportions of asymptomatic infection
pIs = [0.1,0.15,0.2]
for i=1:length(pIs)
    p_i = hcat(p[:,1:6],pIs[i]);
    incA,incI,incP,origin,tmax,K02_10,P02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv = RunSims(df,para,startyr,startmo,endyr,endmo,w,p_i,p1,tAs,tRAs,tEs,tRsANONR,tRsAONR,nsmpls,nsims,r5,p5,pP,pD,hs,ps)
    SaveOutput(incA,incI,incP,K02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv,df,w,tAs,tRAs,tEs,origin,"_" * string(i))
end

# # Save variables needed for plotting - have to use JLD2 due to issues with saving Booleans with JLD
# save("plot_vrbles_w12.jld2","para",para,"startyr",startyr,"w",w,"t",t,"origin",origin,"tmax",tmax,"K02_10",K02_10,"P02_10",P02_10);
