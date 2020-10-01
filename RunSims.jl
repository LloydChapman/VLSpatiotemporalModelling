# Load packages
using FileIO, CSV, DataFrames, Plots, DelimitedFiles, Distributions, StatsBase, LinearAlgebra, SparseArrays, JLD, JLD2

# Load functions
include("VLStm.jl")

# Load data
df = CSV.read("data_final2.csv");
rename!(df,:HHNEWLAT => :latitude);
rename!(df,:HHNEWLNG => :longitude);
# Select which paras to include
para = 1:19; #1; #[1:6;16:19]; #[2:3;5;8:11;13:16;18]; #7:15; #1:4; #1:3; #1:2;

# Subset data
subs = in.(df.PARA,[para]);
df = df[subs,:];

# # Plot longitude and latitude coordinates
# plot(df.longitude,df.latitude)

# Load MCMC output for model parameters and missing data
p_str = "p_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas.csv";
p1_str = "p1_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas.csv";
tAs_str = "tAs_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas.csv";
tRAs_str = "tRAs_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas.csv";
tEs_str = "tEs_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas.csv";
tRsANONR_str = "tRsANONR_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas.csv";
tRsAONR_str = "tRsAONR_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas.csv";
p = readdlm(p_str,',',Float64);
p1 = readdlm(p1_str,',',Float64);
tAs = readdlm(tAs_str,',',Int);
tRAs = readdlm(tRAs_str,',',Int);
tEs = readdlm(tEs_str,',',Int);
tRsANONR = readdlm(tRsANONR_str,',',Int);
tRsAONR = readdlm(tRsAONR_str,',',Int);


# No. of posterior samples to run simulations for
nsmpls = 100; #10; #
# No. simulations for each sample
nsims = 100; #10; #

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
incA, incI, incP, origin, tmax, K02_10, P02_10, t = RunSims(df,para,startyr,startmo,endyr,endmo,w,p,p1,tAs,tRAs,tEs,tRsANONR,tRsAONR,nsmpls,nsims,r5,p5,pP,pD,hs,ps)
# Save variables needed for plotting - have to use JLD2 due to issues with saving Booleans with JLD
save("plot_vrbles_w12.jld2","para",para,"startyr",startyr,"w",w,"t",t,"origin",origin,"tmax",tmax,"K02_10",K02_10,"P02_10",P02_10);
# Save output
save("incA_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_w12.jld","incA",incA)
save("incI_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_w12.jld","incI",incI)
save("incP_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_w12.jld","incP",incP)

# Run simulations with PKDL uninfectious
incA1, incI1, incP1 = RunSims(df,para,startyr,startmo,endyr,endmo,w,p,p1,tAs,tRAs,tEs,tRsANONR,tRsAONR,nsmpls,nsims,r5,p5,pP,pD,zeros(4),ps)
# Save output
save("incA_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_noPKDL_w12.jld","incA",incA1)
save("incI_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_noPKDL_w12.jld","incI",incI1)
save("incP_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_noPKDL_w12.jld","incP",incP1)

# Run simulations with mean PKDL duration of infectiousness halved
r5_1 = r5;
p5_1 = 2*r5_1/((r5*(1/p5-1)) - 1 + 2*r5_1);
incA2, incI2, incP2 = RunSims(df,para,startyr,startmo,endyr,endmo,w,p,p1,tAs,tRAs,tEs,tRsANONR,tRsAONR,nsmpls,nsims,r5_1,p5_1,pP,pD,hs,ps)
# Save output
save("incA_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_shrtrPKDLdur_w12.jld","incA",incA2)
save("incI_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_shrtrPKDLdur_w12.jld","incI",incI2)
save("incP_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_shrtrPKDLdur_w12.jld","incP",incP2)
