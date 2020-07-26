# Load packages
using FileIO, CSV, DataFrames, Plots, DelimitedFiles, Distributions, StatsBase, LinearAlgebra, SparseArrays, JLD, JLD2, Random

# Load functions
include("VLStm.jl")

# Load data
df = CSV.read("data_final2.csv");
rename!(df,:HHNEWLAT => :latitude);
rename!(df,:HHNEWLNG => :longitude);
# Select which paras to include
para = 1:4;
str = "_Paras1to4"

# Subset data
subs = in.(df.PARA,[para]);
df = df[subs,:];

# # Plot longitude and latitude coordinates
# plot(df.longitude,df.latitude)

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
w = 0;

# Set PKDL parameters
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

# Run simulations with different proportions of asymptomatic infection
# Set random seed
Random.seed!(12345)

for i=1:3
    # Load MCMC output for model parameters and missing data
    p = readdlm("p" * str * "_" * string(i) * ".csv",',',Float64);
    p1 = readdlm("p1" * str * "_" * string(i) * ".csv",',',Float64);
    tAs = readdlm("tAs" * str * "_" * string(i) * ".csv",',',Int);
    tRAs = readdlm("tRAs" * str * "_" * string(i) * ".csv",',',Int);
    tEs = readdlm("tEs" * str * "_" * string(i) * ".csv",',',Int);
    tRsANONR = try
        readdlm("tRsANONR" * str * "_" * string(i) * ".csv",',',Int);
    catch
        Matrix(undef,0,size(p,1));
    end
    tRsAONR = try
        readdlm("tRsAONR" * str * "_" * string(i) * ".csv",',',Int);
    catch
        Matrix(undef,0,size(p,1));
    end
    incA,incI,incP,origin,tmax,K02_10,P02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv = RunSims(df,para,startyr,startmo,endyr,endmo,w,p,p1,tAs,tRAs,tEs,tRsANONR,tRsAONR,nsmpls,nsims,r5,p5,pP,pD,hs,ps)
    SaveOutput(incA,incI,incP,K02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv,df,w,tAs,tRAs,tEs,origin,str * "_" * string(i))
end

# Run simulations for testing identifiability of additional within-household transmission with DIC
# Use values from MCMC run with symptomatic proportion of 0.15 for other parameters and missing data
p = readdlm("p" * str * "_2.csv",',',Float64)
p1 = readdlm("p1" * str * "_2.csv",',',Float64)
tAs = readdlm("tAs" * str * "_2.csv",',',Int);
tRAs = readdlm("tRAs" * str * "_2.csv",',',Int);
tEs = readdlm("tEs" * str * "_2.csv",',',Int);
tRsANONR = try
    readdlm("tRsANONR" * str * "_2.csv",',',Int);
catch
    Matrix(undef,0,size(p,1));
end
tRsAONR = try
    readdlm("tRsAONR" * str * "_2.csv",',',Int);
catch
    Matrix(undef,0,size(p,1));
end

# Run simulations with and without additional within-household transmission
# Set random seed
Random.seed!(123456)

# Set additional within-household transmission
deltas = [0,p[4]];

for i=1:length(deltas)
    p_i = hcat(p[:,1:3],deltas[i],p[:,5:7]);
    incA,incI,incP,origin,tmax,K02_10,P02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv = RunSims(df,para,startyr,startmo,endyr,endmo,w,p_i,p1,tAs,tRAs,tEs,tRsANONR,tRsAONR,nsmpls,nsims,r5,p5,pP,pD,hs,ps)
    SaveOutput(incA,incI,incP,K02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv,df,w,tAs,tRAs,tEs,origin,str * "_addtnl_wthnHH_trnsmssn_" * string(i))
end
