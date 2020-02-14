# Load packages
using FileIO, CSV, DataFrames, Plots, DelimitedFiles, Distributions, StatsBase, LinearAlgebra, SparseArrays, JLD, JLD2

# Load plotting functions
include("PlotFns.jl")

# Load data
df = CSV.read("data_final2.csv");
rename!(df,:HHNEWLAT => :latitude);
rename!(df,:HHNEWLNG => :longitude);

# Load variables required for plotting
para = load("plot_vrbles_w12.jld2")["para"];
startyr = load("plot_vrbles_w12.jld2")["startyr"];
w = load("plot_vrbles_w12.jld2")["w"];
t = load("plot_vrbles_w12.jld2")["t"];
origin = load("plot_vrbles_w12.jld2")["origin"];
tmax = load("plot_vrbles_w12.jld2")["tmax"];
K02_10 = load("plot_vrbles_w12.jld2")["K02_10"];
P02_10 = load("plot_vrbles_w12.jld2")["P02_10"];

# Load simulated data
# simulations with PKDL
incI = load("incI_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_w12.jld")["incI"];
incA = load("incA_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_w12.jld")["incA"];
incP = load("incP_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_w12.jld")["incP"];

# simulations without PKDL
incI1 = load("incI_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_noPKDL_w12.jld")["incI"];
incA1 = load("incA_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_noPKDL_w12.jld")["incA"];
incP1 = load("incP_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_noPKDL_w12.jld")["incP"];

# simulations with mean PKDL duration halved
incI2 = load("incI_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_shrtrPKDLdur_w12.jld")["incI"];
incA2 = load("incA_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_shrtrPKDLdur_w12.jld")["incA"];
incP2 = load("incP_DELTA_HGHR_PRESX_ASX_INFCTS_AllParas_shrtrPKDLdur_w12.jld")["incP"];

# Plot simulation output
# Overall epi curves
incIobs = counts(df.KA[K02_10 .& .!ismissing.(df.KA)].-origin,1:tmax);
PlotOvrlEpiCurve(startyr,t,incI,false,:lightgray,"No. VL cases",false,[],"")
PlotOvrlEpiCurve(startyr,t,incI1,true,:lightsalmon,"",true,incIobs,"_diffintvntns_w12")
PlotOvrlEpiCurve(startyr,t,incI2,false,:lightblue,"No. VL cases",true,incIobs,"_shrtrPKDLdur_w12")

# incPobs = counts(df.PKDL[P02_10].-origin,w+1:tmax);

# Para-level epi curves
para_incIobs = zeros(tmax,length(para));
for k=1:length(para)
    i = para[k];
    para_incIobs[:,i] = counts(df.KA[K02_10 .& .!ismissing.(df.KA) .& (df.PARA.==i)].-origin,1:tmax);
    PlotParaEpiCurve(startyr,t,incI,i,false,:lightgray,"No. VL cases",false,[],"")
    PlotParaEpiCurve(startyr,t,incI1,i,true,:lightsalmon,"",true,para_incIobs[:,i],"_diffintvntns_w12")
    # PlotParaLvlEpiCurve(startyr,t,incI1,i,true,:lightsalmon,"",false,[],"")
    # PlotParaLvlEpiCurve(startyr,t,incI2,i,true,:lightblue,"",true,para_incIobs[:,i],"_diffintvntns")
end

# Add incidence for unsimulated window at start of data
incI = PrependIncObsw(incI,w,para_incIobs);
incI1 = PrependIncObsw(incI1,w,para_incIobs);
incI2 = PrependIncObsw(incI2,w,para_incIobs);

# Calculate impact of different PKDL interventions
# Overall
probs = [0.025, 0.975];
mean_diff1, q_diff1, mean_perc_diff1, q_perc_diff1 = CalcPercDiffInc(incI,incI1,probs)
mean_diff2, q_diff2, mean_perc_diff2, q_perc_diff2 = CalcPercDiffInc(incI,incI2,probs)

# Para-level
para_tot_incIobs = sum(para_incIobs,dims=1)[1,:];
df_incI = CalcParaIncSumStats(incI,para,para_tot_incIobs,probs);
df_incI1 = CalcParaIncSumStats(incI1,para,para_tot_incIobs,probs);
df_incI2 = CalcParaIncSumStats(incI2,para,para_tot_incIobs,probs);

para_mean_diff1, para_q_diff1, para_mean_perc_diff1, para_q_perc_diff1 = CalcPercDiffParaInc(incI,incI1,probs);
para_mean_diff2, para_q_diff2, para_mean_perc_diff2, para_q_perc_diff2 = CalcPercDiffParaInc(incI,incI2,probs);

# Construct dataframe of simulation results for all paras
df_inc = join(df_incI,select(df_incI1,Not(:VLcases)),on=:Para,makeunique=true);
df_inc.MeanDiff1 = [para_mean_diff1;mean_diff1];
df_inc.LQ95Diff1 = [para_q_diff1[1,:];q_diff1[1]];
df_inc.UQ95Diff1 = [para_q_diff1[2,:];q_diff1[2]];
df_inc = join(df_inc,select(df_incI2,Not(:VLcases)),on=:Para,makeunique=true);
df_inc.MeanDiff2 = [para_mean_diff2;mean_diff2];
df_inc.LQ95Diff2 = [para_q_diff2[1,:];q_diff2[1]];
df_inc.UQ95Diff2 = [para_q_diff2[2,:];q_diff2[2]];

# Save summary dataframe
CSV.write("sim_vs_obs_VL_cases_w12.csv",df_inc)
