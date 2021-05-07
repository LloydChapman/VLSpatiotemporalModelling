# Define plotting functions for simulation output
function PlotOvrlEpiCurve(startyr,t,inc,hold_on,clr,ylbl,plot_obs_inc,inc_obs,str)
    yr = startyr .+ (t[2:end].-1)/12;

    if hold_on
        plot!(yr,sum(inc[:,rand(1:size(inc,2),100),:],dims=3)[:,:,1],linecolor=clr)
    else
        plot(yr,sum(inc[:,rand(1:size(inc,2),100),:],dims=3)[:,:,1],linecolor=clr,legend=false,grid=false,xlabel="Year",ylabel=ylbl)
    end
    # plot!(yr,mean(sum(incI,dims=3),dims=2)[:,:,1])
    # plot!(yr,median(sum(incI,dims=3),dims=2)[:,:,1])
    plot!(yr,minimum(sum(inc,dims=3),dims=2)[:,:,1],linecolor=clr,linestyle=:dash,lw=2)
    plot!(yr,maximum(sum(inc,dims=3),dims=2)[:,:,1],linecolor=clr,linestyle=:dash,lw=2)

    if plot_obs_inc
        plot!(startyr .+ ((1:t[end]).-1)/12,inc_obs,linecolor=:black,lw=2)
        png("sim_vs_obs_VL_cases_total" * str)
    end
end

function PlotParaEpiCurve(startyr,t,inc,para,hold_on,clr,ylbl,plot_obs_inc,inc_obs,str)
    yr = startyr .+ (t[2:end].-1)/12;
    if hold_on
        plot!(yr,inc[:,rand(1:size(inc,2),100),para],linecolor=clr)
    else
        plot(yr,inc[:,rand(1:size(inc,2),100),para],linecolor=clr,legend=false,grid=false,xlabel="Year",ylabel=ylbl,title="Para " * string(para))
    end
    # plot!(t[2:end],mean(incI[:,:,k],dims=2),linecolor=clr)
    plot!(yr,minimum(inc[:,:,para],dims=2),linecolor=clr,linestyle=:dash,lw=2)
    plot!(yr,maximum(inc[:,:,para],dims=2),linecolor=clr,linestyle=:dash,lw=2)
    if plot_obs_inc
        plot!(startyr .+ ((1:t[end]).-1)/12,inc_obs,linecolor=:black,lw=2)
        # png("sim_vs_obs_VL_cases_para" * string(para) * str)
        savefig("sim_vs_obs_VL_cases_para" * string(para) * str * ".pdf")
    end
end

function CalcPercDiffInc(incI,incI1,probs)
    diff = convert(Array{Int16},incI1) - convert(Array{Int16},incI);
    tot_diff = sum(diff,dims=(1,3))[1,:,1];
    display(plot(tot_diff,seriestype=:histogram))
    mean_diff = mean(tot_diff);
    q_diff = quantile(tot_diff,probs);
    tot_incI = sum(incI,dims=(1,3))[1,:,1];
    perc_diff = 100*tot_diff./tot_incI;
    mean_perc_diff = mean(perc_diff);
    q_perc_diff = quantile(perc_diff,probs);
    return(mean_diff,q_diff,mean_perc_diff,q_perc_diff)
end

function CalcPercDiffParaInc(incI,incI1,probs)
    diff = convert(Array{Int16},incI1) - convert(Array{Int16},incI);
    tot_diff = sum(diff,dims=1)[1,:,:];
    mean_diff = mean(tot_diff,dims=1)[1,:];
    q_diff = mapslices(x->quantile(x,probs),tot_diff,dims=1);
    tot_inc = sum(incI,dims=1)[1,:,:];
    perc_diff = 100*tot_diff./tot_inc;
    mean_perc_diff = mapslices(x->mean(x[.!isnan.(x) .& .!isinf.(x)]),perc_diff,dims=1);
    q_perc_diff = mapslices(x->quantile(x[.!isnan.(x) .& .!isinf.(x)],probs),perc_diff,dims=1);
    return(mean_diff,q_diff,mean_perc_diff,q_perc_diff)
end

function CalcParaIncSumStats(inc,para,para_inc_obs,probs)
    para_tot_inc = sum(inc,dims=1)[1,:,:];
    mean_para_tot_inc = mean(para_tot_inc,dims=1)[1,:];
    median_para_tot_inc = median(para_tot_inc,dims=1)[1,:];
    q_para_tot_inc = mapslices(x -> quantile(x,probs),para_tot_inc,dims=1);
    tot_inc = sum(inc,dims=(1,3))[1,:,1];
    mean_tot_inc = mean(tot_inc);
    median_tot_inc = median(tot_inc);
    q_tot_inc = quantile(tot_inc,probs);
    df_inc = DataFrame([para,para_inc_obs,mean_para_tot_inc,median_para_tot_inc,q_para_tot_inc[1,:],q_para_tot_inc[2,:]],[:Para,:VLcases,:Mean,:Median,:LQ95,:UQ95]);
    push!(df_inc,[0,sum(para_inc_obs),mean_tot_inc,median_tot_inc,q_tot_inc[1],q_tot_inc[2]])
    return(df_inc)
end

function PrependIncObsw(inc,w,para_inc_obs)
    inc=vcat(permutedims(repeat(repeat(para_inc_obs[1:w,:],1,1,1),1,1,size(inc,2)),[1,3,2]),inc);
    return(inc)
end
