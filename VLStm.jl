# Define functions
function CalcHHDists(df)
    # Set approximate radius of Earth in metres for use in Haversine formula
    r = 6371.0088*1000;

    # Find indices of unique HHs
    ia = findfirst.(isequal.(unique(df.HHID)),[df.HHID]);
    ib = reduce(vcat,findall.(isequal.(df.HHID),[unique(df.HHID)]));

    # Convert HH latitudes and longitudes from degrees to radians
    lat = df.latitude[ia]/180*pi;
    lon = df.longitude[ia]/180*pi;

    lat_diff = abs.(lat.-lat');
    lon_diff = abs.(lon.-lon');

    dHH = r*2*asin.(sqrt.((sin.(lat_diff/2)).^2 .+ cos.(lat)*cos.(lat').*(sin.(lon_diff/2)).^2));

    return dHH, ia, ib

end

function KnlHH(dHH,dHHsqrd,alpha,typ,n,nHH,f)

    if typ == "Cauchy"
        KHH = 1 ./(1 .+dHHsqrd/alpha^2);
    elseif typ == "Exp"
        KHH = exp.(-dHH/alpha);
    elseif typ == "Const"
        KHH = ones(nHH,nHH);
    end

    K0 = n/sum(f*f'.*KHH); # calculate normalisation constant
    KHH = K0*KHH;

    return KHH, K0

end

function stata_month(y,m)
    M=(y-1960)*12+m-1
    return(M)
end

function ProbInitStatus(a,lambda,p2)
    probs = [exp.(-lambda*a) (1-exp(-lambda))*((1-p2).^a-exp.(-lambda*a))/(1-p2-exp(-lambda))];
    return(probs)
end

function convertToBool(x)
    x[ismissing.(x)] .= false;
    x = convert(Array{Bool},x);
    return x
end

function randSizeBiasedNegativeBinomial(r,p,n)
    u = rand(n,1);
    mu = r*(1-p)/p+1;
    M = repeat(transpose(1:100),n,1);
    x = sum(cumsum((1 .- cdf.(NegativeBinomial(r,p),M.-2)),dims=2)/mu .<= repeat(u,1,size(M,2)), dims=2) .+ 1;
    return(x)
end

function Iterate(n,Status,rateHH,epsilon,ib,h4,h0,hs,ps,hv,t,tB,tA,tE,tI,tDor,tP,tR,tD,tIM,tEM,p2,r1,p1,r0,p0,r3,p3,r5,p5,pI,pP,pD,infDor,infRec,AsxDor,AsxRec,INTMIG_OUT,INTMIG_IN,IM_OUT,IM_IN)

    Event = zeros(n);
    # Time = Inf*ones(n);
    # TimeD = Inf*ones(n);

    PreB = (Status.==0);
    Sus = (Status.==1) .& (tD .> t) .& (tEM .> t);
    Asx = (Status.==2) .& (tD .> t) .& (tEM .> t);
    Exp = (Status.==3) .& (tD .> t) .& (tEM .> t);
    inf = (Status.==4) .& (tD .> t) .& (tEM .> t);
    Dor = (Status.==5) .& (tD .> t) .& (tEM .> t);
    PKDL = (Status.==6) .& (tD .> t) .& (tEM .> t);
    Rec = (Status.==7) .& (tD .> t) .& (tEM .> t);
    Death = (Status.==8);
    PreDeath = .!PreB .& .!Death;
    Rate = zeros(n);
    # Rate[Sus] = rate[Sus,inf]*ones(sum(inf),1).+epsilon;
    Rate[Sus] = rateHH[ib[Sus],[ib[Asx];ib[Exp];ib[inf];ib[PKDL]]]*[h4*ones(sum(Asx),1);h0*ones(sum(Exp),1);ones(sum(inf),1);hv[PKDL]].+epsilon;

    # # Direct Gillespie algorithm
    # Rate[Exp] .= sigma;
    # Rate[inf] .= gamma;
    # TotalRate = sum(Rate);
    # CumRate = cumsum(Rate);
    # i = findfirst(TotalRate*rand().<CumRate);
    # step = -log(rand())/TotalRate;

    # Gillespie first reaction method
    # Time[PreB] = tB[PreB].-t;
    # # Time[Sus] = -log.(rand(sum(Sus)))./Rate[Sus];
    # Time[Asx] = tR[Asx].-t;
    # Time[Exp] = tI[Exp].-t;
    # Time[inf.&infRec] = tR[inf.&infRec].-t;
    # Time[inf.&infDor] = tDor[inf.&infDor].-t;
    # Time[Dor] = tP[Dor].-t;
    # Time[PKDL] = tR[PKDL].-t;
    # # Time[Asx.|(inf.&infRec).|PKDL] = Time[Asx.|(inf.&infRec).|PKDL].-t;
    # TimeD[PreDeath] = tD[PreDeath].-t;
    # Times = [Time TimeD];
    # step = minimum(Times[Times.>0]);
    # indcs = findall(isequal(step),Times);

    prob = 1 .- exp.(-Rate);

    # Update time
    # t = t+step;
    t = t+1;
    # println(step)
    # println(t)
    indcs = findall((rand(n).<prob) .| (tB.==t) .| (tI.==t) .| (tDor.==t) .| (tP.==t) .| (tR.==t) .| (tD.==t) .| (tIM.==t) .| (tEM.==t));
    # Loop over events (multiple births and deaths occur at the same time as they're only recorded to the nearest month)
    for k=1:length(indcs)
        # Get index of individual for whom event occurs (i)
        i = indcs[k];
        if (tD[i] == t) | (tEM[i] == t) # death or emigration
            Event[i] = 8 - Status[i];
        elseif tIM[i] == t # migration in
            if INTMIG_IN[i] & any(IM_IN.==i) # internal migration in
                # # Make vector of event times
                # event_times = [tA[i],tE[i],tI[i],tDor[i],tP[i],tR[i]];
                # if any(event_times.<tIM[i]) # any event is before migration
                #     # Find last event that happens before migration
                #     Event[i] = findlast(event_times.<tIM[i]) + 1;
                # elseif Status[IM_OUT[IM_IN.==i][1]]!=1
                    Event[i] = Status[IM_OUT[IM_IN.==i][1]] - Status[i];
                # else # no events have happened so still susceptible
                #     Event[i] = 1;
                # end
            else # external migration in
                Event[i] = 1;
            end
        elseif Sus[i]
            if rand()<pI
                tE[i] = t;
                tI[i] = t+rand(NegativeBinomial(r1,p1))+1;
                # I think there's an issue here, as coded like this their PKDL
                # onset could be after their death. Simulate birth and death
                # times? Otherwise I'm assuming death rate is unaffected by KA.
                if rand()<pP
                    tDor[i] = tI[i]+rand(NegativeBinomial(r0,p0))+1;
                    tP[i] = tDor[i]+rand(NegativeBinomial(r3,p3));
                    tR[i] = tP[i]+rand(NegativeBinomial(r5,p5))+1;
                    hv[i] = hs[rand(Categorical(ps))];
                    infDor[i] = true;
                else
                    tR[i] = tI[i]+rand(NegativeBinomial(r0,p0))+1;
                    infRec[i] = true;
                    # Add relapse?
                end
                Event[i] = 2;
            else
                tA[i] = t;
                if rand()<pD
                    tDor[i] = t+rand(Geometric(p2))+1;
                    tP[i] = tDor[i]+rand(NegativeBinomial(r3,p3));
                    tR[i] = tP[i]+rand(NegativeBinomial(r5,p5))+1;
                    hv[i] = hs[rand(Categorical(ps))];
                    AsxDor[i] = true;
                else
                    tR[i] = t+rand(Geometric(p2))+1;
                    AsxRec[i] = true;
                end
                Event[i] = 1;
            end
            # Check if individual internally migrated
            if INTMIG_OUT[i] & any(IM_OUT.==i)
                # Get index of 2nd observation
                i1 = IM_IN[IM_OUT.==i][1];
                # Copy event times and indicators for progression to dormant infection/recovery to 2nd observation
                tA[i1] = tA[i];
                tE[i1] = tE[i];
                tI[i1] = tI[i];
                tDor[i1] = tDor[i];
                tP[i1] = tP[i];
                tR[i1] = tR[i];
                infDor[i1] = infDor[i];
                infRec[i1] = infRec[i];
                AsxDor[i1] = AsxDor[i];
                AsxRec[i1] = AsxRec[i];
            end
        elseif (inf[i] & infRec[i]) | (Asx[i] & AsxDor[i]) # add 3 for direct recovery from KA and for progression from asx infection to dormant infection
            Event[i] = 3;
        elseif Asx[i] & AsxRec[i] # add 5 for direct recovery from asymptomatic infection
            Event[i] = 5;
        elseif !(t<tIM[i]) & !(t>tD[i]) & !(t>tEM[i]) ##if (tIM[i]!=t) (!INTMIG_IN[i] | (INTMIG_IN[i] & all(IM_IN.!=i))) & # individual didn't internally migrate in and hasn't already died/emigrated
            Event[i] = 1; # add 1 for birth, progression to KA, progression from KA to dormant infection, progression to PKDL, and recovery from PKDL
        end
        # elseif (tB[i]==t) & ismissing(tIM[i]) & !INTMIG_IN[i] #if !(tD[i]<t) & !(tEM[i]<t) & #if (tIM[i]!=t) # individual didn't internally migrate in and hasn't already died/emigrated
        #     Event[i] = 1; # add 1 for birth if individual didn't migrate in
        # elseif (tIM[i]==t) & !INTMIG_IN[i] #if !(tD[i]<t) & !(tEM[i]<t) & #if (tIM[i]!=t) # individual didn't internally migrate in and hasn't already died/emigrated
        #     Event[i] = 1; # add 1 for external immigration
        # elseif (tI[i]==t) & !INTMIG_IN[i] #if !(tD[i]<t) & !(tEM[i]<t) & #if (tIM[i]!=t) # individual didn't internally migrate in and hasn't already died/emigrated
        #     Event[i] = 1; # add 1 progression to KA
        # elseif (tDor[i]==t) & infDor[i] & !INTMIG_IN[i] #if !(tD[i]<t) & !(tEM[i]<t) & #if (tIM[i]!=t) # individual didn't internally migrate in and hasn't already died/emigrated
        #     Event[i] = 1; # add 1 for progression from KA to dormant infection
        # elseif (tP[i]==t) & !INTMIG_IN[i] #if !(tD[i]<t) & !(tEM[i]<t) & #if (tIM[i]!=t) # individual didn't internally migrate in and hasn't already died/emigrated
        #     Event[i] = 1; # add 1 for progression to PKDL
        # elseif (tR[i]==t) & (infDor[i] | AsxDor[i]) & !INTMIG_IN[i]
        #     Event[i] = 1; # add 1 for recovery from PKDL
        # end
    end

    # Update status vector
    Status = Status+Event;

    return t, Status, tA, tE, tI, tDor, tP, tR, hv, infDor, infRec, AsxDor, AsxRec
end

function Simulate(n,w,tmax,Status,P02_10,p2,r1,p1,r0,p0,hs,ps,dHH,dHHsqrd,alpha,typ,nHH,f,beta,delta,d0,r3,p3,r5,p5,pI,pP,pD,epsilon,ib,h4,h0,tB,tD,tIM,tEM,INTMIG_OUT,INTMIG_IN,IM_OUT,IM_IN)
# function Simulate(n,w,tmax,Status,P02_10,p2,r1,p1,r0,p0,hs,ps,dHH,dHHsqrd,alpha,typ,nHH,f,beta,delta,d0,r3,p3,r5,p5,pI,pP,epsilon,ib,h0,tB,tD)
    # n,prevStatus,rateHH,epsilon,ib,sigma,gamma,h0,hs,ps,hv,t[i-1],tB,tA,tE,tI,tDor,tP,tR,tD,thetaA,aE,thetaE,aI,thetaI,aD,thetaD,aP,thetaP,pI,pP,infDor,infRec,AsxDor,AsxRec
    # i = 1;
    t = [w];
    PreB = (Status.==0);
    Sus = (Status.==1);
    Asx = (Status.==2);
    Exp = (Status.==3);
    inf = (Status.==4);
    Dor = (Status.==5);
    PKDL = (Status.==6);
    Rec = (Status.==7);
    S = [sum(Sus)];
    A = [sum(Asx)];
    E = [sum(Exp)];
    I = [sum(inf)];
    D = [sum(Dor)];
    P = [sum(PKDL)];
    R = [sum(Rec)];
    tA = Inf*ones(n);
    tE = Inf*ones(n);
    tI = Inf*ones(n);
    tDor = Inf*ones(n);
    tP = Inf*ones(n);
    tR = Inf*ones(n);

    # Check whether it's right to use info about subsequent PKDL or I should sample which KA cases develop PKDL
    infDor = (inf .& P02_10);
    infRec = (inf .& .!P02_10); #falses(n);
    # tprev = t[1];
    AsxDor = falses(n); # all PKDL cases w/o prior KA assumed to have been asx infected during study
    AsxRec = Asx;

    # Draw recovery times for initially asx individuals, and onset and treatment times for initially presx and sx individuals
    # Need to split sx infections into non-PKDL and PKDL
    # Use size-biased distribution for non-geometric waiting times
    tR[Asx] = t .+ rand(Geometric(p2),A[1]) .+ 1;
    tI[Exp] = t .+ randSizeBiasedNegativeBinomial(r1,p1,E[1]);
    tR[Exp] = tI[Exp] + rand(NegativeBinomial(r0[1],p0[1]),E[1]) .+ 1;
    tDor[infDor] = t .+ randSizeBiasedNegativeBinomial(r0[1],p0[1],sum(infDor));
    tP[infDor] = tDor[infDor] + rand(NegativeBinomial(r3,p3),sum(infDor));
    tR[infDor] = tP[infDor] + rand(NegativeBinomial(r5,p5),sum(infDor)) .+ 1;
    tR[infRec] = t .+ randSizeBiasedNegativeBinomial(r0[1],p0[1],sum(infRec));
    tP[Dor] = t .+ randSizeBiasedNegativeBinomial(r3,p3,D[1]);
    tR[Dor] = tP[Dor] + rand(NegativeBinomial(r5,p5),D[1]) .+ 1;
    tR[PKDL] = t .+ randSizeBiasedNegativeBinomial(r5,p5,P[1]);

    hv = zeros(n,1);
    hv[infDor] = hs[rand(Categorical(ps),sum(infDor))];
    hv[Dor] = hs[rand(Categorical(ps),D[1])];
    hv[PKDL] = hs[rand(Categorical(ps),P[1])];

    KHH, K0 = KnlHH(dHH,dHHsqrd,alpha,typ,n,nHH,f);
    rateHH = beta*KHH + delta*d0;

    prevStatus = Status;

    # @time begin
        # while (t[i-1]<tmax) & ((E[i-1]>0) | (I[i-1]>0))
        for i=1:tmax-w
            # global prevStatus, Status, S, A, E, I, D, P, R, hv, t, tA, tE, tI, tDor, tP, tR, i, infDor, infRec, AsxDor, AsxRec
            # tcurr, currStatus, tA, tE, tI, tDor, tP, tR, hv, infDor, infRec, AsxDor, AsxRec = Iterate(n,prevStatus,rateHH,epsilon,ib,h0,hs,ps,hv,t[i],tB,tA,tE,tI,tDor,tP,tR,tD,p2,r1,p1,r0,p0,r3,p3,r5,p5,pI,pP,infDor,infRec,AsxDor,AsxRec);
            tcurr, currStatus, tA, tE, tI, tDor, tP, tR, hv, infDor, infRec, AsxDor, AsxRec = Iterate(n,prevStatus,rateHH,epsilon,ib,h4,h0,hs,ps,hv,t[i],tB,tA,tE,tI,tDor,tP,tR,tD,tIM,tEM,p2,r1,p1,r0[min(max(1,floor(Int64,(w+i-1)/12)),9)],p0[min(max(1,floor(Int64,(w+i-1)/12)),9)],r3,p3,r5,p5,pI,pP,pD,infDor,infRec,AsxDor,AsxRec,INTMIG_OUT,INTMIG_IN,IM_OUT,IM_IN);
            Status = append!(Status,currStatus)
            Sus = (currStatus.==1)
            Asx = (currStatus.==2)
            Exp = (currStatus.==3)
            inf = (currStatus.==4)
            Dor = (currStatus.==5)
            PKDL = (currStatus.==6)
            Rec = (currStatus.==7)
            S = push!(S,sum(Sus))
            A = push!(A,sum(Asx))
            E = push!(E,sum(Exp))
            I = push!(I,sum(inf))
            D = push!(D,sum(Dor))
            P = push!(P,sum(PKDL))
            R = push!(R,sum(Rec))
            t = push!(t,tcurr)
            # i = i+1
            prevStatus = currStatus
            # tprev = tcurr
        end
    # end

    # return(t,Status,S,A,E,I,D,P,R,tB,tA,tE,tI,tDor,tP,tR,tD,hv)
    return(t,Status,S,A,E,I,D,P,R,tA,tE,tI,tDor,tP,tR,hv)
end

function RunSims(df,para,startyr,startmo,endyr,endmo,w,p,p1,tAs,tRAs,tEs,tRsANONR,tRsAONR,nsmpls,nsims,r5,p5,pP,pD,hs,ps)
    n = size(df,1);

    origin = stata_month(startyr,startmo)-1;

    INTMIG_OUT = (df.INTMIG_OUT.==1);
    INTMIG_IN = (df.INTMIG_IN.==1);

    KothrObs = ((df.KA.>=df.MIG_OUT) .& INTMIG_OUT) .| ((df.KA.<df.MIG_IN) .& INTMIG_IN);
    KothrObs = convertToBool(KothrObs);
    PothrObs = ((df.PKDL.>=df.MIG_OUT) .& INTMIG_OUT) .| ((df.PKDL.<df.MIG_IN) .& INTMIG_IN);
    PothrObs = convertToBool(PothrObs);

    # Indicator for active KA at start of study
    actvK = (((df.KA.<=origin+w) .& (df.KARX.>origin+w)) .| (((df.KAYR.==startyr-1+floor(w/12)) .| ((df.KA.>origin) .& (df.KA.<=origin+w))) .& ismissing.(df.KARX))) .& ismissing.(df.MIG_IN); # OK for mod(w,12)=0, but not sure about mod(w,12)!=0
    actvK = convertToBool(actvK);

    # Indicator for previous KA at start of study
    prevK = ((df.KAYR.<startyr) .| (df.KA.<=origin+w));
    prevK[ismissing.(prevK)] .= false;
    prevK = convert(Array{Bool},prevK);
    prevK = (prevK .& .!actvK);

    # Indicator for VL during study
    K02_10 = (df.KA_1210.==1) .& (df.KAYR.>=startyr) .& .!KothrObs;
    K02_10 = convertToBool(K02_10);

    # Indicator for PKDL during study
    P02_10 = (df.PKDL_1210.==1) .& (df.PKDLYR.>=startyr) .& .!PothrObs;
    P02_10[ismissing.(P02_10)] .= false;
    P02_10 = convert(Array{Bool},P02_10);

    ANONR = (actvK .& ismissing.(df.KA) .& ismissing.(df.KARX));
    AONR = (actvK .& .!ismissing.(df.KA) .& ismissing.(df.KARX));

    tB = zeros(n); # dummy birth month of 0 for individuals missing DOB
    tB[.!ismissing.(df.DOB)] = df.DOB[.!ismissing.(df.DOB)] .- origin;
    tD = Inf*ones(n); # set death month for individuals that didn't die to infinity
    tD[.!ismissing.(df.DEATH)] = df.DEATH[.!ismissing.(df.DEATH)] .- origin;
    tIM = zeros(n); # dummy immigration month of 0 for individuals missing immigration month
    tIM[.!ismissing.(df.MIG_IN)] = df.MIG_IN[.!ismissing.(df.MIG_IN)] .- origin;
    tEM = Inf*ones(n); # set death month for individuals that didn't emigrate to infinity
    tEM[.!ismissing.(df.MIG_OUT)] = df.MIG_OUT[.!ismissing.(df.MIG_OUT)] .- origin;

    B = (tB.>w);
    # age = max.(-tB,0);
    IM = (tIM.>w);
    # EM = (.!isinf.(tEM));
    IM_OUT = findall(in.(df.RESP_ID,[df.ORIG_ID[INTMIG_IN]]) .& INTMIG_OUT);
    IM_IN = findall(in.(df.ORIG_ID,[df.RESP_ID[INTMIG_OUT]]) .& INTMIG_IN);

    # Set parameter values
    p2 = 1/5;
    r1 = 3;
    # Julia doesn't have an inbuilt function for fitting Neg Bin so use parameters from fitting in MATLAB
    OTparams = readdlm(joinpath(@__DIR__,"OTparams.csv"),',',Float64);
    r0 = OTparams[:,1];
    p0 = OTparams[:,2];
    # r0 = 1.369639756715030*ones(9,1);
    # p0 = 0.380308479571719*ones(9,1);
    r3 = 1.729749102087970;
    p3 = 0.064318501868660;
    # # MATLAB code: pars4=nbinfit(tRP(~isnan(tRP)&~PothrObs)-tP(~isnan(tRP)&~PothrObs)-1)
    # r5 = 1.183947478584665;
    # p5 = 0.065759912826616;

    # pI0 = 0.15;
    # pP = 0; #0.17;
    # pD = 16/((1-pI0)/pI0*1018); # estimate proportion of asx individuals who develop PKDL as ratio of asx incidence to incidence of PKDL w/o prior KA

    tmax = stata_month(endyr,endmo)-origin;
    # h0 = 0.02;
    # h1 = 9/26/(10/15);
    # h3 = 18/21/(10/15);
    # h2 = (h1+h3)/2;
    # hs = [h1;h2;h3];
    # h4 = 0.02;
    # ps = [101/138;31/138;6/138];
    # hu = sum(ps.*hs);
    # hs = [h1;h2;h3;hu];
    # ps = [101/190;31/190;6/190;52/190];

    # Calculate HH pairwise distances and spatial kernel
    typ = "Exp";
    dHH, ia, ib = CalcHHDists(df);
    nHH = length(ia);
    f = counts(ib,1:nHH);
    if typ == "Cauchy"
        dHHsqrd = dHH.^2;
    else
        dHHsqrd = [];
    end
    d0 = sparse(Matrix{Float64}(LinearAlgebra.I,nHH,nHH));

    # Process MCMC output for missing data to determine initial conditions
    nsmpls1 = size(tAs,2);
    tEs1 = Matrix{Union{Missing,Int64}}(missing,n,nsmpls1);
    tEs1[K02_10,:] = tEs;
    tRsANONR1 = Matrix{Union{Missing, Int64}}(missing,n,nsmpls1);
    if sum(ANONR)!=0
        tRsANONR1[ANONR,:] = tRsANONR;
    end
    tRsAONR1 = Matrix{Union{Missing, Int64}}(missing,n,nsmpls1);
    if sum(AONR)!=0
        tRsAONR1[AONR,:] = tRsAONR;
    end

    Asx0s = (tAs.<=w) .& (tRAs.>w) .& (tD .> t) .& (tEM .> t);
    Exp0s = (tEs1.<=w) .& (df.KA.>origin+w) .& .!IM .& (tD .> t) .& (tEM .> t);
    Exp0s = convertToBool(Exp0s);
    inf0s = (((df.KA.<=origin+w) .& (df.KARX.>origin+w) .& ismissing.(df.MIG_IN)) .| (tRsANONR1.>w) .| (tRsAONR1.>w)) .& .!INTMIG_IN .& (tD .> t) .& (tEM .> t);
    inf0s = convertToBool(inf0s);
    Dor0s = (prevK  .| (tRsANONR1.<=w) .| (tRsAONR1.<=w)) .& P02_10 .& .!INTMIG_IN .& (tD .> t) .& (tEM .> t);
    Dor0s = convertToBool(Dor0s);
    Rec0s = (((prevK .| (tRsANONR1.<=w) .| (tRsAONR1.<=w)) .& .!P02_10) .| ((tAs.<=w) .& (tRAs.<=w))) .& .!INTMIG_IN .& (tD .> t) .& (tEM .> t);
    Rec0s = convertToBool(Rec0s);

    Status0s = ones(n,nsmpls1);
    # Status0s[B,:] .= 0;
    Status0s[B .| IM,:] .= 0;
    Status0s[Asx0s] .= 2;
    Status0s[Exp0s] .= 3;
    Status0s[inf0s] .= 4;
    if pP!=0 # treated VL cases can develop PKDL
        Status0s[Dor0s] .= 5;
    else # no VL cases progress to PKDL
        Status0s[Dor0s] .=7;
    end
    Status0s[Rec0s] .= 7;

    # plot(vec(Status0s),bins=1:8,seriestype=:histogram)

    incA = Array{UInt16}(undef,tmax-w,nsmpls*nsims,length(para));
    incI = Array{UInt16}(undef,tmax-w,nsmpls*nsims,length(para));
    incP = Array{UInt16}(undef,tmax-w,nsmpls*nsims,length(para));
    # dfo = Array{Float64,3}(undef,n,8,nsims);
    for j=1:nsmpls
        println(j)
        for l=1:nsims
            global t, Status, S, A, E, I, D, P, R, tA, tE, tI, tDor, tP, tR, hv = Simulate(n,w,tmax,Status0s[:,j],P02_10,p2,r1,p1[j],r0,p0,hs,ps,dHH,dHHsqrd,p[j,2],typ,nHH,f,p[j,1],p[j,4],d0,r3,p3,r5,p5,p[j,7],pP,pD,p[j,3],ib,p[j,6],p[j,6],tB,tD,tIM,tEM,INTMIG_OUT,INTMIG_IN,IM_OUT,IM_IN);
            # global t, Status, S, A, E, I, D, P, R, tB, tA, tE, tI, tDor, tP, tR, tD, hv = Simulate(n,w,tmax,Status0s[:,j],P02_10,p2,tEs1[:,j],r1,p1,r0,p0,ANONR,AONR,tRs,tRsANONR1[:,j],tRsAONR1[:,j],tPs,tRPs,hs,ps,dHH,dHHsqrd,alpha,typ,nHH,f,beta,delta,d0,r3,p3,r5,p5,pI,pP,epsilon,ib,h0,tB,tD);
            # local Status, S, A, E, I, D, P, R, tA, tE, tI, tDor, tP, tR
            # global t, Status, S, A, E, I, D, P, R, tB, tA, tE, tI, tDor, tP, tR, tD, hv = Simulate(n,w,tmax,Status0s[:,j],P02_10,p2,r1,p1,r0,p0,hs,ps,dHH,dHHsqrd,alpha,typ,nHH,f,beta,delta,d0,r3,p3,r5,p5,pI,pP,epsilon,ib,h0,tB,tD);
            Status_wide = reshape(Status,n,length(t));
            # println(size(Status_wide))
            for i=1:length(t)-1
                for k=1:length(para)
                    incA[i,(j-1)*nsims+l,k]=sum((Status_wide[:,i].==1) .& (Status_wide[:,i+1].==2) .& (df.PARA.==para[k]));
                    incI[i,(j-1)*nsims+l,k]=sum((Status_wide[:,i].==3) .& (Status_wide[:,i+1].==4) .& (df.PARA.==para[k]));
                    incP[i,(j-1)*nsims+l,k]=sum((Status_wide[:,i].==5) .& (Status_wide[:,i+1].==6) .& (df.PARA.==para[k]));
                end
            end
        end
    end

    return(incA,incI,incP,origin,tmax,K02_10,P02_10,t,Status,tA,tE,tI,tDor,tP,tR,hv)

    # # Plot output
    # plot(t,S) # susceptibles
    # plot(t,[A,E,I]) # asx, presx and KA
    # plot(t,[A,E,I,R]) # asx, presx, KA and recovered
    # plot(t,[S,A,E,I,D,P,R]) # all
    # plot(t,sum([S A E I D P R],dims=2)) # total population

end

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
    P02_10sim = vec(any(Status_wide.==6,dims=2));
    DorOrRecw = ((df.KARX.>origin) .& (df.KARX.<=origin+w) .& ismissing.(df.MIG_IN)) .& .!(df.INTMIG_IN.==1);
    Dorw = DorOrRecw .& P02_10sim;
    Dorw = convertToBool(Dorw);
    tDor[Dorw] = df.KARX[Dorw].-origin;
    Recw = DorOrRecw .& .!P02_10sim;
    Recw = convertToBool(Recw)
    tR[Recw] = df.KARX[Recw].-origin;

    # Write event times into dataframe
    dfo = DataFrame([tA,tE,tI,tDor,tP,tR],[:tA,:tE,:tI,:tDor,:tP,:tR]);

    # Save simulation output
    CSV.write("sim_output" * str * ".csv",dfo)
    writedlm("status_mtrx" * str * ".csv",Status_wide,",")
    writedlm("init_status" * str * ".csv",Status_wide[:,1],",")
    writedlm("PKDL_infctsnss" * str * ".csv",hv,",")

    save("incA" * str * ".jld","incA",incA)
    save("incI" * str * ".jld","incI",incI)
    save("incP" * str * ".jld","incP",incP)

end
