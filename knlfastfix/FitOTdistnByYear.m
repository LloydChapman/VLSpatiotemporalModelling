load('data_final2.mat')

startyr=2002;
endyr=2010;
nyrs=endyr-startyr+1;
r0=NaN(nyrs,1);
p0=NaN(nyrs,1);
for i=1:nyrs
    [r0(i),p0(i)]=FitOTdistn(data.KA(data.KAYR==startyr+i-1),data.KARX(data.KAYR==startyr+i-1));
    if isinf(r0(i))
        r0(i)=1;
        p0(i)=mle(data.KARX(data.KAYR==startyr+i-1)-data.KA(data.KAYR==startyr+i-1)-1,'distribution','geo');
    end
end
dlmwrite('OTparams.csv',[r0,p0])
