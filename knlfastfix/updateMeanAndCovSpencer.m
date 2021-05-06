function [mean_new,cov_new]=updateMeanAndCovSpencer(mean_old,cov_old,p,i,i0)
% Add 1 to p row index as can't index from 0
fi=floor(i/2);
if fi==floor((i-1)/2)+1
    mean_new=mean_old+(p(i+1,:)-p(fi,:))/(i-fi+1);
    cov_new=cov_old+(p(i+1,:)'*p(i+1,:)-p(fi,:)'*p(fi,:)+(i-fi+1)*(mean_old'*mean_old-mean_new'*mean_new))/(i-fi+i0);
else
    mean_new=((i-fi)*mean_old+p(i+1,:))/(i-fi+1);
    cov_new=((i-1-fi+i0)*cov_old+p(i+1,:)'*p(i+1,:)+(i-fi)*(mean_old'*mean_old)-(i-fi+1)*(mean_new'*mean_new))/(i-fi+i0);
end