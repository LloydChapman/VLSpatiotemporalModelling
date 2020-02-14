function saveas2(h,name,flag,frmt)
if flag && nargin==3
    saveas(h,name)
elseif flag && ~isempty(frmt)
    saveas(h,name,frmt)
end