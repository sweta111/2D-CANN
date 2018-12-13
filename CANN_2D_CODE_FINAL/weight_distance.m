
function y = weight_distance(loc_x,loc_y,mm,nn,dx,sig_w)
y=zeros(mm,nn);
a = dx;
for i_x=1:mm;
    for i_y=1:nn;
        d_x=min(abs((i_x-loc_x)*a));
        d_y=min(abs((i_y-loc_y)*a));
        d = sqrt(d_x^2 + d_y^2);
        y(i_x,i_y)=exp(-d^2/sig_w.^2);
        size(y);
    end
end
return