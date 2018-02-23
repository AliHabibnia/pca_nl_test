function out_1 = f_1(Ss,w,Er,pc_n)

vect = (svd(Ss + Er/w))';
out_1 = min(sum(vect(:,(pc_n+1):end),2));

end