function D_img=div_img(data,Gx_oprt,Gy_oprt,Gz_oprt,step_size)
Gx=imfilter(data,Gx_oprt);
Gy=imfilter(data,Gy_oprt);
Gz=imfilter(permute(data,[3 2 1]),Gz_oprt);
Dx=imfilter(Gx,Gx_oprt);
Dy=imfilter(Gy,Gy_oprt);
Dz=imfilter(Gz,Gz_oprt);
if step_size<300
    test=-(Dx+Dy+permute(Dz,[3 2 1]));
else
    test=-(Dx+Dy);
end
test(test<0)=0;
D_img=double((test));
end