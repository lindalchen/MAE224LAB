% Square pyramid
a = 6;
b = 1;
p = 0.5895;
v = sqrt(3)/12*b*a^2;

v_disp = v*p;

h = (4*v_disp*b^2/(sqrt(3)*a^2))^(1/3);

i_xx = (a/b)^4*h^4*(3*sqrt(3)/8)/36;

l_cg = 0.75*(b-h);

lgm = i_xx/v_disp-l_cg;