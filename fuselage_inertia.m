function [Ixf, Iyf, Izf, Ixzf] = fuselage_inertia(fuselageradius,noseconeweight,tailconeweight,mainfuselageweight,xyplanetocenterline,noseconelength,tailconelength,mainfuselagelength)

R = fuselageradius;
W_n = noseconeweight;
W_t = tailconeweight;
W_c = mainfuselageweight;
z_b = xyplanetocenterline; % z distance from bottom of fuselage to its centerline
l_n = noseconelength;
l_t = tailconelength;
l_c = mainfuselagelength;


W_s = W_n+W_c+W_t;

Ixf = R^2/2*(W_n+2*W_c+W_t)+W_s*(z_b)^2;
Iyf = R^2/4*(W_n+2*W_c+W_t)+l_n^2*(W_n/2+W_c+W_t)+l_c^2*(W_c/3+W_t)+l_t^2*W_t/6+l_c*l_n*(W_c+2*W_t)+2/3*l_t*l_c*W_t+2/3*l_t*l_n*W_t+W_s*(z_b)^2;
Izf = Iyf-W_s*(z_b)^2;
Ixzf = W_n*(3/4*l_n*z_b)+W_c*z_b*(l_n+l_c/2)+W_t*z_b(l_n+l_c+l_t/4);

end