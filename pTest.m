function [P U V]  = pTest(Xu,Yu,Xv,Yv,Xp,Yp,R,L)

dx = Yu(2,1) - Yu(1,1);

P = zeros(size(Xp));
U = zeros(size(Xu));
V = zeros(size(Xv));

%%%%%%%%%%%%%%The value of E might change%%%%%%
E = .5 * R;
%E = dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = sqrt( (Xp-(L+R)).^2 + Yp.^2 ) - R;

id = -E <= z & z <= E;
id2 = z < -E;

P(id) = -1/(2*R)*(1 - z(id)/E - sin(pi*z(id)/E)/pi);
P(id2) = -1/(R);

end