function dxdt = twobody_stateAndSTM(t,x,mu)
R = x(1:3);
v = x(4:6);
STM = reshape(x(7:42), [6,6]);
[Ur,r] = unit(R);
g = mu./r.^2;
I = eye(3);
A = zeros(6,6);
A(4:6,1:3) = -g./r*(I - 3*Ur*Ur');
A(1:3,4:6) = I;
a    = -mu*R/r^3;
dSTMdt = A * STM;
dxdt = [v;a;reshape(dSTMdt,36,1)];
end