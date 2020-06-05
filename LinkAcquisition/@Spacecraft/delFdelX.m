% EKFのための伝搬式の微分を求める
function A  = delFdelX(xv,mu)
% 初期化
A = zeros(7);
A(2,5) = 1; % 位置の微分は速度 
A(3,6) = 1;
A(4,7) = 1;
R = (xv(1)^2+xv(2)^2+xv(3)^2)^0.5;

A(5,2) = -mu/R^3+3*mu*xv(1)^2/R^5; %uのx微分
A(5,3) = 3*mu*xv(1)*xv(2)/R^5;     %uのy微分
A(5,4) = 3*mu*xv(1)*xv(3)/R^5;     %uのz微分
A(6,2) = 3*mu*xv(1)*xv(2)/R^5;     %vのx微分
A(6,3) = -mu/R^3+3*mu*xv(2)^2/R^5; %vのy微分
A(6,4) = 3*mu*xv(2)*xv(3)/R^5;     %vのz微分
A(7,2) = 3*mu*xv(1)*xv(3)/R^5;     %wのx微分
A(7,3) = 3*mu*xv(2)*xv(3)/R^5;     %wのy微分
A(7,4) = -mu/R^3+3*mu*xv(3)^2/R^5; %wのz微分

end