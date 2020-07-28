% リファレンスの状態量の場合の観測値を出力
% dt は探査機がuplinkを受けて，downlinkするまでの時間
function  Y_bar = calcG2(xvs1,xve0,xvg0,xve2,xvg2,dt, constant)
c = constant.lightSpeed;
% RTLT
l1 = norm(xvs1(1:3) - xvs1(4:6)*dt - xve0(1:3) - xvg0(1:3) );
l2 = norm(xvs1(1:3)                - xve2(1:3) - xvg2(1:3) );
RTLT = (l1 + l2)/c;
azm = atan2(c * ( xvs1(2)-xve2(2)-xvg2(2) ) + l2 * (xve2(5)+xvg2(5)),...
            c * ( xvs1(1)-xve2(1)-xvg2(1) ) + l2 * (xve2(4)+xvg2(4)));
elv = atan( ( c * ( xvs1(3)-xve2(3)-xvg2(3) ) + l2 * (xve2(6)+xvg2(6)) )...
           / ( ( c * ( xvs1(1)-xve2(1)-xvg2(1) ) + l2 * (xve2(4)+xvg2(4)) )^2 ...
           +   ( c * ( xvs1(2)-xve2(2)-xvg2(2) ) + l2 * (xve2(5)+xvg2(5)) )^2 ) ^0.5  ); 
Y_bar = [RTLT;azm;elv];
    
end