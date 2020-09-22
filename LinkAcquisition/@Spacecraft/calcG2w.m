% リファレンスの状態量の場合の観測値(2way)
% 伝搬時間の運動は等速直線運動を仮定
% 出力: T(2way), theta(方位角), phi(仰角), velAccel(4:6)(加速度)

function  Y_bar2w = calcG2w(X_bar2w,xve,xvg,xveDr,xvgDr,dt,constant)
c = constant.lightSpeed;
xvs = X_bar2w;
% まず，RTLTを求める．二次方程式を解く際によく出てくるものを記号でおく
sA = c;
sB = - c * dt - norm(xvs(1:3) - xve(1:3) - xvg(1:3) );
sC = xvs(4:6);
sD = xveDr(1:3) + xvgDr(1:3) - xvs(1:3);
% sCとsDの内積(inner product)
ipCd = sC.' * sD;
sC2 = sC.'*sC;
sD2 = sD.'*sD;

T = (- (sA * sB - ipCd) + sqrt(sA^2 * sD2  +  sC2 * sB^2  - 2 * sA * sB * ipCd + ipCd^2 - sC2 * sD2)) / (sA^2 - sC2);

ltdUp = norm( xve(1:3) + xvg(1:3) - xvs(1:3) )/c;

% 角度の計算
theta = atan2(xve(2) + xvg(2) - xvs(2) + xvs(5)*ltdUp...
                                , xve(1) + xvg(1) - xvs(1) + xvs(4)*ltdUp);
phi = atan( (xve(3) + xvg(3) - xvs(3) + xvs(6)*ltdUp)/...
                                   ((xve(2) + xvg(2) - xvs(2) + xvs(5) *ltdUp)^2 ...
                                  + (xve(1) + xvg(1) - xvs(1) + xvs(4) *ltdUp)^2)^0.5 );
%
velAccel = Spacecraft.twobody(xvs, constant.sunMu, 0);

Y_bar2w = [T;theta;phi; velAccel(4:6)];

end