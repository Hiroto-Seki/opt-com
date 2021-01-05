function Y_star = calcG2w_ur(X_star,xvet,xvgt,xver,xvgr,dtAtGs, dt2w, constant,time)
deltaT = X_star(1);
xvsr    = X_star(2:7);

% ここを書き換える
xvst   = Spacecraft.timeUpdate_sc(xvsr,constant.sunMu, -dt2w, time.simDt);
xst = xvst(1:3);
% xst    = X_star(2:4) - dt2w * X_star(5:7);
%% 測角(受信)
direction_ur = xvet(1:3) + xvgt(1:3) - xvsr(1:3)...
                + norm(xvet(1:3) + xvgt(1:3) - xvsr(1:3)) * xvsr(4:6)/ constant.lightSpeed;
azm_ur = atan2(direction_ur(2), direction_ur(1));
elv_ur = atan(direction_ur(3)/(direction_ur(1)^2 +direction_ur(2)^2)^0.5);
%% 加速度
velAccel = CelestialBody.twobody(xvsr,constant.sunMu,0);
accel    = velAccel(4:6);
%% 測角(送信)
direction_ut = xvsr(1:3) - xvet(1:3) - xvgt(1:3);
azm_ut = atan2(direction_ut(2), direction_ut(1));
elv_ut = atan(direction_ut(3)/(direction_ut(1)^2 +direction_ut(2)^2)^0.5);
% 測距(1way)
length1w = norm(xvet(1:3) + xvgt(1:3) - xvsr(1:3)) + deltaT * constant.lightSpeed;
% 測距(2way)
length2w = norm(xvet(1:3) + xvgt(1:3) - xvsr(1:3)) + norm(xver(1:3) + xvgr(1:3) - xst(1:3)) + dtAtGs * constant.lightSpeed;

Y_star = [azm_ur; elv_ur; accel; azm_ut; elv_ut;length1w;length2w];
end