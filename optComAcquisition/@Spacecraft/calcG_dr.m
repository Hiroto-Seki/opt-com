% 宇宙機がuplinkを受信する時の，推定値に対応する観測量を計算する
% 入力
%     % uplinkを送信した時の地球+地上局の位置・速度
%     xv_ut = earth.state_ut(:,dr_counter) +  gsTrue.state_ut(:,dr_counter);
%     % downlinkを受信した時刻の地球+地上局の位置・速度
%     xv_dr = earth.state_dr(:,dr_counter) +  gsTrue.state_dr(:,dr_counter);
%     % 宇宙機が受信してから送信するまでの時間
%     dtAtSc = gsTrue.durationAtSc(dr_counter);

function Y_star = calcG_dr(X_star,xv_ut,xv_dr,dtAtSc,constant,mu)
deltaT = X_star(1);
xvs_dt    = X_star(2:7);
% xs_ur(宇宙機がuplinkを受信する時の位置)
xs_ur = xvs_dt(1:3) - dtAtSc * xvs_dt(4:6);

%% 測角(受信)
direction_dr = xvs_dt(1:3) - xv_dr(1:3)...
                + norm(xvs_dt(1:3) - xv_dr(1:3)) * xv_dr(4:6)/ constant.lightSpeed;
azm_dr = atan2(direction_dr(2), direction_dr(1));
elv_dr = atan(direction_dr(3)/(direction_dr(1)^2 +direction_dr(2)^2)^0.5);
%% 加速度
velAccel = CelestialBody.twobody(xvs_dt,mu,0);
accel    = velAccel(4:6);
% 測距(1way)
length1w = norm(xvs_dt(1:3) - xv_dr(1:3)) - deltaT * constant.lightSpeed;
% 測距(2way)
length2w = norm(xvs_dt(1:3) - xv_dr(1:3)) + norm(xs_ur - xv_ut(1:3)) + dtAtSc * constant.lightSpeed;

Y_star = [azm_dr; elv_dr; accel;length1w;length2w];
end