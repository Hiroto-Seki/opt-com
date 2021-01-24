% 宇宙機がuplinkを送信する時の，推定値に対応する観測量を計算する
% 入力
%     % uplinkを送信した時の地球+地上局の位置・速度
%     xv_ut = earth.state_ut(:,dr_counter) +  gsTrue.state_ut(:,dr_counter);
%     % downlinkを受信した時刻の地球+地上局の位置・速度
%     xv_dr = earth.state_dr(:,dr_counter) +  gsTrue.state_dr(:,dr_counter);
%     % 宇宙機が受信してから送信するまでの時間
%     dtAtSc = gsTrue.durationAtSc(dr_counter);

function Y_star = calcG_dr(X_star,xv_ut,xv_dr,dtAtSc,constant)
deltaT = X_star(1);
xvs_dt    = X_star(2:7);
% xs_ur(宇宙機がuplinkを受信する時の位置). 速度は短時間で変わっていないと考える
xvs_ur = [xvs_dt(1:3) - dtAtSc * xvs_dt(4:6);xvs_dt(4:6)] ;


%% 測角(地上局受信)
direction_dr = xvs_dt(1:3) - xv_dr(1:3)...
                + norm(xvs_dt(1:3) - xv_dr(1:3)) * xv_dr(4:6)/ constant.lightSpeed;
% azm_dr = atan2(direction_dr(2), direction_dr(1));
% elv_dr = atan(direction_dr(3)/(direction_dr(1)^2 +direction_dr(2)^2)^0.5);
direction_dr = direction_dr/norm(direction_dr);

%% 加速度
velAccel = CelestialBody.twobody(xvs_dt,constant.sunMu,0);
accel    = velAccel(4:6);
% 測距(1way)
length1w = norm(xvs_dt(1:3) - xv_dr(1:3)) - deltaT * constant.lightSpeed;
% 測距(2way)
length2w = norm(xvs_dt(1:3) - xv_dr(1:3)) + norm(xvs_ur(1:3) - xv_ut(1:3)) + dtAtSc * constant.lightSpeed;

%% downlinkに載っている情報に対応させる
%% 測角(uplinkの宇宙機受信)
direction_ur = xv_ut(1:3) - xvs_ur(1:3)...
                + norm(xv_ut(1:3) - xvs_ur(1:3)) * xvs_ur(4:6)/ constant.lightSpeed;
% azm_ur = atan2(direction_ur(2), direction_ur(1));
% elv_ur = atan(direction_ur(3)/(direction_ur(1)^2 +direction_ur(2)^2)^0.5);
direction_ur = direction_ur/norm(direction_ur);


%% 測角(uplinkの地上局送信)
direction_ut = xvs_ur(1:3) - xv_ut(1:3);
% azm_ut = atan2(direction_ut(2), direction_ut(1));
% elv_ut = atan(direction_ut(3)/(direction_ut(1)^2 +direction_ut(2)^2)^0.5);
% Y_star = [azm_ur;elv_ur;azm_ut;elv_ut; accel;azm_dr; elv_dr; length1w;length2w];

direction_ut = direction_ut/norm(direction_ut);


% Y_star.azm_ur      = azm_ur;
% Y_star.elv_ur      = elv_ur;
% Y_star.azm_ut      = azm_ut;
% Y_star.elv_ut      = elv_ut;
Y_star.accel_ur    = accel;
% Y_star.azm_dr      = azm_dr;
% Y_star.elv_dr      = elv_dr;
Y_star.length1w_dr =length1w;
Y_star.length2w_dr =length2w;

Y_star.direction_ur = direction_ur;
Y_star.direction_ut = direction_ut;
Y_star.direction_dr = direction_dr;

end