% 宇宙機がuplinkを受信する時の，推定値に対応する観測量を計算する
% 入力
% X_star: リファレンスの状態量
% xve: uplinkを送信した時刻の地球
% xvg: uplinkを送信した時刻の地上局


function Y_star = calcG1w_ur(X_star,xve,xvg,constant)
deltaT = X_star(1);
xvs    = X_star(2:7);
%% 測角(受信)
direction_ur = xve(1:3) + xvg(1:3) - xvs(1:3)...
                + norm(xve(1:3) + xvg(1:3) - xvs(1:3)) * xvs(4:6)/ constant.lightSpeed;
azm_ur = atan2(direction_ur(2), direction_ur(1));
elv_ur = atan(direction_ur(3)/(direction_ur(1)^2 +direction_ur(2)^2)^0.5);
%% 加速度
velAccel = CelestialBody.twobody(xvs,constant.sunMu,0);
accel    = velAccel(4:6);
%% 測角(送信)
direction_ut = xvs(1:3) - xve(1:3) - xvg(1:3);
azm_ut = atan2(direction_ut(2), direction_ut(1));
elv_ut = atan(direction_ut(3)/(direction_ut(1)^2 +direction_ut(2)^2)^0.5);
% 測距(1way)
length1w = norm(xve(1:3) + xvg(1:3) - xvs(1:3)) + deltaT * constant.lightSpeed;

Y_star.azm_ur      = azm_ur;
Y_star.elv_ur      = elv_ur;
Y_star.azm_ut      = azm_ut;
Y_star.elv_ut      = elv_ut;
Y_star.accel_ur    = accel;
Y_star.length1w_ur =length1w;

% Y_star = [azm_ur; elv_ur; azm_ut; elv_ut; accel;length1w];
end

