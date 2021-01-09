% 宇宙機がuplinkを受信する時の，推定値に対応する観測量を計算する
% 入力
% X_star: リファレンスの状態量
% xve: uplinkを送信した時刻の地球
% xvg: uplinkを送信した時刻の地上局


function Y_star = calcG_ur(X_star,xvet,xvgt,xver,xvgr,dtAtGs, dt2w, constant,time,type)
deltaT = X_star(1);
xvsr    = X_star(2:7);
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

Y_star.azm_ur      = azm_ur;
Y_star.elv_ur      = elv_ur;
Y_star.azm_ut      = azm_ut;
Y_star.elv_ut      = elv_ut;
Y_star.accel_ur    = accel;
Y_star.length1w_ur =length1w;

%% 2wayのみ計算する

if strcmp(type,"2way")
    xvst   = Spacecraft.timeUpdate_sc(xvsr,constant.sunMu, -dt2w, time.simDt);
    xst = xvst(1:3);
    length2w = norm(xvet(1:3) + xvgt(1:3) - xvsr(1:3)) + norm(xver(1:3) + xvgr(1:3) - xst(1:3)) + dtAtGs * constant.lightSpeed;
    % 地上局の観測(ダウンリンクの受信)
    direction_dr = xvst(1:3) - xver(1:3) - xvgr(1:3)...
                + norm(xvst(1:3) - xver(1:3) - xvgr(1:3)) * (xver(4:6) + xvgr(4:6))/ constant.lightSpeed;
    azm_dr   =  atan2(direction_dr(2), direction_dr(1));
    elv_dr   =  atan(direction_dr(3)/(direction_dr(1)^2 +direction_dr(2)^2)^0.5);
    Y_star.length2w_ur = length2w;
    Y_star.azm_dr      = azm_dr;
    Y_star.elv_dr      = elv_dr;
end



% Y_star = [azm_ur; elv_ur; azm_ut; elv_ut; accel;length1w];
end

