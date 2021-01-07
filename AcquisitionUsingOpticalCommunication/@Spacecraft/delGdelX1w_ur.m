% 探査機が推定するEKFのための観測式の微分を求める


function H = delGdelX1w_ur(X_star,xve,xvg,constant)
    deltaT = X_star(1);
    xvs    = X_star(2:7);
    c = constant.lightSpeed;

    %% 頻繁に出てくるもの
    % 距離
    dr = xvs(1:3) - xve(1:3) - xvg(1:3);
    L = norm(dr);
    % 太陽からの距離
    R = norm(xvs(1:3));
    % uplinkを受信する方向 (Direction uplink receive)
    DurX = xve(1) + xvg(1) - xvs(1) + xvs(4) * L/c;
    DurY = xve(2) + xvg(2) - xvs(2) + xvs(5) * L/c;
    DurZ = xve(3) + xvg(3) - xvs(3) + xvs(6) * L/c;
    % uplinkを送信した方向 (Direction uplink transmit)
    DutX = xvs(1) - xve(1) - xvg(1);
    DutY = xvs(2) - xve(2) - xvg(2);
    DutZ = xvs(3) - xve(3) - xvg(3);

%% 頻繁に出てくるものの状態量微分
    % Lの微分
    delL(1) = 0;
    delL(2:4) = dr.'/L;
    delL(5:7) = 0;
    % Rの微分
    delR(1) = 0;
    delR(2) = xvs(1)/R;
    delR(3) = xvs(2)/R;
    delR(4) = xvs(3)/R;
    delR(5) = 0;
    delR(6) = 0;
    delR(7) = 0;
    % DurXの微分
    delDurX(1) =       xvs(4) * delL(1)/c;
    delDurX(2) = -1 +  xvs(4) * delL(2)/c;
    delDurX(3) =       xvs(4) * delL(3)/c;
    delDurX(4) =       xvs(4) * delL(4)/c;
    delDurX(5) = L/c + xvs(4) * delL(5)/c;
    delDurX(6) =       xvs(4) * delL(6)/c;
    delDurX(7) =       xvs(4) * delL(7)/c;
    % DurYの微分
    delDurY(1) =       xvs(5) * delL(1)/c;
    delDurY(2) =       xvs(5) * delL(2)/c;
    delDurY(3) = -1 +  xvs(5) * delL(3)/c;
    delDurY(4) =       xvs(5) * delL(4)/c;
    delDurY(5) =       xvs(5) * delL(5)/c;
    delDurY(6) = L/c + xvs(5) * delL(6)/c;
    delDurY(7) =       xvs(5) * delL(7)/c;
    % DurZの微分
    delDurZ(1) =       xvs(6) * delL(1)/c;
    delDurZ(2) =       xvs(6) * delL(2)/c;
    delDurZ(3) =       xvs(6) * delL(3)/c;
    delDurZ(4) = -1 +  xvs(6) * delL(4)/c;
    delDurZ(5) =       xvs(6) * delL(5)/c;
    delDurZ(6) =       xvs(6) * delL(6)/c;
    delDurZ(7) = L/c + xvs(6) * delL(7)/c;
    % DutXの微分
    delDutX(1) =  0;
    delDutX(2) =  1;
    delDutX(3) =  0;
    delDutX(4) =  0;
    delDutX(5) =  0;
    delDutX(6) =  0;
    delDutX(7) =  0;
    % DutYの微分
    delDutY(1) =  0;
    delDutY(2) =  0;
    delDutY(3) =  1;
    delDutY(4) =  0;
    delDutY(5) =  0;
    delDutY(6) =  0;
    delDutY(7) =  0;
    % DutZの微分
    delDutZ(1) =  0;
    delDutZ(2) =  0;
    delDutZ(3) =  0;
    delDutZ(4) =  1;
    delDutZ(5) =  0;
    delDutZ(6) =  0;
    delDutZ(7) =  0;
    % 宇宙機の位置ベクトルの微分
    Rsc    = xvs(1:3);
    delRsc = [zeros(3,1),eye(3),zeros(3,3)];

    %% 観測式の微分を求める
    %  測角(受信)
    delAzm_ur = (delDurY * DurX - delDurX * DurY)/(DurX^2 + DurY^2);
    delElv_ur = ( (DurX^2 + DurY^2)*delDurZ - (DurX * delDurX + DurY * delDurY)*DurZ )...
                /( (DurX^2 + DurY^2 + DurZ^2) * (DurX^2 + DurY^2)^0.5 ) ;
    % 加速度
    delAccel  = -constant.sunMu * ( delRsc * R^-3 - 3 * Rsc * R^-4 * delR);
    % 送信方向
    delAzm_ut = (delDutY * DutX - delDutX * DutY)/(DutX^2 + DutY^2);
    delElv_ut = ( (DutX^2 + DutY^2)*delDutZ - (DutX * delDutX + DutY * delDutY)*DutZ )...
                /( (DutX^2 + DutY^2 + DutZ^2) * (DutX^2 + DutY^2)^0.5 ) ;
    % 測距(クロックのオフセットがのっている)
    delL1w    = [c,0,0,0,0,0,0] + delL;
    
    %% まとめる
%     H = [delAzm_ur;delElv_ur;delAzm_ut;delElv_ut;delAccel;delL1w];
    H.azm_ur      = delAzm_ur;
    H.elv_ur      = delElv_ur;
    H.azm_ut      = delAzm_ut;
    H.elv_ut      = delElv_ut;
    H.accel_ur    = delAccel;
    H.length1w_ur = delL1w;

end