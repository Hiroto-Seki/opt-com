% 探査機が推定するEKFのための観測式の微分を求める


function H = delGdelX2w_ur(X_star,xvet,xvgt,xver,xvgr, dt2w, constant,mu)
    deltaT = X_star(1);
    xvsr    = X_star(2:7);
    xst     = X_star(2:4) - dt2w* X_star(5:7);
    c = constant.lightSpeed;

    %% 頻繁に出てくるもの
    % 距離(uplink)
    Lu = ((xvet(1) + xvgt(1) - xvsr(1))^2 + (xvet(2) + xvgt(2) - xvsr(2))^2 +(xvet(3) + xvgt(3) - xvsr(3))^2)^0.5;
    % 距離(downlink)
    Ld = ((xver(1) + xvgr(1) - xst(1))^2 + (xver(2) + xvgr(2) - xst(2))^2 +(xver(3) + xvgr(3) - xst(3))^2)^0.5;
    
    % 太陽からの距離
    R = (xvsr(1)^2 + xvsr(2)^2 + xvsr(3)^2)^0.5;
    % uplinkを受信する方向 (Direction uplink receive)
    DurX = xvet(1) + xvgt(1) - xvsr(1) + xvsr(4) * Lu/c;
    DurY = xvet(2) + xvgt(2) - xvsr(2) + xvsr(5) * Lu/c;
    DurZ = xvet(3) + xvgt(3) - xvsr(3) + xvsr(6) * Lu/c;
    % uplinkを送信した方向 (Direction uplink transmit)
    DutX = xvsr(1) - xvet(1) - xvgt(1);
    DutY = xvsr(2) - xvet(2) - xvgt(2);
    DutZ = xvsr(3) - xvet(3) - xvgt(3);

%% 頻繁に出てくるものの状態量微分
    % Lの微分
    delLu(1) = 0;
    delLu(2) = (xvsr(1) - xvet(1) - xvgt(1))/Lu;
    delLu(3) = (xvsr(2) - xvet(2) - xvgt(2))/Lu;
    delLu(4) = (xvsr(3) - xvet(3) - xvgt(3))/Lu;
    delLu(5) = 0;
    delLu(6) = 0;
    delLu(7) = 0;
    % Rの微分
    delR(1) = 0;
    delR(2) = xvsr(1)/R;
    delR(3) = xvsr(2)/R;
    delR(4) = xvsr(3)/R;
    delR(5) = 0;
    delR(6) = 0;
    delR(7) = 0;
    % DurXの微分
    delDurX(1) =       xvsr(4) * delLu(1)/c;
    delDurX(2) = -1 +  xvsr(4) * delLu(2)/c;
    delDurX(3) =       xvsr(4) * delLu(3)/c;
    delDurX(4) =       xvsr(4) * delLu(4)/c;
    delDurX(5) = Lu/c + xvsr(4) * delLu(5)/c;
    delDurX(6) =       xvsr(4) * delLu(6)/c;
    delDurX(7) =       xvsr(4) * delLu(7)/c;
    % DurYの微分
    delDurY(1) =       xvsr(5) * delLu(1)/c;
    delDurY(2) =       xvsr(5) * delLu(2)/c;
    delDurY(3) = -1 +  xvsr(5) * delLu(3)/c;
    delDurY(4) =       xvsr(5) * delLu(4)/c;
    delDurY(5) =       xvsr(5) * delLu(5)/c;
    delDurY(6) = Lu/c + xvsr(5) * delLu(6)/c;
    delDurY(7) =       xvsr(5) * delLu(7)/c;
    % DurZの微分
    delDurZ(1) =       xvsr(6) * delLu(1)/c;
    delDurZ(2) =       xvsr(6) * delLu(2)/c;
    delDurZ(3) =       xvsr(6) * delLu(3)/c;
    delDurZ(4) = -1 +  xvsr(6) * delLu(4)/c;
    delDurZ(5) =       xvsr(6) * delLu(5)/c;
    delDurZ(6) =       xvsr(6) * delLu(6)/c;
    delDurZ(7) = Lu/c + xvsr(6) * delLu(7)/c;
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
    Rsc    = xvsr(1:3);
    delRsc = [zeros(3,1),eye(3),zeros(3,3)];

    %% 観測式の微分を求める
    %  測角(受信)
    delAzm_ur = (delDurY * DurX - delDurX * DurY)/(DurX^2 + DurY^2);
    delElv_ur = ( (DurX^2 + DurY^2)*delDurZ - (DurX * delDurX + DurY * delDurY)*DurZ )...
                /( (DurX^2 + DurY^2 + DurZ^2) * (DurX^2 + DurY^2)^0.5 ) ;
    % 加速度
    delAccel  = -mu * ( delRsc * R^-3 - 3 * Rsc * R^-4 * delR);
    % 送信方向
    delAzm_ut = (delDutY * DutX - delDutX * DutY)/(DutX^2 + DutY^2);
    delElv_ut = ( (DutX^2 + DutY^2)*delDutZ - (DutX * delDutX + DutY * delDutY)*DutZ )...
                /( (DutX^2 + DutY^2 + DutZ^2) * (DutX^2 + DutY^2)^0.5 ) ;
    % 測距(クロックのオフセットがのっている)
    delL1w    = [c,0,0,0,0,0,0] + delLu;
    
    
    %% 2way用
    % xstの微分
    delXst = [zeros(3,1),eye(3),- dt2w*eye(3)];
    % Ldの微分
    delLd(1) = 1/Ld * ((xst(1)-xver(1)-xvgr(1))*delXst(1,1) + (xst(2)-xver(2)-xvgr(2))*delXst(2,1) + (xst(3)-xver(3)-xvgr(3))*delXst(3,1));
    delLd(2) = 1/Ld * ((xst(1)-xver(1)-xvgr(1))*delXst(1,2) + (xst(2)-xver(2)-xvgr(2))*delXst(2,2) + (xst(3)-xver(3)-xvgr(3))*delXst(3,2));
    delLd(3) = 1/Ld * ((xst(1)-xver(1)-xvgr(1))*delXst(1,3) + (xst(2)-xver(2)-xvgr(2))*delXst(2,3) + (xst(3)-xver(3)-xvgr(3))*delXst(3,3));
    delLd(4) = 1/Ld * ((xst(1)-xver(1)-xvgr(1))*delXst(1,4) + (xst(2)-xver(2)-xvgr(2))*delXst(2,4) + (xst(3)-xver(3)-xvgr(3))*delXst(3,4));
    delLd(5) = 1/Ld * ((xst(1)-xver(1)-xvgr(1))*delXst(1,5) + (xst(2)-xver(2)-xvgr(2))*delXst(2,5) + (xst(3)-xver(3)-xvgr(3))*delXst(3,5));
    delLd(6) = 1/Ld * ((xst(1)-xver(1)-xvgr(1))*delXst(1,6) + (xst(2)-xver(2)-xvgr(2))*delXst(2,6) + (xst(3)-xver(3)-xvgr(3))*delXst(3,6));
    delLd(7) = 1/Ld * ((xst(1)-xver(1)-xvgr(1))*delXst(1,7) + (xst(2)-xver(2)-xvgr(2))*delXst(2,7) + (xst(3)-xver(3)-xvgr(3))*delXst(3,7));
    
    %2wayの測距の微分
    delL2w = delLd + delLu;

    %% まとめる
    H = [delAzm_ur;delElv_ur;delAccel;delAzm_ut;delElv_ut;delL1w;delL2w];
    

end