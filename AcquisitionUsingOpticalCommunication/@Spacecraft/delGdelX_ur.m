% 探査機が推定するEKFのための観測式の微分を求める


function H = delGdelX_ur(X_star,xvet,xvgt,xver,xvgr, dt2w, constant,type,time)
    deltaT = X_star(1);
    xvs    = X_star(2:7);
    c = constant.lightSpeed;

    %% 頻繁に出てくるもの
    % 距離
    dru = xvs(1:3) - xvet(1:3) - xvgt(1:3);
    Lu = norm(dru);
    % 太陽からの距離
    R = norm(xvs(1:3));
    % uplinkを受信する方向 (Direction uplink receive)
    DurX = xvet(1) + xvgt(1) - xvs(1) + xvs(4) * Lu/c;
    DurY = xvet(2) + xvgt(2) - xvs(2) + xvs(5) * Lu/c;
    DurZ = xvet(3) + xvgt(3) - xvs(3) + xvs(6) * Lu/c;
    % uplinkを送信した方向 (Direction uplink transmit)
    DutX = xvs(1) - xvet(1) - xvgt(1);
    DutY = xvs(2) - xvet(2) - xvgt(2);
    DutZ = xvs(3) - xvet(3) - xvgt(3);

%% 頻繁に出てくるものの状態量微分
    % Luの微分
    delLu(1) = 0;
    delLu(2:4) = dru.'/Lu;
    delLu(5:7) = 0;
    % Rの微分
    delR(1) = 0;
    delR(2) = xvs(1)/R;
    delR(3) = xvs(2)/R;
    delR(4) = xvs(3)/R;
    delR(5) = 0;
    delR(6) = 0;
    delR(7) = 0;
    % DurXの微分
    delDurX(1) =       xvs(4) * delLu(1)/c;
    delDurX(2) = -1 +  xvs(4) * delLu(2)/c;
    delDurX(3) =       xvs(4) * delLu(3)/c;
    delDurX(4) =       xvs(4) * delLu(4)/c;
    delDurX(5) = Lu/c + xvs(4) * delLu(5)/c;
    delDurX(6) =       xvs(4) * delLu(6)/c;
    delDurX(7) =       xvs(4) * delLu(7)/c;
    % DurYの微分
    delDurY(1) =       xvs(5) * delLu(1)/c;
    delDurY(2) =       xvs(5) * delLu(2)/c;
    delDurY(3) = -1 +  xvs(5) * delLu(3)/c;
    delDurY(4) =       xvs(5) * delLu(4)/c;
    delDurY(5) =       xvs(5) * delLu(5)/c;
    delDurY(6) = Lu/c + xvs(5) * delLu(6)/c;
    delDurY(7) =       xvs(5) * delLu(7)/c;
    % DurZの微分
    delDurZ(1) =       xvs(6) * delLu(1)/c;
    delDurZ(2) =       xvs(6) * delLu(2)/c;
    delDurZ(3) =       xvs(6) * delLu(3)/c;
    delDurZ(4) = -1 +  xvs(6) * delLu(4)/c;
    delDurZ(5) =       xvs(6) * delLu(5)/c;
    delDurZ(6) =       xvs(6) * delLu(6)/c;
    delDurZ(7) = Lu/c + xvs(6) * delLu(7)/c;
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
%     delAzm_ur = (delDurY * DurX - delDurX * DurY)/(DurX^2 + DurY^2);
%     delElv_ur = ( (DurX^2 + DurY^2)*delDurZ - (DurX * delDurX + DurY * delDurY)*DurZ )...
%                 /( (DurX^2 + DurY^2 + DurZ^2) * (DurX^2 + DurY^2)^0.5 ) ;
    % 加速度
    delAccel  = -constant.sunMu * ( delRsc * R^-3 - 3 * Rsc * R^-4 * delR);
    % 送信方向
    delAzm_ut = (delDutY * DutX - delDutX * DutY)/(DutX^2 + DutY^2);
    delElv_ut = ( (DutX^2 + DutY^2)*delDutZ - (DutX * delDutX + DutY * delDutY)*DutZ )...
                /( (DutX^2 + DutY^2 + DutZ^2) * (DutX^2 + DutY^2)^0.5 ) ;
    % 測距(クロックのオフセットがのっている)
    delL1w    = [c,0,0,0,0,0,0] + delLu;
    
    % 3軸にした場合
    dur_norm = (DurX^2 + DurY^2 + DurZ^2)^0.5;
    delDirectionX_ur = 1/dur_norm^3 * ...
        (delDurX * dur_norm^2 - DurX*(DurX *delDurX + DurY *delDurY + DurZ *delDurZ )  );
    delDirectionY_ur = 1/dur_norm^3 * ...
        (delDurY * dur_norm^2 - DurY*(DurX *delDurX + DurY *delDurY + DurZ *delDurZ )  );
    delDirectionZ_ur = 1/dur_norm^3 * ...
        (delDurZ * dur_norm^2 - DurZ*(DurX *delDurX + DurY *delDurY + DurZ *delDurZ )  );
 
    dut_norm = (DutX^2 + DutY^2 + DutZ^2)^0.5;
    delDirectionX_ut = 1/dut_norm^3 * ...
        (delDutX * dut_norm^2 - DutX*(DutX *delDutX + DutY *delDutY + DutZ *delDutZ )  );
    delDirectionY_ut = 1/dut_norm^3 * ...
        (delDutY * dut_norm^2 - DutY*(DutX *delDutX + DutY *delDutY + DutZ *delDutZ )  );
     delDirectionZ_ut = 1/dut_norm^3 * ...
        (delDutZ * dut_norm^2 - DutZ*(DutX *delDutX + DutY *delDutY + DutZ *delDutZ )  );   
    
    
    %% まとめる
%     H.azm_ur      = delAzm_ur;
%     H.elv_ur      = delElv_ur;
    H.direction_ur  = [delDirectionX_ur;delDirectionY_ur;delDirectionZ_ur];
%     H.azm_ut      = delAzm_ut;
%     H.elv_ut      = delElv_ut;
    H.direction_ut  = [delDirectionX_ut;delDirectionY_ut;delDirectionZ_ut];
    H.accel_ur    = delAccel;
    H.length1w_ur = delL1w;
    
    %% 2wayの場合
    if strcmp(type,"2way")
        xvst   = Spacecraft.timeUpdate_sc(X_star(2:7),constant.sunMu, -dt2w, time.simDt);
        xst = xvst(1:3);
%         xst     = X_star(2:4) - dt2w* X_star(5:7);
%         delXst = [zeros(3,1),eye(3),- dt2w*eye(3)];
        delXvst = Spacecraft.getSTM(X_star(2:7),constant.sunMu,0,-dt2w);
        delXst  = [zeros(3,1),delXvst(1:3,:)];
        drd = (xst(1:3) - xver(1:3) - xvgr(1:3));
        Ld  = norm(drd);
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

        %% 地上局の観測(ダウンリンクの受信方向)の微分を求める
        DdrX = (xst(1) - xver(1) - xvgr(1)) + (xver(4) + xvgr(4)) * Ld/c;
        DdrY = (xst(2) - xver(2) - xvgr(2)) + (xver(5) + xvgr(5)) * Ld/c;
        DdrZ = (xst(3) - xver(3) - xvgr(3)) + (xver(6) + xvgr(6)) * Ld/c;
        % DdrXの微分
        delDdrX(1) = delXst(1,1) + (xver(4) + xvgr(4)) * delLd(1)/c;
        delDdrX(2) = delXst(1,2) + (xver(4) + xvgr(4)) * delLd(2)/c;
        delDdrX(3) = delXst(1,3) + (xver(4) + xvgr(4)) * delLd(3)/c;
        delDdrX(4) = delXst(1,4) + (xver(4) + xvgr(4)) * delLd(4)/c;
        delDdrX(5) = delXst(1,5) + (xver(4) + xvgr(4)) * delLd(5)/c;
        delDdrX(6) = delXst(1,6) + (xver(4) + xvgr(4)) * delLd(6)/c;
        delDdrX(7) = delXst(1,7) + (xver(4) + xvgr(4)) * delLd(7)/c;
        % DdrYの微分
        delDdrY(1) = delXst(2,1) + (xver(5) + xvgr(5)) * delLd(1)/c;
        delDdrY(2) = delXst(2,2) + (xver(5) + xvgr(5)) * delLd(2)/c;
        delDdrY(3) = delXst(2,3) + (xver(5) + xvgr(5)) * delLd(3)/c;
        delDdrY(4) = delXst(2,4) + (xver(5) + xvgr(5)) * delLd(4)/c;
        delDdrY(5) = delXst(2,5) + (xver(5) + xvgr(5)) * delLd(5)/c;
        delDdrY(6) = delXst(2,6) + (xver(5) + xvgr(5)) * delLd(6)/c;
        delDdrY(7) = delXst(2,7) + (xver(5) + xvgr(5)) * delLd(7)/c;   
        % DdrZの微分
        delDdrZ(1) = delXst(3,1) + (xver(6) + xvgr(6)) * delLd(1)/c;
        delDdrZ(2) = delXst(3,2) + (xver(6) + xvgr(6)) * delLd(2)/c;
        delDdrZ(3) = delXst(3,3) + (xver(6) + xvgr(6)) * delLd(3)/c;
        delDdrZ(4) = delXst(3,4) + (xver(6) + xvgr(6)) * delLd(4)/c;
        delDdrZ(5) = delXst(3,5) + (xver(6) + xvgr(6)) * delLd(5)/c;
        delDdrZ(6) = delXst(3,6) + (xver(6) + xvgr(6)) * delLd(6)/c;
        delDdrZ(7) = delXst(3,7) + (xver(6) + xvgr(6)) * delLd(7)/c;

%         delAzm_dr = (delDdrY * DdrX - delDdrX * DdrY)/(DdrX^2 + DdrY^2);
%         delElv_dr = ( (DdrX^2 + DdrY^2)*delDdrZ - (DdrX * delDdrX + DdrY * delDdrY)*DdrZ )...
%                     /( (DdrX^2 + DdrY^2 + DdrZ^2) * (DdrX^2 + DdrY^2)^0.5 ) ;
        delDirectionX_dr = 1/(DdrX^2 + DdrY^2 + DdrZ^2)^1.5 * ...
                        (delDdrX * (DdrX^2 + DdrY^2 + DdrZ^2)  + DdrX * (delDdrX * DdrX + delDdrY * DdrY + delDdrZ * DdrZ) );
        delDirectionY_dr = 1/(DdrX^2 + DdrY^2 + DdrZ^2)^1.5 * ...
                        (delDdrY * (DdrX^2 + DdrY^2 + DdrZ^2)  + DdrY * (delDdrX * DdrX + delDdrY * DdrY + delDdrZ * DdrZ) );
        delDirectionZ_dr = 1/(DdrX^2 + DdrY^2 + DdrZ^2)^1.5 * ...
                        (delDdrZ * (DdrX^2 + DdrY^2 + DdrZ^2)  + DdrZ * (delDdrX * DdrX + delDdrY * DdrY + delDdrZ * DdrZ) );

        H.length2w_ur = delL2w;
%         H.azm_dr      = delAzm_dr;
%         H.elv_dr      = delElv_dr;
        H.direction_dr = [delDirectionX_dr;delDirectionY_dr;delDirectionZ_dr];
    end
    

end