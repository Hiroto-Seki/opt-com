% 地上局が推定するEKFのための観測式の微分を求める

function H = delGdelX_dr(X_star,xv_ut,xv_dr,dtAtSc,constant)
    deltaT = X_star(1);
    xvs    = X_star(2:7);
    c = constant.lightSpeed;

    %% 頻繁に出てくるもの
    % 距離(ダウンリンク, クロック誤差なし)
    ld = ((xv_dr(1) - xvs(1))^2 + (xv_dr(2) - xvs(2))^2 + (xv_dr(3) - xvs(3))^2)^0.5 ;
    % 距離(アップリンク, クロック誤差なし)
    lu = ((xvs(1) - xvs(4) * dtAtSc - xv_ut(1))^2 + (xvs(2) - xvs(5) * dtAtSc - xv_ut(2))^2 + (xvs(3) - xvs(6) * dtAtSc - xv_ut(3))^2)^0.5 ;
    % 太陽からの距離
    R = (xvs(1)^2 + xvs(2)^2 + xvs(3)^2)^0.5;
    % downlinkを受信する方向
    d_drX = (xvs(1) - xv_dr(1)) + xv_dr(4) * ld/c;
    d_drY = (xvs(2) - xv_dr(2)) + xv_dr(5) * ld/c;
    d_drZ = (xvs(3) - xv_dr(3)) + xv_dr(6) * ld/c;
     
    %% 頻繁に出てくるものの状態量微分
    % ldの微分
    delLd(1) = 0;
    delLd(2) = (xvs(1) - xv_dr(1))/ld;
    delLd(3) = (xvs(2) - xv_dr(2))/ld;
    delLd(4) = (xvs(3) - xv_dr(3))/ld;
    delLd(5) = 0;
    delLd(6) = 0;
    delLd(7) = 0;
    % luの微分
    delLu(1) = 0;
    delLu(2) = (xvs(1) - xvs(4) * dtAtSc - xv_ut(1))/lu;
    delLu(3) = (xvs(2) - xvs(5) * dtAtSc - xv_ut(2))/lu;
    delLu(4) = (xvs(3) - xvs(6) * dtAtSc - xv_ut(3))/lu;
    delLu(5) = (xvs(1) - xvs(4) * dtAtSc - xv_ut(1))/lu * (-dtAtSc);
    delLu(6) = (xvs(2) - xvs(5) * dtAtSc - xv_ut(2))/lu * (-dtAtSc);
    delLu(7) = (xvs(3) - xvs(6) * dtAtSc - xv_ut(3))/lu * (-dtAtSc);   
    % Rの微分
    delR(1) = 0;
    delR(2) = xvs(1)/R;
    delR(3) = xvs(2)/R;
    delR(4) = xvs(3)/R;
    delR(5) = 0;
    delR(6) = 0;
    delR(7) = 0;
    % d_drXの微分
    delD_drX(1)  =       xv_dr(4) * delLd(1)/c;
    delD_drX(2)  =  1 +  xv_dr(4) * delLd(2)/c;
    delD_drX(3)  =       xv_dr(4) * delLd(3)/c;
    delD_drX(4)  =       xv_dr(4) * delLd(4)/c;
    delD_drX(5)  =       xv_dr(4) * delLd(5)/c;
    delD_drX(6)  =       xv_dr(4) * delLd(6)/c;
    delD_drX(7)  =       xv_dr(4) * delLd(7)/c;
    % d_drYの微分
    delD_drY(1)  =       xv_dr(5) * delLd(1)/c;
    delD_drY(2)  =       xv_dr(5) * delLd(2)/c;
    delD_drY(3)  =  1 +  xv_dr(5) * delLd(3)/c;
    delD_drY(4)  =       xv_dr(5) * delLd(4)/c;
    delD_drY(5)  =       xv_dr(5) * delLd(5)/c;
    delD_drY(6)  =       xv_dr(5) * delLd(6)/c;
    delD_drY(7)  =       xv_dr(5) * delLd(7)/c;   
    % d_drZの微分
    delD_drZ(1)  =       xv_dr(6) * delLd(1)/c;
    delD_drZ(2)  =       xv_dr(6) * delLd(2)/c;
    delD_drZ(3)  =       xv_dr(6) * delLd(3)/c;
    delD_drZ(4)  =  1 +  xv_dr(6) * delLd(4)/c;
    delD_drZ(5)  =       xv_dr(6) * delLd(5)/c;
    delD_drZ(6)  =       xv_dr(6) * delLd(6)/c;
    delD_drZ(7)  =       xv_dr(6) * delLd(7)/c;
    % 宇宙機の位置ベクトルの微分
    Rsc    = xvs(1:3);
    delRsc = [zeros(3,1),eye(3),zeros(3,3)];

    %% 観測式の微分を求める
    %  測角(受信)
%     delAzm_dr = (delD_drY * d_drX - delD_drX * d_drY)/(d_drX^2 + d_drY^2);
%     delElv_dr = (delD_drZ * (d_drX^2 + d_drY^2) - d_drZ * ( d_drX * delD_drX + d_drY * delD_drY  ))...
%         / ((d_drX^2 + d_drY^2 + d_drZ^2) * (d_drX^2 + d_drY^2)^0.5 );

    delDirectionX_dr = 1/(d_drX^2 + d_drY^2 + d_drZ^2)^1.5 * ...
                        (delD_drX * (d_drX^2 + d_drY^2 + d_drZ^2)  + d_drX * (delD_drX * d_drX + delD_drY * d_drY + delD_drZ * d_drZ) );
    delDirectionY_dr = 1/(d_drX^2 + d_drY^2 + d_drZ^2)^1.5 * ...
                        (delD_drY * (d_drX^2 + d_drY^2 + d_drZ^2)  + d_drY * (delD_drX * d_drX + delD_drY * d_drY + delD_drZ * d_drZ) );
    delDirectionZ_dr = 1/(d_drX^2 + d_drY^2 + d_drZ^2)^1.5 * ...
                        (delD_drZ * (d_drX^2 + d_drY^2 + d_drZ^2)  + d_drZ * (delD_drX * d_drX + delD_drY * d_drY + delD_drZ * d_drZ) );                
    
    % 加速度
    delAccel  = - constant.sunMu * ( delRsc * R^-3 - 3 *  R^-4 * Rsc  * delR);
    % 1wayの測距(クロックのオフセットがのっている)
    delL1w    = [-c,0,0,0,0,0,0] + delLd;
    % 2wayの測距
    delL2w    = delLd + delLd;
       
    %% uplinkに関するもの
    % uplinkの送信方向
    d_utX  = xvs(1) - xvs(4) * dtAtSc - xv_ut(1);
    d_utY  = xvs(2) - xvs(5) * dtAtSc - xv_ut(2);
    d_utZ  = xvs(3) - xvs(6) * dtAtSc - xv_ut(3);
    % d_utXの微分
    delD_utX(1) = 0;
    delD_utX(2) = 1;
    delD_utX(3) = 0;
    delD_utX(4) = 0;
    delD_utX(5) = -dtAtSc;
    delD_utX(6) = 0;
    delD_utX(7) = 0;
    % d_utYの微分
    delD_utY(1) = 0;
    delD_utY(2) = 0;
    delD_utY(3) = 1;
    delD_utY(4) = 0;
    delD_utY(5) = 0;
    delD_utY(6) = -dtAtSc;
    delD_utY(7) = 0;
    % d_utZの微分
    delD_utZ(1) = 0;
    delD_utZ(2) = 0;
    delD_utZ(3) = 0;
    delD_utZ(4) = 1;
    delD_utZ(5) = 0;
    delD_utZ(6) = 0;
    delD_utZ(7) = -dtAtSc; 
%     % 方位角
%     delAzm_ut = (delD_utY * d_utX - delD_utX * d_utY)/(d_utX^2 + d_utY^2);
%     % 仰角
%     delElv_ut = (delD_utZ * (d_utX^2 + d_utY^2) - d_utZ * ( d_utX * delD_utX + d_utY * delD_utY  ))...
%         / ((d_utX^2 + d_utY^2 + d_utZ^2) * (d_utX^2 + d_utY^2)^0.5 );
    
    delDirectionX_ut = 1/(d_utX^2 + d_utY^2 + d_utZ^2)^1.5 * ...
                        (delD_utX * (d_utX^2 + d_utY^2 + d_utZ^2)  + d_utX * (delD_utX * d_utX + delD_utY * d_utY + delD_utZ * d_utZ) );
    delDirectionY_ut = 1/(d_utX^2 + d_utY^2 + d_utZ^2)^1.5 * ...
                        (delD_utY * (d_utX^2 + d_utY^2 + d_utZ^2)  + d_utY * (delD_utX * d_utX + delD_utY * d_utY + delD_utZ * d_utZ) );
    delDirectionZ_ut = 1/(d_utX^2 + d_utY^2 + d_utZ^2)^1.5 * ...
                        (delD_utZ * (d_utX^2 + d_utY^2 + d_utZ^2)  + d_utZ * (delD_utX * d_utX + delD_utY * d_utY + delD_utZ * d_utZ) );               
    % uplinkの受信方向
    d_urX  =  xv_ut(1) - (xvs(1) - xvs(4) * dtAtSc) + xvs(4)/c * lu;
    d_urY  =  xv_ut(2) - (xvs(2) - xvs(5) * dtAtSc) + xvs(5)/c * lu;
    d_urZ  =  xv_ut(3) - (xvs(3) - xvs(6) * dtAtSc) + xvs(6)/c * lu;
    % d_urXの微分
    delD_urX(1) =                     xvs(4) * delLu(1)/c;
    delD_urX(2) = -1 +                xvs(4) * delLu(2)/c;
    delD_urX(3) =                     xvs(4) * delLu(3)/c;
    delD_urX(4) =                     xvs(4) * delLu(4)/c;
    delD_urX(5) = dtAtSc + 1/c * lu + xvs(4) * delLu(5)/c;
    delD_urX(6) =                     xvs(4) * delLu(6)/c;
    delD_urX(7) =                     xvs(4) * delLu(7)/c;
    % d_urYの微分
    delD_urY(1) =                     xvs(5) * delLu(1)/c;
    delD_urY(2) =                     xvs(5) * delLu(2)/c;
    delD_urY(3) = -1 +                xvs(5) * delLu(3)/c;
    delD_urY(4) =                     xvs(5) * delLu(4)/c;
    delD_urY(5) =                     xvs(5) * delLu(5)/c;
    delD_urY(6) = dtAtSc + 1/c * lu + xvs(5) * delLu(6)/c;
    delD_urY(7) =                     xvs(5) * delLu(7)/c;
    % d_urZの微分
    delD_urZ(1) =                     xvs(6) * delLu(1)/c;
    delD_urZ(2) =                     xvs(6) * delLu(2)/c;
    delD_urZ(3) =                     xvs(6) * delLu(3)/c;
    delD_urZ(4) = -1+                 xvs(6) * delLu(4)/c;
    delD_urZ(5) =                     xvs(6) * delLu(5)/c;
    delD_urZ(6) =                     xvs(6) * delLu(6)/c;
    delD_urZ(7) = dtAtSc + 1/c * lu + xvs(6) * delLu(7)/c;
    % 方位角
%     delAzm_ur = (delD_urY * d_urX - delD_urX * d_urY)/(d_urX^2 + d_urY^2);
%     % 仰角
%     delElv_ur = (delD_urZ * (d_urX^2 + d_urY^2) - d_urZ * ( d_urX * delD_urX + d_urY * delD_urY  ))...
%         / ((d_urX^2 + d_urY^2 + d_urZ^2) * (d_urX^2 + d_urY^2)^0.5 );
    delDirectionX_ur = 1/(d_urX^2 + d_urY^2 + d_urZ^2)^1.5 * ...
                        (delD_urX * (d_urX^2 + d_urY^2 + d_urZ^2)  + d_urX * (delD_urX * d_urX + delD_urY * d_urY + delD_urZ * d_urZ) );
    delDirectionY_ur = 1/(d_urX^2 + d_urY^2 + d_urZ^2)^1.5 * ...
                        (delD_urY * (d_urX^2 + d_urY^2 + d_urZ^2)  + d_urY * (delD_urX * d_urX + delD_urY * d_urY + delD_urZ * d_urZ) );
    delDirectionZ_ur = 1/(d_urX^2 + d_urY^2 + d_urZ^2)^1.5 * ...
                        (delD_urZ * (d_urX^2 + d_urY^2 + d_urZ^2)  + d_urZ * (delD_urX * d_urX + delD_urY * d_urY + delD_urZ * d_urZ) );   
%     H.azm_ur      = delAzm_ur;
%     H.elv_ur      = delElv_ur;
%     H.azm_ut      = delAzm_ut;
%     H.elv_ut      = delElv_ut;
    H.accel_ur    = delAccel;
%     H.azm_dr      = delAzm_dr;
%     H.elv_dr      = delElv_dr;
    H.length1w_dr = delL1w;
    H.length2w_dr = delL2w;
    
    H.direction_ur = [delDirectionX_ur;delDirectionY_ur;delDirectionZ_ur];
    H.direction_ut = [delDirectionX_ut;delDirectionY_ut;delDirectionZ_ut];
    H.direction_dr = [delDirectionX_dr;delDirectionY_dr;delDirectionZ_dr];
end