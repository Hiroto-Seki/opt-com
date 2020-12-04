% 地上局が推定するEKFのための観測式の微分を求める

function H = delGdelX_dr(X_star,xv_ut,xv_dr,dtAtSc,constant,mu)
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
    delAzm_dr = (delD_drY * d_drX - delD_drX * d_drY)/(d_drX^2 + d_drY^2);
    delElv_dr = (delD_drZ * (d_drX^2 + d_drY^2) - d_drZ * ( d_drX * delD_drX + d_drY * delD_drY  ))...
        / ((d_drX^2 + d_drY^2 + d_drZ^2) * (d_drX^2 + d_drY^2)^0.5 );
    % 加速度
    delAccel  = - mu * ( delRsc * R^-3 - 3 *  R^-4 * Rsc  * delR);
    % 1wayの測距(クロックのオフセットがのっている)
    delL1w    = [-c,0,0,0,0,0,0] + delLd;
    % 2wayの測距
    delL2w    = delLd + delLd;
    

    H = [delAzm_dr;delElv_dr;delAccel;delL1w;delL2w];
    

end