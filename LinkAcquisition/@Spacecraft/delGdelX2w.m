% 探査機が推定するEKFのための観測式の微分を求める(観測値が2wayの場合)
function H_childa2w = delGdelX2w(xvs,xve,xvg,xveDr,xvgDr,dt,constant)
    c = constant.lightSpeed;
    mu = constant.sunMu;
    
    %% 2wayのRTLTの状態量微分を求める
    % よく使うものを記号でおく
    sA  = c;
    sB  = -c * dt - norm (xvs(1:3) - xve(1:3) - xvg(1:3));
    sC2 = xvs(4:6).' * xvs(4:6);
    sD2 = ( xveDr(1:3) + xvgDr(1:3) - xvs(1:3) ).' * ( xveDr(1:3) + xvgDr(1:3) - xvs(1:3) );
    sE  = xvs(4:6).' * ( xveDr(1:3) + xvgDr(1:3) - xvs(1:3) );
    
    sF = sA^2 - sC2;
    sG = - (sA * sB - sE);
    sH = sA^2 * sD2 + sC2*sB^2 - 2*sA * sB *sE + sE^2 - sC2*sD2;
    
    sI = (xvs(1:3) - xve(1:3) - xvg(1:3)).' * (xvs(1:3) - xve(1:3) - xvg(1:3));
    
    % sA ~ sHの xvs微分
    dA  = [0,0,0,0,0,0];
    dB  = - 1/sI^0.5 * [xvs(1) - xve(1) - xvg(1), xvs(2) - xve(2) - xvg(2), xvs(3) - xve(3) - xvg(3), 0,0,0];
    dC2 = [0,0,0, 2*xvs(4),2*xvs(5),2*xvs(6)];
    dD2 =  2* [ xvs(1) - xveDr(1) - xvgDr(1), xvs(2) - xveDr(2) - xvgDr(2), xvs(3) - xveDr(3) - xvgDr(3), 0,0 0];
    dE  = [-xvs(4),-xvs(5),-xvs(6), - xvs(1) + xveDr(1) + xvgDr(1), - xvs(2) + xveDr(2) + xvgDr(2) , - xvs(3) + xveDr(3) + xvgDr(3)];
    dF  = 2 * sA * dA - dC2;
    dG = -dA * sB - sA * dB + dE;
    dH =   2 * sA * sD2 * dA    +    sA^2 * dD2   +  2 * sB * sC2 * dB   +   sB^2 * dC2 ...
          - 2 * sB *  sE * dA    -  2 * sA * sE * dB - 2 * sA * sB * dE   +   2 * sE * dE ...
          - sD2 * dC2 - sC2 * dD2;
    
    % 2wayのRTLTの時間微分    
    dRTLT = 1/sF^2 * ( (dG + 1/2/sH^0.5*dH)*sF  - (sG + sH^0.5)*dF);
    
    %% 測角の時間微分も求めていく
    % uplinkの距離
    L = ((xve(1) + xvg(1) - xvs(1))^2 + (xve(2) + xvg(2) - xvs(2))^2 +(xve(3) + xvg(3) - xvs(3))^2)^0.5;

    % 頻繁に出てくるもの
    ax = xve(1) + xvg(1) - xvs(1) + xvs(4) * (L/c);
    ay = xve(2) + xvg(2) - xvs(2) + xvs(5) * (L/c);
    az = xve(3) + xvg(3) - xvs(3) + xvs(6) * (L/c);
    
    % Lのxvsの微分
    dL = 1/L * [xvs(1) - xve(1) - xvg(1), xvs(2) - xve(2) - xvg(2), xvs(3) - xve(3) - xvg(3), 0, 0, 0];
    
    % ax,ay,azの微分
    dAx = [ -1 + xvs(4) * dL(1)/c,        xvs(4) * dL(2)/c,        xvs(4) * dL(3)/c,...
           L/c + xvs(4) * dL(4)/c,        xvs(4) * dL(5)/c,        xvs(4) * dL(6)/c ];
    dAy = [      xvs(5) * dL(1)/c,  -1  + xvs(5) * dL(2)/c,        xvs(5) * dL(3)/c,...
                 xvs(5) * dL(4)/c,  L/c + xvs(5) * dL(5)/c,        xvs(5) * dL(6)/c ];    
    dAz = [      xvs(6) * dL(1)/c,        xvs(6) * dL(2)/c,   -1 + xvs(6) * dL(3)/c,...
                 xvs(6) * dL(4)/c,        xvs(6) * dL(5)/c,  L/c + xvs(6) * dL(6)/c ];     
             
    % 方位角の微分
    dTheta = 1/(ax^2 + ay^2) * ( dAy * ax - dAx * ay);
    % 仰角の微分
    dPhi   = 1/(ax^2 + ay^2 + az^2)/(ax^2 + ay^2)^0.5 ...
                    * ( (ax^2 + ay^2)* dAz - (az * ax) * dAx - (az * ay) *dAy  );
                
    % 加速度の微分
    dAccelX = [ -mu* (-3 * L^(-4)* dL(1) * xvs(1) + L^(-3)),...
                -mu* (-3 * L^(-4)* dL(2) * xvs(1)),...
                -mu* (-3 * L^(-4)* dL(3) * xvs(1)), 0,0,0];
    dAccelY = [ -mu* (-3 * L^(-4)* dL(1) * xvs(2) ),...
                -mu* (-3 * L^(-4)* dL(2) * xvs(2) + L^(-3)),...
                -mu* (-3 * L^(-4)* dL(3) * xvs(2)), 0,0,0];            
    dAccelZ = [ -mu* (-3 * L^(-4)* dL(1) * xvs(3) ),...
                -mu* (-3 * L^(-4)* dL(2) * xvs(3) ),...
                -mu* (-3 * L^(-4)* dL(3) * xvs(3) + L^(-3)), 0,0,0];     
    
    H_childa2w = [dRTLT;dTheta;dPhi;dAccelX;dAccelY;dAccelZ];
    
    
      
end