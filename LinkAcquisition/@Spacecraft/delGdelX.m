% 探査機が推定するEKFのための観測式の微分を求める
function H_childa = delGdelX(X_bar,xve,xvg,constant)
    deltaT =  X_bar(1);
    xvs = X_bar(2:7);
    
    c = constant.lightSpeed;
    mu = constant.sunMu;
    H_childa = zeros(3,7);
    % D:距離
    D = ((xve(1) + xvg(1) - xvs(1))^2 + (xve(2) + xvg(2) - xvs(2))^2 +(xve(3) + xvg(3) - xvs(3))^2)^0.5;
    T = D/c + deltaT;
    
    % 頻繁に出てくるもの
    ax = xve(1) + xvg(1) - xvs(1) + xvs(4) * (T-deltaT);
    ay = xve(2) + xvg(2) - xvs(2) + xvs(5) * (T-deltaT);
    az = xve(3) + xvg(3) - xvs(3) + xvs(6) * (T-deltaT);
    
    % Dのxs微分
    delDdelDt = 0;
    delDdelXs = (xvs(1) - xve(1) - xvg(1))/D;
    delDdelYs = (xvs(2) - xve(2) - xvg(2))/D;
    delDdelZs = (xvs(3) - xve(3) - xvg(3))/D;
    delDdelUs = 0;
    delDdelVs = 0;
    delDdelWs = 0;
    
    delDdelXVs = [delDdelDt, delDdelXs, delDdelYs, delDdelZs, delDdelUs, delDdelVs, delDdelWs ];
    delDtdelDt = [1,0,0,0,0,0,0];
    
    % Tの微分
    delT = delDdelXVs/c + delDtdelDt;
    
    % Tの微分
    %TのdeltaT微分
    H_childa(1,:) = delT;
    
    %% axの微分
    delAxdelDt =       xvs(4) * (delT(1) - 1);
    delAxdelXs = - 1 + xvs(4) *  delT(2);
    delAxdelYs =       xvs(4) *  delT(3);
    delAxdelZs =       xvs(4) *  delT(4);
    delAxdelUs =   T + xvs(4) *  delT(5);
    delAxdelVs =       xvs(4) *  delT(6);
    delAxdelWs =       xvs(4) *  delT(7);
    
    delAx = [delAxdelDt,delAxdelXs,delAxdelYs,delAxdelZs,delAxdelUs,delAxdelVs,delAxdelWs];

    %% ayの微分
    delAydelDt =       xvs(5) * (delT(1) - 1);
    delAydelXs =       xvs(5) *  delT(2);
    delAydelYs = - 1 + xvs(5) *  delT(3);
    delAydelZs =       xvs(5) *  delT(4);
    delAydelUs =       xvs(5) *  delT(5);
    delAydelVs =   T + xvs(5) *  delT(6);
    delAydelWs =       xvs(5) *  delT(7);
    
    delAy = [delAydelDt,delAydelXs,delAydelYs,delAydelZs,delAydelUs,delAydelVs,delAydelWs];

    %% azの微分
    delAzdelDt =       xvs(6) * (delT(1) - 1);
    delAzdelXs =       xvs(6) *  delT(2);
    delAzdelYs =       xvs(6) *  delT(3);
    delAzdelZs =  -1 + xvs(6) *  delT(4);
    delAzdelUs =       xvs(6) *  delT(5);
    delAzdelVs =       xvs(6) *  delT(6);
    delAzdelWs =   T + xvs(6) *  delT(7);
    
    delAz = [delAzdelDt,delAzdelXs,delAzdelYs,delAzdelZs,delAzdelUs,delAzdelVs,delAzdelWs];
    
    % thetaの微分
    H_childa(2,:) = 1/(ax^2 + ay^2) * ( delAy * ax - delAx * ay);
    
    % phiの微分
    H_childa(3,:) = 1/(ax^2 + ay^2 + az^2)/(ax^2 + ay^2)^0.5 ...
                    * ( (ax^2 + ay^2)* delAz - (az * ax) * delAx - (az * ay) *delAy  );
                
   % 加速度の微分
   delAccelXdelXs = -mu* (-3 * D^(-4)* delDdelXs * xvs(1) + D^(-3));
   delAccelXdelYs = -mu* (-3 * D^(-4)* delDdelYs * xvs(1));
   delAccelXdelZs = -mu* (-3 * D^(-4)* delDdelZs * xvs(1));
   H_childa(4,:) = [0, delAccelXdelXs, delAccelXdelYs, delAccelXdelZs,0,0,0 ];
   delAccelYdelXs = -mu* (-3 * D^(-4)* delDdelXs * xvs(2));
   delAccelYdelYs = -mu* (-3 * D^(-4)* delDdelYs * xvs(2) + D^(-3));
   delAccelYdelZs = -mu* (-3 * D^(-4)* delDdelZs * xvs(2));
   H_childa(5,:) = [0, delAccelYdelXs, delAccelYdelYs, delAccelYdelZs,0,0,0 ];
   delAccelZdelXs = -mu* (-3 * D^(-4)* delDdelXs * xvs(3));
   delAccelZdelYs = -mu* (-3 * D^(-4)* delDdelYs * xvs(3));
   delAccelZdelZs = -mu* (-3 * D^(-4)* delDdelZs * xvs(3) + D^(-3));
   H_childa(6,:) = [0, delAccelZdelXs, delAccelZdelYs, delAccelZdelZs,0,0,0 ];
      
end