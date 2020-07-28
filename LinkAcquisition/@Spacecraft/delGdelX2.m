% 地上局が推定するEKFのための観測式の微分を求める
function H_childa = delGdelX2(xvs1,xve0,xvg0,xve2,xvg2,dt, constant)
    c = constant.lightSpeed;
    H_childa = zeros(3,6);

    % 頻繁に出てくるもの
    X1 = xvs1(1) - xvs1(4)* dt - xve0(1) - xvg0(1);
    Y1 = xvs1(2) - xvs1(5)* dt - xve0(2) - xvg0(2);
    Z1 = xvs1(3) - xvs1(6)* dt - xve0(3) - xvg0(3); 
    X2 = xvs1(1) - xve2(1) - xvg2(1);
    Y2 = xvs1(2) - xve2(2) - xvg2(2);
    Z2 = xvs1(3) - xve2(3) - xvg2(3);
    % D1:アップリンク距離
    D1 = (X1^2 + Y1^2 + Z1^2)^0.5;
    % D2:ダウンリンク距離
    D2 = (X2^2 + Y2^2 + Z2^2)^0.5;
    % よく出てくるもの
    ax = c * X2 +  (xve2(4) + xvg2(4)) * D2;
    ay = c * Y2 +  (xve2(5) + xvg2(5)) * D2;
    az = c * Z2 +  (xve2(6) + xvg2(6)) * D2;
    
    % X1~Z2の微分(xvs1=[x,y,z,u,v,w])
    delX1_delXvs1 = [   1,   0,   0, -dt,   0,   0];
    delY1_delXvs1 = [   0,   1,   0,   0, -dt,   0];
    delZ1_delXvs1 = [   0,   0,   1,   0,   0, -dt];
    delX2_delXvs1 = [   1,   0,   0,   0,   0,   0];
    delY2_delXvs1 = [   0,   1,   0,   0,   0,   0];
    delZ2_delXvs1 = [   0,   0,   1,   0,   0,   0]; 
    % D1の微分
    delD1_delXvs1 = 1/D1 * (delX1_delXvs1 + delY1_delXvs1 + delZ1_delXvs1);
    % D2の微分
    delD2_delXvs1 = 1/D2 * (delX2_delXvs1 + delY2_delXvs1 + delZ2_delXvs1);
    % ax,ay,azの微分
    delAx_delXVs1 = c * delX2_delXvs1 + (xve2(4) + xvg2(4)) * delD2_delXvs1;
    delAy_delXVs1 = c * delY2_delXvs1 + (xve2(5) + xvg2(5)) * delD2_delXvs1;
    delAz_delXVs1 = c * delZ2_delXvs1 + (xve2(6) + xvg2(6)) * delD2_delXvs1;
    
    % RTLTの偏微分
    H_childa(1,:) = 1/c * delD1_delXvs1 + 1/c * delD2_delXvs1;
    
    % 方位角の微分
    H_childa(2,:) = (delAy_delXVs1 * ax - ay * delAx_delXVs1)/(ax^2 + ay^2);
    
    % 仰角の微分
    H_childa(3,:) = ( delAz_delXVs1 * (ax^2 + ay^2) - az * (ax * delAx_delXVs1 + ay * delAy_delXVs1) )...
                    / ( (ax^2 + ay^2)^0.5 * (ax^2 + ay^2 + az^2) );
      
end