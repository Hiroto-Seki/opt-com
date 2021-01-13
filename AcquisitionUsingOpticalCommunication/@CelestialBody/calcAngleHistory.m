function time = calcAngleHistory(time,earth,gsTrue,scTrue)

time.SpeAngle  = zeros(1,length(time.list)); % Sun-Probe-Earth angle
time.upAvail   = zeros(1,length(time.list)); % Sun-Probe-Earth angleが要求角度(3度)より大きいかどうか→uplinkできる
time.SepAngle  = zeros(1,length(time.list)); % Sun-Earth-Probe angle
time.downAvail = zeros(1,length(time.list)); % Sun-Earth-Probe angleが要求角度(12度)より大きいかどうか→uplinkできる
time.elvHorizon = zeros(1,length(time.list)); % 地上局から見て宇宙機が地平線に隠れていないかの計算
time.elvAvail   = zeros(1,length(time.list)); % 地上局から見て宇宙機が地平線に隠れていないかの計算
time.comAvail  = zeros(1,length(time.list));

for i = 1:length(time.list)
    % 地球の位置(太陽中心慣性座標系)
    xeS = earth.state(1:3,i);
    % 宇宙機の位置
    xsS = scTrue.state(1:3,i);
    % 地上局の位置(地球中心慣性座標系)
    xgE = gsTrue.state(1:3,i);
    % 地上局の位置(太陽中心慣性座標系)
    xgS = xeS + xgE;
    % SPE角
    spe = acos( - xsS.' * (xeS-xsS) / (norm(xsS) * norm ((xeS-xsS)) ))/pi * 180;
    % SEP角
    sep = acos( - xeS.' * (xsS-xeS) / (norm(xeS) * norm ((xsS-xeS)) ))/pi * 180 ;
    % elvHorizon
    elvHorizon = 90 - acos( xgE.' * (xsS - xgS)/(norm(xgE) * norm(xsS - xgS) ))/pi * 180; 
    
    % 記録する
    time.SpeAngle(i)   = spe;
    time.SepAngle(i)   = sep;
    time.elvHorizon(i) = elvHorizon;
    
    % 通信ができるかの判定
    if abs(spe) > 3
        time.upAvail(i) = 1;
    end
    if abs(sep) > 12
        time.downAvail(i) = 1;
    end
    if elvHorizon > 20 %だいたい1パスが6時間くらい?
        time.elvAvail(i) = 1;
    end
    
    % 通信できるかを記録
    if time.upAvail(i) * time.downAvail(i) * time.elvAvail(i) == 1
        time.comAvail(i) = 1;
    end
    
end

end