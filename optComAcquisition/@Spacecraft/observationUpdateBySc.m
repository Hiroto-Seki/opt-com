%---- 中間パラメーター -----
% X_star : アプリオリな推定値
% P_bar  : 誤差共分散行列(観測前)
% Y      : 観測値
% Y_star : アプリオリな推定値に対応する観測値
% y      : Y - Y_star
% H      : 観測方程式を状態量(=アプリオリな推定値)で微分したもの
% R      : 観測誤差共分散行列
% K      : カルマンゲイン
% x      : X - X_star
% X      : 更新した推定値
% P      : 更新した誤差共分散行列

function observationUpdateBySc(obj,scTrue,constant,type)
    % 何度目の観測か
    ur_counter = scTrue.ur_counter;
    % X_starとP_barを取得
    X_star = obj.X;
    P_bar  = obj.P;
    % Rを取得．QD
    R = obj.R;
    R(1,1) = scTrue.directionAccuracy_ur(ur_counter)^2;
    R(2,2) = R(1,1);
    
%     %% とりあえずRを書き換える　
%     R = [1e-7*eye(2),              zeros(2,6);          % 測角 (受信電力で書き換える)
%            zeros(3,2), 1e-14*eye(3),zeros(3,3);         % 加速度計
%            zeros(2,5), 1e-12*eye(2),zeros(2,1);         % uplinkの送信方向
%            zeros(1,7), 5e-1];  
    
    
    
    
    
    %% YとY_star,H取得. 観測によって場合分け
    if type == 1 %1wayの観測
        Y = [scTrue.directionObserved_ur(:,ur_counter);... % 測角
            scTrue.accelObseved_ur(:,ur_counter);...       % 加速度計
            scTrue.transDirection_ur(:,ur_counter);...     % uplink方向
            scTrue.lengthObserved_ur(ur_counter)];         % 測距
            Y_star = Spacecraft.calcG1w_ur(X_star,scTrue.eState_ur(:,ur_counter),scTrue.gsState_ur(:,ur_counter),constant,obj.mu);
            H = Spacecraft.delGdelX1w_ur(X_star,scTrue.eState_ur(:,ur_counter),scTrue.gsState_ur(:,ur_counter),constant,obj.mu);
    else         %2wayの観測

    end
    y = Y - Y_star;
    % ここで，方位角の測角の部分で2piの不連続を回避するためにの計算をする．
    for y_i = [1,6] %1番目は受信側の方位角の測角，6番目は送信側の方位角の測角
        y(y_i) = mod(y(y_i) + pi, 2*pi) - pi;
    end
        
%     %% Kを計算. うまくスケーリングすることで，計算精度をあげる
%     % 単位長さ(位置・測距), 単位時間(クロックオフセット)，単位角度(測角), 単位速度, 単位加速度を決める
%     u_length = 1e2; %[km]
%     u_time   = 1e-4; %[s]
%     u_angle  = 1e-6; %
%     u_speed  = 1e-2;
%     u_accel  = 1e-11;
% 
%     if type == 1
%         scaleY = [ (1/u_angle) * eye(2),                                  zeros(2,6); 
%                                 zeros(3,2), (1/u_accel) * eye(3),            zeros(3,3);
%                                 zeros(2,5), (1/u_angle) * eye(2),            zeros(2,1);
%                                                        zeros(1,7),(1/u_length)*eye(1,1)];
%     else
%         
%     end                                    
%     scaleX = [ (1/u_time) * eye(1),                                  zeros(1,6); 
%                            zeros(3,1), (1/u_length) * eye(3),           zeros(3,3);
%                                                   zeros(3,4), (1/u_speed) * eye(3)];                                            
%     % まず，y, P_bar, H,Rをスケーリングする→ys, P_bars, Hs,Rs
%     ys = scaleY * y;
%     P_bars = scaleX * P_bar * scaleX;
%  
%     Hs = scaleY * H /(scaleX);
%     Rs = scaleY^2 * obj.R;
%     % スケーリングしたys, P_bars, HsからスケーリングされたKsを得る
%     Ks = P_bars * Hs.'/(Hs*P_bars*Hs.' + Rs);
%     xs = Ks * ys;
%     x = scaleX \ xs;
%     % KsをKに直す 
%     K = scaleX\  Ks * scaleY;
    
    K = P_bar * H.'/(H*P_bar*H.' + R);
    x = K * y;
    % XとPを計算
    X = X_star +x;
    P = (eye(7) - K * H)*P_bar;
    obj.X = X;
    obj.P = P;

end
