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

function observationUpdateBySc(obj,scTrue,earth, gsTrue,constant,type)
    % 何度目の観測か
    ur_counter = scTrue.ur_counter;
    ur2w_counter = scTrue.ur2w_counter;
    % X_starとP_barを取得
    X_star = obj.X;
    P_bar  = obj.P;
    
    % 観測方程式に出てくるもの
    xve_ut = earth.state_ut(:,ur_counter);
    xvg_ut = gsTrue.state_ut(:,ur_counter);
    
        
    %% YとY_star,H取得. 観測によって場合分け
    if type == 1 %1wayの観測
        Y = [scTrue.directionObserved_ur(:,ur_counter);... % 測角
            scTrue.accelObseved_ur(:,ur_counter);...       % 加速度計
            scTrue.transDirection_ur(:,ur_counter);...     % uplink方向
            scTrue.lengthObserved_ur(ur_counter)];         % 測距
        Y_star = Spacecraft.calcG1w_ur(X_star,xve_ut,xvg_ut,constant,obj.mu);
        H = Spacecraft.delGdelX1w_ur(X_star,xve_ut,xvg_ut,constant,obj.mu);
        R = obj.R1wSc;
    else %2wayの観測
        Y = [scTrue.directionObserved_ur(:,ur_counter);... % 測角
            scTrue.accelObseved_ur(:,ur_counter);...       % 加速度計
            scTrue.transDirection_ur(:,ur_counter);...     % uplink方向
            scTrue.lengthObserved_ur(ur_counter);          % 測距1way
            scTrue.length2wObserved_ur(ur_counter);];      % 測距2way
        % 観測方程式に出てくるもの
        xve_dr = earth.state_dr(:,ur2w_counter);
        xvg_dr = gsTrue.state_dr(:,ur2w_counter);
        dtAtGs = gsTrue.t_ut(ur_counter) - gsTrue.t_dr(ur2w_counter);
        dt2w  = scTrue.t_ur(ur_counter) - scTrue.t_dt(ur2w_counter);
        Y_star = Spacecraft.calcG2w_ur(X_star,xve_ut,xvg_ut,...
                                        xve_dr,xvg_dr,...
                                        dtAtGs, dt2w,...
                                        constant,obj.mu);
        H = Spacecraft.delGdelX2w_ur(X_star,xve_ut,xvg_ut,...
                                        xve_dr,xvg_dr,...
                                        dt2w,...
                                        constant,obj.mu);
        R = obj.R2wSc;
    end
    y = Y - Y_star;
    obj.y = y;
    % ここで，方位角の測角の部分で2piの不連続を回避するためにの計算をする．
    for y_i = [1,6] %1番目は受信側の方位角の測角，6番目は送信側の方位角の測角
        y(y_i) = mod(y(y_i) + pi, 2*pi) - pi;
    end
    R(1,1) = scTrue.directionAccuracy_ur(ur_counter)^2;
    R(2,2) = R(1,1);    

% %     % 観測を一部無視する場合(2wayの測距を無視する)
    y = y([1 2 3 4 5 6 7 8]);
    H = H([1 2 3 4 5 6 7 8],:);
    R = R([1 2 3 4 5 6 7 8],[1 2 3 4 5 6 7 8]);

    
    K = P_bar * H.'/(H*P_bar*H.' + R);
    x = K * y;
    % XとPを計算
    X = X_star +x;
    P = (eye(7) - K * H)*P_bar;
    obj.X = X;
    obj.P = P;
    obj.H = H;
% %     if type == 1
% %         obj.R1wSc = R;
% %     else
% %         obj.R2wSc = R;
% %     end
    
    
%     if ur_counter == 494
%         disp("stop")
%     end

end
