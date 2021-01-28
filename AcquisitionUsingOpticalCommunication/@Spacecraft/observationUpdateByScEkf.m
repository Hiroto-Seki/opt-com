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

function observationUpdateByScEkf(obj,scTrue,earth, gsTrue,constant,type,time,ekf)
    % 何度目の観測か
    ur_counter = scTrue.ur_counter;
    ur2w_counter = scTrue.ur2w_counter;
    % X_starとP_barを取得
    X_star = obj.X;
    P_bar  = obj.P;
    
    % 観測方程式に出てくるもの
    xve_ut = earth.state_ut(:,ur_counter);
    xvg_ut = gsTrue.state_ut(:,ur_counter);
    
    % Rの書き換え
    % 1,2要素目の観測は，1wayでも2wayでもuplinkの測角になっている．
    obj.R.direction_ur = scTrue.directionAccuracy_ur(ur_counter)^2;
    obj.R.direction_ut = gsTrue.directionAccuracy_ut(ur_counter)^2;
    
    
    
    %% YとY_star,H取得. 観測によって場合分け
    if type == 1 %1wayの観測 (1,2,3,4,5,6,7,8)の観測を得られる
        Y.direction_ur = scTrue.directionObserved_ur(:,ur_counter);
        Y.direction_ut = scTrue.transDirection_ur(:,ur_counter);
        Y.accel_ur     = scTrue.accelObseved_ur(:,ur_counter);
        Y.length1w_ur  = scTrue.lengthObserved_ur(ur_counter);
        Y.rangeRate_ur = scTrue.rangeRateObserved_ur(ur_counter);
        Y_star = Spacecraft.calcG_ur(X_star,xve_ut,xvg_ut,[],[],[],[],constant,[],"1way");
        H = Spacecraft.delGdelX_ur(X_star,xve_ut,xvg_ut,[],[],[],constant,"1way",time);
        
        % 観測できているかの判定を使う
        if scTrue.ur_observability(ur_counter) == 1
            obsType = "1u_noObs";
        elseif scTrue.ur_observability(ur_counter) == 2
            obsType = "1u_lowSnr";
        elseif scTrue.ur_observability(ur_counter) == 3
            obsType = "1u_Obs";
        else
            disp("observation type is not defined correctly ")
        end
        % Y, Y_star, H, Rから必要な要素だけ取り出す
        [y,Hm,Rm,~,estNoUseList] = Spacecraft.alignReqInfo4Est(Y,Y_star,H,obj.R,obsType,"ekf",obj.useObs,P_bar);
    else %2wayの観測 (1,2,3,4,5,6,7,8,9,10,11)の観測を得られる．(10,11はuplink内容に含まれる) 
        Y.direction_ur = scTrue.directionObserved_ur(:,ur_counter); % 測角
        Y.direction_ut = scTrue.transDirection_ur(:,ur_counter);...     % uplink方向
        Y.accel_ur     = scTrue.accelObseved_ur(:,ur_counter);...       % 加速度計
        Y.length1w_ur  = scTrue.lengthObserved_ur(ur_counter);          % 測距1way
        Y.length2w_ur  = scTrue.length2wObserved_ur(ur_counter);        % 測距2way
        Y.direction_dr = scTrue.recDownAngle_ur(:,ur_counter);         % uplinkされる，地上局での観測
        Y.length1w_dr  = gsTrue.lengthObserved_dr(:,ur2w_counter);     %地上局側の観測量
        Y.length2w_dr  = gsTrue.length2wObserved_dr(:,ur2w_counter);     %地上局側の観測量
        Y.rangeRate_ur = scTrue.rangeRateObserved_ur(ur_counter);
        % 観測方程式に出てくるもの
        xve_dr = earth.state_dr(:,ur2w_counter);
        xvg_dr = gsTrue.state_dr(:,ur2w_counter);
        dtAtGs = gsTrue.t_ut(ur_counter) - gsTrue.t_dr(ur2w_counter);
        dt2w  = scTrue.t_ur(ur_counter) - scTrue.t_dt(ur2w_counter);
        % 地上局の2wayの観測を交換する場合に出てくる
        xve_ut3w = earth.state_ut(:,ur2w_counter); 
        xvg_ut3w = gsTrue.state_ut(:,ur2w_counter);
        dtAtSc3w = scTrue.t_dt(ur2w_counter) - scTrue.t_ur(ur2w_counter);
        
        Y_star = Spacecraft.calcG_ur(X_star,xve_ut,xvg_ut,...
                                        xve_dr,xvg_dr,...
                                        dtAtGs, dt2w,...
                                        constant,time,"2way",xve_ut3w,xvg_ut3w,dtAtSc3w);
        H = Spacecraft.delGdelX_ur(X_star,xve_ut,xvg_ut,...
                                        xve_dr,xvg_dr,...
                                        dt2w,...
                                        constant,"2way",time,xve_ut3w,xvg_ut3w,dtAtSc3w);
        % Rの書き換え
        obj.R.direction_dr = scTrue.recDownAngleAccuracy_ur(ur_counter)^2;
        %% 観測できているのかの判定
        % downlinkで観測できていたか: gsTrue.dr_observability(ur2w_counter)
        % uplinkが観測できていたか:scTrue.ur_observability(ur_counter)
        switch gsTrue.dr_observability(ur2w_counter)
            case 1 %downlinkが観測されていない
                switch scTrue.ur_observability(ur_counter)
                    case 1
                        obsType = "2u_noObsD_noObsU";
                    case 2
                        obsType = "2u_noObsD_lowU";
                    case 3
                        obsType = "2u_noObsD_obsU";
                    otherwise
                        disp("obsType is not set correctly")
                end  
            case 2 %downlinkのSNRが低い
                switch scTrue.ur_observability(ur_counter)
                    case 1
                        obsType = "2u_lowD_noObsU";
                    case 2
                        obsType = "2u_lowD_lowU";
                    case 3
                        obsType = "2u_lowD_obsU";
                    otherwise
                        disp("obsType is not set correctly")
                end  
            case 3 %downlinkが観測されている
                switch scTrue.ur_observability(ur_counter)
                    case 1
                        obsType = "2u_obsD_noObsU";
                    case 2
                        obsType = "2u_obsD_lowU";
                    case 3
                        obsType = "2u_obsD_obsU";
                    otherwise
                        disp("obsType is not set correctly")
                end
            otherwise
                disp("obsType is not set correctly")
        end
        
        
        [y,Hm,Rm,~,estNoUseList] = Spacecraft.alignReqInfo4Est(Y,Y_star,H,obj.R,obsType,"ekf",obj.useObs,P_bar);
        
    end
    
    % estNoUseListとscEstBySc.estNoUseを比較する．同じ要素が含まれていれば，誤差共分散を大きくする
    noUseObsAgain = [];
    
    if isempty(obj.estNoUse) || isempty(estNoUseList)
        % 特に何もしない
    else
        for i = 1:length(obj.estNoUse)
            for j = 1:length(estNoUseList)
                if obj.estNoUse(i) == estNoUseList(j)
                    noUseObsAgain = [noUseObsAgain,obj.estNoUse(i)];
                end
            end
        end
    end
    
    % noUseObsAgainの要素によって，Pの大きさを変化させる
    if isempty(noUseObsAgain)
        % 特に何もしない
    else
        tempDt = P_bar(1,1)^0.5;
        tempX  = P_bar(2,2)^0.5;
        tempY  = P_bar(3,3)^0.5;
        tempZ  = P_bar(4,4)^0.5;
        tempU  = max(P_bar(5,5)^0.5); %経験値 
        tempV  = max(P_bar(6,6)^0.5);
        tempW  = max(P_bar(7,7)^0.5);
        if sum(strcmp(noUseObsAgain,"direction_urX"))== 1
            tempX  = max(tempX , abs(Y.direction_ur(1) - Y_star.direction_ur(1))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_ur(1) - Y_star.direction_ur(1))*constant.AU * 10/constant.lightSpeed);
        end
        if sum(strcmp(noUseObsAgain,"direction_urY")) == 1
            tempY  = max(tempY , abs(Y.direction_ur(2) - Y_star.direction_ur(2))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_ur(2) - Y_star.direction_ur(2))*constant.AU * 10/constant.lightSpeed);      
        end
        if sum(strcmp(noUseObsAgain,"direction_urZ")) == 1
            tempZ  = max(tempZ , abs(Y.direction_ur(3) - Y_star.direction_ur(3))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_ur(3) - Y_star.direction_ur(3))*constant.AU * 10/constant.lightSpeed);                   
        end
        if sum(strcmp(noUseObsAgain,"direction_utX")) == 1
            tempX  = max(tempX , abs(Y.direction_ut(1) - Y_star.direction_ut(1))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_ut(1) - Y_star.direction_ut(1))*constant.AU * 10/constant.lightSpeed);
        end
        if sum(strcmp(noUseObsAgain,"direction_utY")) == 1
            tempY  = max(tempY , abs(Y.direction_ut(2) - Y_star.direction_ut(2))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_ut(2) - Y_star.direction_ut(2))*constant.AU * 10/constant.lightSpeed);
        end
        if sum(strcmp(noUseObsAgain,"direction_utZ")) == 1
            tempZ  = max(tempZ , abs(Y.direction_ut(3) - Y_star.direction_ut(3))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_ut(3) - Y_star.direction_ut(3))*constant.AU * 10/constant.lightSpeed);
        end
        if sum(strcmp(noUseObsAgain,"length1w_ur")) == 1
            tempX  = max(tempX , abs(Y.length1w_ur - Y_star.length1w_ur));
            tempY  = max(tempY , abs(Y.length1w_ur - Y_star.length1w_ur));
            tempZ  = max(tempZ , abs(Y.length1w_ur - Y_star.length1w_ur));
            tempDt = max(tempDt, abs(Y.length1w_ur - Y_star.length1w_ur)/constant.lightSpeed);            
        end
        if sum(strcmp(noUseObsAgain,"length2w_ur")) == 1
            tempX  = max(tempX , abs(Y.length2w_ur - Y_star.length2w_ur));
            tempY  = max(tempY , abs(Y.length2w_ur - Y_star.length2w_ur));
            tempZ  = max(tempZ , abs(Y.length2w_ur - Y_star.length2w_ur));
            tempDt = max(tempDt, abs(Y.length2w_ur - Y_star.length2w_ur)/constant.lightSpeed);  
        end
        if sum(strcmp(noUseObsAgain,"direction_drX")) == 1
            tempX  = max(tempX , abs(Y.direction_dr(1) - Y_star.direction_dr(1))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_dr(1) - Y_star.direction_dr(1))*constant.AU * 10/constant.lightSpeed);
        end
        if sum(strcmp(noUseObsAgain,"direction_drY")) == 1
            tempY  = max(tempY , abs(Y.direction_dr(2) - Y_star.direction_dr(2))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_dr(2) - Y_star.direction_dr(2))*constant.AU * 10/constant.lightSpeed);
        end
        if sum(strcmp(noUseObsAgain,"direction_drZ")) == 1
            tempZ  = max(tempZ , abs(Y.direction_dr(3) - Y_star.direction_dr(3))*constant.AU * 10);
            tempDt = max(tempDt, abs(Y.direction_dr(3) - Y_star.direction_dr(3))*constant.AU * 10/constant.lightSpeed);
        end
        % 速度の誤差は，公転半径の差由来だと考えて，
        tempDr = (tempX^2 + tempY^2 + tempZ^2)^0.5;
        tempR  = (X_star(2)^2 + X_star(3)^2 + X_star(4)^2)^0.5;
        tempSpeed  = (((tempR+tempDr)/tempR)^0.5 -1) * (X_star(5)^2 + X_star(6)^2 + X_star(7)^2)^0.5;
        tempU  = max(P_bar(5,5)^0.5,tempSpeed);
        tempV  = max(P_bar(6,6)^0.5,tempSpeed);
        tempW  = max(P_bar(7,7)^0.5,tempSpeed);
    end
    
   % 特定の観測が連続で除外されていた場合はPを書き換え，Xはそのまま
   if isempty(noUseObsAgain) == 0
       X = X_star;
       P = [  tempDt^2, zeros(1,6);
            zeros(1,1),    tempX^2, zeros(1,5);
            zeros(1,2),    tempY^2, zeros(1,4);
            zeros(1,3),    tempZ^2, zeros(1,3);
            zeros(1,4),    tempU^2, zeros(1,2);
            zeros(1,5),    tempV^2, zeros(1,1);
            zeros(1,6),    tempW^2];
    % 有効な観測がない時は例外処理
   elseif isempty(y)
       X = X_star;
       P = P_bar;
   else
%         K = P_bar * Hm.'/(Hm*P_bar*Hm.' + Rm);
        K = P_bar * Hm.'* pinv(Hm*P_bar*Hm.' + Rm);
        x = K * y;
        % XとPを計算
        X = X_star +x;
%         P = (eye(7) - K * Hm)*P_bar;
        P = (eye(7) - K * Hm)*P_bar * (eye(7) - K * Hm).' + K * Rm * K.';
   end
    
    obj.X = X;
    obj.P = P;
    
    obj.estNoUse = estNoUseList;

end
