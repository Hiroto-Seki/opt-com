
function observationUpdateByScUkf(obj,scTrue,earth, gsTrue,constant,type,ukf,time)
    % 何度目の観測か
    ur_counter = scTrue.ur_counter;
    % 何度目のダウンリンクを2wayに使っているか(0の場合は2wayはできていない)
    ur2w_counter = scTrue.ur2w_counter;
    % シグマ点列を取得
    x_sp   = obj.x_sp;
    % シグマ点列の重み平均を取得
    x_mean = obj.X;
    
    % 観測方程式に出てくるもの
    xve_ut = earth.state_ut(:,ur_counter);
    xvg_ut = gsTrue.state_ut(:,ur_counter);
    
    % Rの書き換え
    % 1,2要素目の観測は，1wayでも2wayでもuplinkの測角になっている．
    obj.R.direction_ur = scTrue.directionAccuracy_ur(ur_counter)^2;
    obj.R.direction_ut = gsTrue.directionAccuracy_ut(ur_counter)^2;
    
    
    if type ==2
        xve_dr = earth.state_dr(:,ur2w_counter);
        xvg_dr = gsTrue.state_dr(:,ur2w_counter);
        dtAtGs = scTrue.durationAtGs(ur_counter);
        dt2w  = scTrue.t_ur(ur_counter) - scTrue.t_dt(ur2w_counter);        
    end
    %% 観測値を取得
    if type == 1 %1wayの観測
        Y.direction_ur = scTrue.directionObserved_ur(:,ur_counter);
        Y.direction_ut = scTrue.transDirection_ur(:,ur_counter);
        Y.accel_ur     = scTrue.accelObseved_ur(:,ur_counter);
        Y.length1w_ur  = scTrue.lengthObserved_ur(ur_counter);
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
    else %2wayの観測
        Y.direction_ur = scTrue.directionObserved_ur(:,ur_counter); % 測角
        Y.direction_ut = scTrue.transDirection_ur(:,ur_counter);...     % uplink方向
        Y.accel_ur     = scTrue.accelObseved_ur(:,ur_counter);...       % 加速度計
        Y.length1w_ur  = scTrue.lengthObserved_ur(ur_counter);          % 測距1way
        Y.length2w_ur  = scTrue.length2wObserved_ur(ur_counter);        % 測距2way
        Y.direction_dr = scTrue.recDownAngle_ur(:,ur_counter);         % uplinkされる，地上局での観測
        % Rの書き換え
        obj.R.direction_dr = scTrue.recDownAngleAccuracy_ur(ur_counter)^2;
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
    end
    
    %% 観測があるかどうかで場合わけ
    [~,~,~,~,obsNum] = Spacecraft.alignReqInfo4Est([],[],[],[],obsType,"obsCheck",obj.useObs); % 観測があるかどうかだけ見ている
   if obsNum == 0
       obj.X = x_mean;
       obj.P = obj.P;
   else
    
    
        %% シグマ点列による観測ベクトルの計算
        % 初期化
        y_mean.azm_ur = 0;
        y_mean.elv_ur = 0;
        y_mean.azm_ut = 0;
        y_mean.elv_ut = 0;
        y_mean.accel_ur = zeros(3,1);
        y_mean.length1w_ur = 0;
        y_mean.length2w_ur = 0;
        y_mean.azm_dr = 0;
        y_mean.elv_dr = 0;
        for i = 1:size(x_sp,2)
            x_spi = x_sp(:,i);
            % 観測によって場合分け
            if type == 1
                y_sp(i) = Spacecraft.calcG_ur(x_spi,xve_ut,xvg_ut,[],[],[],[],constant,[],"1way");
                %ベクトル化
                [~,y_spM(:,i),~,~] = Spacecraft.alignReqInfo4Est(Y,y_sp(i),[],obj.R,obsType,"ukf",obj.useObs);
            else
                y_sp(i) = Spacecraft.calcG_ur(x_spi,xve_ut,xvg_ut,xve_dr,xvg_dr,...
                                            dtAtGs, dt2w,constant,time,"2way");
                % ベクトル化
                [~,y_spM(:,i),~,~] = Spacecraft.alignReqInfo4Est(Y,y_sp(i),[],obj.R,obsType,"ukf",obj.useObs);
            end
            if i ==1
                wm = ukf.w0_m;
            else
                wm = ukf.wi_m;
            end
            y_mean.azm_ur = y_mean.azm_ur + wm * y_sp(i).azm_ur;
            y_mean.elv_ur = y_mean.elv_ur + wm * y_sp(i).elv_ur;
            y_mean.azm_ut = y_mean.azm_ut + wm * y_sp(i).azm_ut;
            y_mean.elv_ut = y_mean.elv_ut + wm * y_sp(i).elv_ut;
            y_mean.accel_ur = y_mean.accel_ur  + wm * y_sp(i).accel_ur;
            y_mean.length1w_ur = y_mean.length1w_ur + wm * y_sp(i).length1w_ur;
            if type ==2
                y_mean.length2w_ur = y_mean.length2w_ur + wm * y_sp(i).length2w_ur;
                y_mean.azm_dr = y_mean.azm_dr + wm * y_sp(i).azm_dr;
                y_mean.elv_dr = y_mean.elv_dr + wm * y_sp(i).elv_dr;
            end
        end

        % Y, Y_meanをベクトル化，Rを行列に
        if type == 1
            [Yv,y_meanV,~,Rm] = Spacecraft.alignReqInfo4Est(Y,y_mean,[],obj.R,obsType,"ukf",obj.useObs);
        else
            [Yv,y_meanV,~,Rm] = Spacecraft.alignReqInfo4Est(Y,y_mean,[],obj.R,obsType,"ukf",obj.useObs);
        end




        %% 共分散行列及び相互共分散行列
        for j = 1:size(x_sp,2)
            if j == 1
                Pvv = ukf.w0_c * (y_spM(:,j) - y_meanV) * (y_spM(:,j) - y_meanV).';
                Pxy = ukf.w0_c * (x_sp(:,j) - x_mean) * (y_spM(:,j) - y_meanV).';
            else
                Pvv = Pvv + ukf.wi_c * (y_spM(:,j) - y_meanV) * (y_spM(:,j) - y_meanV).';
                Pxy = Pxy + ukf.wi_c * (x_sp(:,j) - x_mean) * (y_spM(:,j) - y_meanV).';           
            end
        end

        Pvv = Pvv + Rm;


        %% 観測残差及び残差検定
        v = Yv - y_meanV;

        % デバッグ用に記録
        obj.y = v;

        for k = length(v):-1:1
            if ukf.sigmaN < abs(v(k))/sqrt(Pvv(k,k))
                % 観測を棄却する
                v(k) = [];
                Pvv(k,:) = [];
                Pvv(:,k) = [];
                Pxy(:,k) = [];
            end
        end

        %有効な観測がない時は例外処理
       if isempty(v)
           obj.X = x_mean;
           obj.P = obj.P;
       else
            %% カルマンゲインの計算
            K = Pxy /Pvv;

            %% 状態ベクトルと共分散行列の観測更新
            obj.X = x_mean + K * v;
            obj.P = obj.P - K *Pxy.';
       end      
   end

end