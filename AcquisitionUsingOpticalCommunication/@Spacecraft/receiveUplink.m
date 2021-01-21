% 宇宙機がuplinkから受信する内容を求める
% obj = scTrue
function [obj,gsTrue] = receiveUplink(obj,gsTrue,earth,constant,time) 
    %% gsTrue.ut_counter番目のuplinkを宇宙機が受信する時刻と状態量を求める
    ut_counter = gsTrue.ut_counter;
    %% 
    
    [obj.t_ur(ut_counter),obj.state_ur(:,ut_counter),dtlt] ...
        = GroundStation.calcTarget(gsTrue.t_ut(ut_counter),gsTrue.state_ut(:,ut_counter),earth.state_ut(:,ut_counter),[],obj,time,constant,"true value"); 
%     % 確認する
%     distance_t = (obj.t_ur(ut_counter) - gsTrue.t_ut(ut_counter))*constant.lightSpeed;
%     distance_r = norm(gsTrue.state_ut(1:3,ut_counter) + earth.state_ut(1:3,ut_counter) - obj.state_ur(1:3,ut_counter));
%     disp(distance_t - distance_r )
    
    %% 受信する内容
    obj.gsT_ur(ut_counter)               =  gsTrue.t_ut(ut_counter);                 %uplinkに載っている，送信時刻
    obj.eState_ur(:,ut_counter)          =  earth.state_ut(:,ut_counter);            %uplinkに載っている，送信時刻の地球の状態量   
    obj.gsState_ur(:,ut_counter)         =  gsTrue.state_ut(:,ut_counter);           %uplinkに載っている, 送信時刻の地上局の状態量
    obj.transDirection_ur(:,ut_counter)  =  gsTrue.direction_ut(:,ut_counter);       %uplinkに載っている，送信方向(誤差を持つ)
    obj.lengthTrue_ur(ut_counter)        =  dtlt*constant.lightSpeed;                %クロック誤差が乗っていないいない測距情報
    
    
    
    
    
    %% 送信機のpointing誤差の計算
    trueDirection_ut = obj.state_ur(:,ut_counter) - obj.eState_ur(:,ut_counter) - obj.gsState_ur(:,ut_counter);
    
    gsTrue.pointingError_ut(ut_counter) = ((atan2(trueDirection_ut(2),trueDirection_ut(1)) - obj.transDirection_ur(1,ut_counter))^2 ...
                                                + (atan(trueDirection_ut(3) /(trueDirection_ut(1)^2 + trueDirection_ut(2)^2)^0.5) - obj.transDirection_ur(2,ut_counter))^2)^0.5;
    
    
end