classdef observeState < handle
    properties
          t                       %時刻
          length
          azimuth
          elevation
          ltd                     % light time delay
          clockError
          clockErrorLog
          xve                     %伝搬時間前の地球の位置速度
          xvg                     %伝搬時間前の地上局の位置速度
    end
    methods 
        function obj = observeState(time) %コンストラクタ 
                     obj.t = time.list;
                     obj.length = zeros(1,length(time.list));
                     obj.azimuth = zeros(1,length(time.list));
                     obj.elevation = zeros(1,length(time.list));
                     obj.ltd = zeros(1,length(time.list));
                     obj.clockError = zeros(1,length(time.list));
                     obj.xve = zeros(6,length(time.list));
                     obj.xvg = zeros(6,length(time.list));
        end   
        fromSc2Gs(obj,time,earthState,gsState,scState,constant,error,i) % 探査機から地上局を観測する
        calcClockOffset(obj,time,earthState,gsState,scState,constant,i,clockErrorCorrection)
        calcObservedValue(obj,time,ephemData,scState,constant,error,i)
    end
     methods(Static)
        Y = calcG(X_bar,xve,xvg,constant)
        H_childa = delGdelX(X_bar,xve,xvg,constant)
     end
    
end