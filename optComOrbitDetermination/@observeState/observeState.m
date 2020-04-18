classdef observeState < handle
    properties
          t                       %時刻
          length
          azimuth
          elevation
          ltd                     % light time delay
          clockError
          clockErrorLog
    end
    methods 
        function obj = observeState(time) %コンストラクタ 
                     obj.t = time.list;
                     obj.length = zeros(1,length(time.list));
                     obj.azimuth = zeros(1,length(time.list));
                     obj.elevation = zeros(1,length(time.list));
                     obj.ltd = zeros(1,length(time.list));
                     obj.clockError = zeros(1,length(time.list));
        end   
        fromSc2Gs(obj,time,earthState,gsState,scState,constant,error,i) % 探査機から地上局を観測する
        calcClockOffset(obj,time,earthState,gsState,scState,constant,i,clockErrorCorrection)
    end
end