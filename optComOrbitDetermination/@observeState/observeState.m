classdef observeState < handle
    properties
          t                       %����
          length
          azimuth
          elevation
          ltd                     % light time delay
          clockError
          clockErrorLog
    end
    methods 
        function obj = observeState(time) %�R���X�g���N�^ 
                     obj.t = time.list;
                     obj.length = zeros(1,length(time.list));
                     obj.azimuth = zeros(1,length(time.list));
                     obj.elevation = zeros(1,length(time.list));
                     obj.ltd = zeros(1,length(time.list));
                     obj.clockError = zeros(1,length(time.list));
        end   
        fromSc2Gs(obj,time,earthState,gsState,scState,constant,error,i) % �T���@����n��ǂ��ϑ�����
        calcClockOffset(obj,time,earthState,gsState,scState,constant,i,clockErrorCorrection)
    end
end