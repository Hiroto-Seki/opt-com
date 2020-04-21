classdef observeState < handle
    properties
          t                       %����
          length
          azimuth
          elevation
          ltd                     % light time delay
          clockError
          clockErrorLog
          xve                     %�`�����ԑO�̒n���̈ʒu���x
          xvg                     %�`�����ԑO�̒n��ǂ̈ʒu���x
    end
    methods 
        function obj = observeState(time) %�R���X�g���N�^ 
                     obj.t = time.list;
                     obj.length = zeros(1,length(time.list));
                     obj.azimuth = zeros(1,length(time.list));
                     obj.elevation = zeros(1,length(time.list));
                     obj.ltd = zeros(1,length(time.list));
                     obj.clockError = zeros(1,length(time.list));
                     obj.xve = zeros(6,length(time.list));
                     obj.xvg = zeros(6,length(time.list));
        end   
        fromSc2Gs(obj,time,earthState,gsState,scState,constant,error,i) % �T���@����n��ǂ��ϑ�����
        calcClockOffset(obj,time,earthState,gsState,scState,constant,i,clockErrorCorrection)
        calcObservedValue(obj,time,ephemData,scState,constant,error,i)
    end
     methods(Static)
        Y = calcG(X_bar,xve,xvg,constant)
        H_childa = delGdelX(X_bar,xve,xvg,constant)
     end
    
end