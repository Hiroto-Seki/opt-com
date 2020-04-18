classdef orbitalState < handle
    properties
          name                    %name
          dataType                %��ԗʂ̎��(����l, �^�l)
          t                       %����
          centralBody             %�_�C�i�~�N�X�̒��S�V��
%           a                       %semi-major axis           (km)
%           e                       %Eccentricity Magnitude
%           i                       %inclination[rad]
%           AN                      %ascending node            (rad)
%           AP                      %argument of periapses     (rad)
%           f                       %true anomaly angle (rad)
          mu                      %Gravitational Constant of central body
          pos                     %position realtive to the central body [km]
          vel                     %velocity realtive to the central body [km]
          pos0
          vel0
    end
    methods 
        function obj = orbitalState(name,dataType,centralBody) %�R���X�g���N�^ 
            obj.name = name;
            obj.dataType = dataType;
            obj.centralBody = centralBody;
                     
        end              
        elem2rv(obj);
        getOrbitTwoBody(obj,time,constant)
        getEphemData(obj,time,error)
    end
    methods(Static)
        xv = twobody(xv, mu)
    end
end