% ���Θ_���ʂ̍��͓���Ă��Ȃ��̂ŁC���Ƃœ����\��
% �^�l�̌v�Z

function fromSc2Gs(obj,time,earthState,gsState,scState,constant,clockError,i)


     % �`���x���덷���e�l
     relTimeDelayError = 1e-4;
     timeDelayErrorTemp = 100;
     % �`���x���̌v�Z 
     timeDelayTemp = 1/constant.lightSpeed * ...
         ((earthState.pos(1,i) + gsState.pos(1,i) - scState.pos(1,i))^2 + ...
          (earthState.pos(2,i) + gsState.pos(2,i) - scState.pos(2,i))^2 +....
          (earthState.pos(3,i) + gsState.pos(3,i) - scState.pos(3,i))^2)^0.5;
     while timeDelayErrorTemp > relTimeDelayError
         % �`�����Ԃ����k�������Ԃ̒n���̈ʒu���x�����߂�         
         % state�\���̂̏�ԗʂ����ԑє͈͓���timeDelayTemp���ԑk�������������܂��Ă��鎞
         if (time.list(i) - timeDelayTemp) > time.list(1)
             % �\���̂̏�ԗʂ̓���ԋ߂������̂��̂�T��
            closeTimeIndex = i - round((time.list(i) - timeDelayTemp)/time.simDt);                    % ��ԋ߂�������t(closeTimeIndex)
            closeTimeOffset = (time.list(i) -timeDelayTemp) - time.list(closeTimeIndex);      % ��ԋ߂���������`�����Ȃ���΂����Ȃ�����
             % earth�ɂ��ē`�d�x��(temp)���ԑO�̏�ԗʂ𓾂�D
            xve = [earthState.pos(:,closeTimeIndex);earthState.vel(:,closeTimeIndex)];
            k1e = orbitalState.twobody(xve,constant.sunMu,0);
            k2e = orbitalState.twobody(xve+0.5*closeTimeOffset*k1e,constant.sunMu,0);
            k3e = orbitalState.twobody(xve+0.5*closeTimeOffset*k2e,constant.sunMu,0);
            k4e = orbitalState.twobody(xve+closeTimeOffset*k3e,constant.sunMu,0);
            xve = xve + closeTimeOffset/6*(k1e+2*k2e+2*k3e+k4e); 
            
%             % ground station�ɂ��ē`�d�x��(temp)���ԑO�̏�ԗʂ𓾂�D
%             [xg,vg] = groundState.earthRotation(gsState.pos(:,closeTimeIndex), closeTimeOffset, constant);
%             xvg = [xg;vg];
         % state�\���̂̏�ԗʂ����ԑє͈͓���timeDelayTemp���ԑk�������������܂��Ă��Ȃ���
         else
            % state�\���̂̏��������Ƃ̍���.������������`�����Ȃ���΂Ȃ�Ȃ�����(�}�C�i�X���t�����`���ɂȂ�)
            goBackTime = (time.list(i) - timeDelayTemp) -  time.list(1); 
            % ������������k��
            xve =  [earthState.pos(:,1);earthState.vel(:,1)];
            % time.simDt���t�����`������D
            while goBackTime < 0
                if abs(goBackTime) > time.simDt
                    backTimeStep = -time.simDt;
                else
                    backTimeStep = goBackTime;
                end
                k1e = orbitalState.twobody(xve,constant.sunMu,0);
                k2e = orbitalState.twobody(xve+0.5*backTimeStep*k1e,constant.sunMu,0);
                k3e = orbitalState.twobody(xve+0.5*backTimeStep*k2e,constant.sunMu,0);
                k4e = orbitalState.twobody(xve+backTimeStep*k3e,constant.sunMu,0);
                xve = xve + backTimeStep/6*(k1e+2*k2e+2*k3e+k4e);
                goBackTime = goBackTime + time.simDt;
            end           
         end
         % ground sation�̓`�����ԑO�̈ʒu���x�����߂�
         [xg,vg] = groundState.earthRotation(gsState.pos(:,i), -timeDelayTemp, constant);
         xvg = [xg;vg];
         % RLT�����܂܂Ȃ��`���x��
         timeDelayNew = 1/constant.lightSpeed * ...
         ((xve(1) + xvg(1) - scState.pos(1,i))^2 + ...
          (xve(2) + xvg(2) - scState.pos(2,i))^2 +....
          (xve(3) + xvg(3)- scState.pos(3,i))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     obj.ltd(i) = timeDelayTemp + clockError; %���v�덷�̕������ϑ��ʂ͂����
     obj.length(i) = obj.ltd(i) * constant.lightSpeed;
    % �p�x���v�Z
    obj.azimuth(i) = atan2(xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed...
                                , xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed);
    obj.elevation(i) = atan( (xve(3) + xvg(3) - scState.pos(3,i) + scState.vel(3,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)/...
                                   ((xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)^2 ...
                                  + (xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)^2)^0.5 );
    
%  end
end