% �n��ǁ��T���@��one-way�̒ʐM�ŁC���v�덷�����ς���

function calcClockOffset(obj,time,earthState,gsState,scState,constant,i,clockErrorCorrection)
relDt = 1;
relTolDt = 0.0001;
iterNum = 1;
while (relDt > relTolDt) && (iterNum < 100)
% for i = 1:timeNumber
    t3 = time.list(i) ;
    t2 = t3 - obj.ltd(i) - clockErrorCorrection;
    %%  ����t2�ł̐���l�̈ʒu, ���x�����߂� 
    if t2 > time.list(1)
            % �\���̂̏�ԗʂ̓���ԋ߂������̂��̂�T��
            closeTimeIndex = i - round(t2/time.simDt);                    % ��ԋ߂�������t(closeTimeIndex)
            closeTimeOffset = t2 - time.list(closeTimeIndex);      % ��ԋ߂���������`�����Ȃ���΂����Ȃ�����
             % earth�ɂ��ē`�d�x��(temp)���ԑO�̏�ԗʂ𓾂�D
            xve = [earthState.pos(:,closeTimeIndex);earthState.vel(:,closeTimeIndex)];
            k1e = orbitalState.twobody(xve,constant.sunMu);
            k2e = orbitalState.twobody(xve+0.5*closeTimeOffset*k1e,constant.sunMu);
            k3e = orbitalState.twobody(xve+0.5*closeTimeOffset*k2e,constant.sunMu);
            k4e = orbitalState.twobody(xve+closeTimeOffset*k3e,constant.sunMu);
            xve = xve + closeTimeOffset/6*(k1e+2*k2e+2*k3e+k4e); 
         % state�\���̂̏�ԗʂ����ԑє͈͓���timeDelayTemp���ԑk�������������܂��Ă��Ȃ���
     else
            % state�\���̂̏��������Ƃ̍���.������������`�����Ȃ���΂Ȃ�Ȃ�����(�}�C�i�X���t�����`���ɂȂ�)
            goBackTime = t2 -  time.list(1); 
            % ������������k��
            xve =  [earthState.pos(:,1);earthState.vel(:,1)];
            % time.simDt���t�����`������D
            while goBackTime < 0
                if abs(goBackTime) > time.simDt
                    backTimeStep = -time.simDt;
                else
                    backTimeStep = goBackTime;
                end
                k1e = orbitalState.twobody(xve,constant.sunMu);
                k2e = orbitalState.twobody(xve+0.5*backTimeStep*k1e,constant.sunMu);
                k3e = orbitalState.twobody(xve+0.5*backTimeStep*k2e,constant.sunMu);
                k4e = orbitalState.twobody(xve+backTimeStep*k3e,constant.sunMu);
                xve = xve + backTimeStep/6*(k1e+2*k2e+2*k3e+k4e);
                goBackTime = goBackTime + time.simDt;
            end           
     end
     % ground sation�̓`�����ԑO�̈ʒu���x�����߂�
     [xg,rg] = groundState.earthRotation(gsState.pos(:,i), -obj.ltd(i), constant);
     
     % ����t3�̉F���@�̎���t2�̒n��ǂɑ΂��鑊�Έʒu�����߂�
     r2 = xve(1:3) + xg;
     v2 = xve(4:6) + rg;
     r23 = scState.pos(:,i) - r2;
%      r2r2Dot = r2.'*v2;
%      r3r3Dot = scState.pos(:,i).' * scState.vel(:,i);
     pDot23 = r23.' * v2/norm(r23);    
    % ���v�덷�����܂�
     obj.clockError(i) = - (t3 - t2 - norm(r23)/constant.lightSpeed)/(1-pDot23/constant.lightSpeed) + clockErrorCorrection;
     relDt = abs(obj.clockError(i) - clockErrorCorrection);
     clockErrorCorrection = obj.clockError(i);
     obj.clockErrorLog(i,iterNum) = clockErrorCorrection; 
     iterNum = iterNum +1;
end


end