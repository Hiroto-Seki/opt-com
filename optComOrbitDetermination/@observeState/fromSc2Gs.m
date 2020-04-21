% ‘Š‘Î˜_Œø‰Ê‚Ì€‚Í“ü‚ê‚Ä‚¢‚È‚¢‚Ì‚ÅC‚ ‚Æ‚Å“ü‚ê‚é—\’è
% ^’l‚ÌŒvZ

function fromSc2Gs(obj,time,earthState,gsState,scState,constant,clockError,i)


     % “`”À’x‰„Œë·‹–—e’l
     relTimeDelayError = 1e-4;
     timeDelayErrorTemp = 100;
     % “`”À’x‰„‚ÌŒvZ 
     timeDelayTemp = 1/constant.lightSpeed * ...
         ((earthState.pos(1,i) + gsState.pos(1,i) - scState.pos(1,i))^2 + ...
          (earthState.pos(2,i) + gsState.pos(2,i) - scState.pos(2,i))^2 +....
          (earthState.pos(3,i) + gsState.pos(3,i) - scState.pos(3,i))^2)^0.5;
     while timeDelayErrorTemp > relTimeDelayError
         % “`”ÀŠÔ‚¾‚¯‘k‚Á‚½ŠÔ‚Ì’n‹…‚ÌˆÊ’u‘¬“x‚ğ‹‚ß‚é         
         % state\‘¢‘Ì‚Ìó‘Ô—Ê‚ğŠÔ‘Ñ”ÍˆÍ“à‚ÉtimeDelayTempŠÔ‘k‚Á‚½‚ªû‚Ü‚Á‚Ä‚¢‚é
         if (time.list(i) - timeDelayTemp) > time.list(1)
             % \‘¢‘Ì‚Ìó‘Ô—Ê‚Ì“àˆê”Ô‹ß‚¢‚Ì‚à‚Ì‚ğ’T‚·
            closeTimeIndex = i - round((time.list(i) - timeDelayTemp)/time.simDt);                    % ˆê”Ô‹ß‚¢‚ªt(closeTimeIndex)
            closeTimeOffset = (time.list(i) -timeDelayTemp) - time.list(closeTimeIndex);      % ˆê”Ô‹ß‚¢‚©‚ç“`”À‚µ‚È‚¯‚ê‚Î‚¢‚¯‚È‚¢ŠÔ
             % earth‚É‚Â‚¢‚Ä“`”d’x‰„(temp)ŠÔ‘O‚Ìó‘Ô—Ê‚ğ“¾‚éD
            xve = [earthState.pos(:,closeTimeIndex);earthState.vel(:,closeTimeIndex)];
            k1e = orbitalState.twobody(xve,constant.sunMu,0);
            k2e = orbitalState.twobody(xve+0.5*closeTimeOffset*k1e,constant.sunMu,0);
            k3e = orbitalState.twobody(xve+0.5*closeTimeOffset*k2e,constant.sunMu,0);
            k4e = orbitalState.twobody(xve+closeTimeOffset*k3e,constant.sunMu,0);
            xve = xve + closeTimeOffset/6*(k1e+2*k2e+2*k3e+k4e); 
            
%             % ground station‚É‚Â‚¢‚Ä“`”d’x‰„(temp)ŠÔ‘O‚Ìó‘Ô—Ê‚ğ“¾‚éD
%             [xg,vg] = groundState.earthRotation(gsState.pos(:,closeTimeIndex), closeTimeOffset, constant);
%             xvg = [xg;vg];
         % state\‘¢‘Ì‚Ìó‘Ô—Ê‚ğŠÔ‘Ñ”ÍˆÍ“à‚ÉtimeDelayTempŠÔ‘k‚Á‚½‚ªû‚Ü‚Á‚Ä‚¢‚È‚¢
         else
            % state\‘¢‘Ì‚Ì‰Šú‚Æ‚Ì·•ª.‰Šú‚©‚ç“`”À‚µ‚È‚¯‚ê‚Î‚È‚ç‚È‚¢ŠÔ(ƒ}ƒCƒiƒX¨‹t•ûŒü“`”À‚É‚È‚é)
            goBackTime = (time.list(i) - timeDelayTemp) -  time.list(1); 
            % ‰Šú‚©‚ç‘k‚é
            xve =  [earthState.pos(:,1);earthState.vel(:,1)];
            % time.simDt‚¸‚Â‹t•ûŒü“`”À‚·‚éD
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
         % ground sation‚Ì“`”ÀŠÔ‘O‚ÌˆÊ’u‘¬“x‚ğ‹‚ß‚é
         [xg,vg] = groundState.earthRotation(gsState.pos(:,i), -timeDelayTemp, constant);
         xvg = [xg;vg];
         % RLT€‚ğŠÜ‚Ü‚È‚¢“`”À’x‰„
         timeDelayNew = 1/constant.lightSpeed * ...
         ((xve(1) + xvg(1) - scState.pos(1,i))^2 + ...
          (xve(2) + xvg(2) - scState.pos(2,i))^2 +....
          (xve(3) + xvg(3)- scState.pos(3,i))^2)^0.5;
         timeDelayErrorTemp = abs(timeDelayNew - timeDelayTemp);
         timeDelayTemp = timeDelayNew;
     end
     obj.ltd(i) = timeDelayTemp + clockError; %ŒvŒë·‚Ì•ª‚¾‚¯ŠÏ‘ª—Ê‚Í‚¸‚ê‚é
     obj.length(i) = obj.ltd(i) * constant.lightSpeed;
    % Šp“x‚àŒvZ
    obj.azimuth(i) = atan2(xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed...
                                , xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed);
    obj.elevation(i) = atan( (xve(3) + xvg(3) - scState.pos(3,i) + scState.vel(3,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)/...
                                   ((xve(2) + xvg(2) - scState.pos(2,i) + scState.vel(2,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)^2 ...
                                  + (xve(1) + xvg(1) - scState.pos(1,i) + scState.vel(1,i)/constant.lightSpeed *timeDelayTemp*constant.lightSpeed)^2)^0.5 );
    
%  end
end