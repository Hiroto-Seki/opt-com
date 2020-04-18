% time.list‚É‘Î‰‚·‚é‚Ì“V‘Ì‚ÌˆÊ’u‘¬“x‚ğ“¾‚é
function getEphemData(obj,time,error)
 %‚Ü‚¸‚Íobj‚ğ‰Šú‰»
    obj.t = time.list;
    obj.pos = zeros(3,length(time.list));
    obj.vel = zeros(3,length(time.list));
    if strcmp(obj.dataType, 'true') || strcmp(obj.dataType, 'True')
        offset = error.time0;
    else
        offset = 0;
    end
   for i = 1:length(time.list)
       day = time.t0 + time.list(i)/60/60/24 + offset;
       [pos,vel] = planetEphemeris(day,obj.centralBody,obj.name);
       obj.pos(:,i) = reshape(pos,[3, 1]);
       obj.vel(:,i) = reshape(vel,[3, 1]);
   end
end