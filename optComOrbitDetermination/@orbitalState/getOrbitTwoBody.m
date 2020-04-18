% �T�v: time.list���̒T���@�̈ʒu��`�����ċ��߂Ă���
function getOrbitTwoBody(obj,time,constant)
    mu = constant.sunMu;
    %�܂���obj��������
    obj.t = time.list;
    obj.pos = zeros(3,length(time.list));
    obj.vel = zeros(3,length(time.list));
    % �����l�̓���
    obj.pos(:,time.refId) = reshape(obj.pos0,[3,1]);
    obj.vel(:,time.refId) = reshape(obj.vel0,[3,1]);
    
    %���`��
    xv = [obj.pos(:,time.refId);obj.vel(:,time.refId)] ;
    timestep = time.simDt;
    for j = time.refId+1:length(time.list)
        k1 = obj.twobody(xv,mu);
        k2 = obj.twobody(xv+0.5*timestep*k1,mu);
        k3 = obj.twobody(xv+0.5*timestep*k2,mu);
        k4 = obj.twobody(xv+timestep*k3,mu);
        xv = xv + timestep/6*(k1+2*k2+2*k3+k4); 
        obj.pos(:,j) = xv(1:3);
        obj.vel(:,j) = xv(4:6);
    end
    %�t�`��
    xv =[obj.pos(:,time.refId);obj.vel(:,time.refId)];
    timestep = -time.simDt;
    for j = time.refId-1 : -1:  1
        k1 = obj.twobody(xv,mu);
        k2 = obj.twobody(xv+0.5*timestep*k1,mu);
        k3 = obj.twobody(xv+0.5*timestep*k2,mu);
        k4 = obj.twobody(xv+timestep*k3,mu);
        xv = xv + timestep/6*(k1+2*k2+2*k3+k4); 
        obj.pos(:,j) = xv(1:3);
        obj.vel(:,j) = xv(4:6);
    end
end



