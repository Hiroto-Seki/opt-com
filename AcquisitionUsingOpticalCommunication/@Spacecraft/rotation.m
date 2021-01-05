function rotationMatrix = rotation(roll,pitch,yaw,order) 
    %X軸周り(roll)の回転行列
    Rx = [1,         0,         0;
          0, cos(roll),-sin(roll);
          0, sin(roll), cos(roll)];
    %Y軸周り(pitch)の回転行列
    Ry = [ cos(pitch),         0, sin(pitch);
                    0,         1,          0;
          -sin(pitch),         0, cos(pitch)];
    %Z軸周りの
    Rz = [   cos(yaw), -sin(yaw),          0;
             sin(yaw),  cos(yaw),          0;
                    0,         0,          1];
                
    % 回転の順番によって場合分け
    if order == 1 %X->Y->Zの順番に回転する場合
        rotationMatrix = Rz * Ry * Rx;
    elseif order ==2  %Z->Y->Xの順番に回転する場合
        rotationMatrix = Rx * Ry * Rz;        
    end
end