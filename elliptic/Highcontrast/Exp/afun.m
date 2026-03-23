% function a=afun(t,s)
% if abs(t-0.4)<0.015
%     a=10^4;
% else
%     a=1;
% end

function a=afun(t,s)
    pts1=0.2:0.1:0.8;
    pts2=0.2:0.1:0.8;
    [pts_x,pts_y]=meshgrid(pts1,pts2);
    distance=sqrt((pts_x-t).^2+(pts_y-s).^2);
    if min(distance)>0.015
        a=1;
    else
        a=2^6;
    end
end
