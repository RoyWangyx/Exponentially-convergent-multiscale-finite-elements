function [f] =vfun(x,y)
    pts1=0.2:0.1:0.8;
    pts2=0.2:0.1:0.8;
    [pts_x,pts_y]=meshgrid(pts1,pts2);
    distance=max(abs(pts_x-x),abs(pts_y-y));
    if min(min(distance))>0.025
        f=1;
    else
        f=sqrt(2);
    end
end