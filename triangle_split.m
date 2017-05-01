function [vo,fo] = triangle_split(v,f)

%set angle thresh
angle_thresh = 110*pi/180;
angles = get_angles(v,f);

for i = 1:size(f,1)
    %decide whether to split or not
    triangle_split = false;
    split_angle = false(3,1);
    for j = 1:3
        split_angle(j) = (angles(i,j)>=angle_thresh) ;
        triangle_split = triangle_split || split_angle(j);
    end
    
    %split triangles
    if triangle_split
        %find the angle to be splited
        col = find(split_angle>0);
        %take 3 vertices to derive the perpendicular point d
        idx_v1 = f(i,col);
        idx_v2 = f(i,rem(col,3)+1);
        idx_v3 = f(i,rem(rem(col,3)+1,3)+1);
        v1 = v(idx_v1,:);
        v2 = v(idx_v2,:);
        v3 = v(idx_v3,:);
        p = (v3 - v2)/norm(v3 - v2);
        q = (v1 - v2)/norm(v1 - v2);
        tao = cross(p,cross(p,q));
        tao = tao/norm(tao);
        d = v1 + tao * dot(tao,(v2 - v1));
        %adding new vertex and 2 new faces, delete the splited faces
        idx_d = size(v,1) + 1;
        f(i,:) = [];
        v(idx_d,:) = d;
        f = [f ; idx_v1 idx_v2 idx_d ; idx_v1 idx_v3 idx_d];
    end    
end
vo = v;
fo = f;