function [vo,fo] = triangle_split(v,f)
%set angle thresh
angle_thresh = 110*pi/180;
angles = get_angles(v,f);
idx = 1;

while(idx <= size(f,1))
    %decide whether to split or not
    split_angle = double(angles(idx,:)>= angle_thresh);
    triangle_split = (sum(split_angle) > 0);
    
    %split triangles
    if triangle_split
        %find the angle to be splited
        col = sum([1 2 3].*split_angle);
        %take 3 vertices to derive the perpendicular point d
        idx_v1 = f(idx,col);
        idx_v2 = f(idx,rem(col,3)+1);
        idx_v3 = f(idx,rem(rem(col,3)+1,3)+1);
        toc;
        v1 = v(idx_v1,:);
        v2 = v(idx_v2,:);
        v3 = v(idx_v3,:);
        p = (v3 - v2)/norm(v3 - v2);
        q = (v1 - v2)/norm(v1 - v2);
        tau = cross(p,cross(p,q));
        tau = tau/norm(tau);
        d = v1 + tau * dot(tau,(v2 - v1));
        %adding new vertex and 2 new faces, delete the splited faces
        idx_d = size(v,1) + 1;
        f(idx,:) = [];
        v(idx_d,:) = d;
        new_faces = [idx_v1 idx_v2 idx_d ; idx_v1 idx_v3 idx_d];
        f = [f ; new_faces];
        % update angle matrix
        angles(idx,:) = [];
        angles = [angles;get_angles(v,new_faces)];
        idx = idx-1;
    end
    idx = idx+1;    
end
vo = v;
fo = f;