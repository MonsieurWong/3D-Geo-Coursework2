function angles = get_angles(v,f)
i = 1:size(f,1);
v1 = v(f(i,1),:);
v2 = v(f(i,2),:);
v3 = v(f(i,3),:);
p1 = v2 - v1;
q1 = v3 - v1;
p2 = v1 - v2;
q2 = v3 - v2;
p3 = v2 - v3;
q3 = v1 - v3;
% fun = @(e1,e2) atan2(norm(cross(e1,e2)),dot(e1,e2));
angles = zeros(size(f));
for i = 1:size(f,1)
    angles(i,1) = atan2(norm(cross(p1(i,:),q1(i,:))),dot(p1(i,:),q1(i,:)));
    angles(i,2) = atan2(norm(cross(p2(i,:),q2(i,:))),dot(p2(i,:),q2(i,:)));
    angles(i,3) = atan2(norm(cross(p3(i,:),q3(i,:))),dot(p3(i,:),q3(i,:)));
end

