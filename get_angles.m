function angles = get_angles(v,f)
i = 1:size(f,1);
v1 = v(f(i,1),:);
v2 = v(f(i,2),:);
v3 = v(f(i,3),:);

p1 = v2 - v1;
q1 = v3 - v1;
p2 = -p1; %v1 - v2;
q2 = v3 - v2;
p3 = -q2; %v2 - v3;
q3 = -q1; %v1 - v3;

% the following looped code is slow (factor > 30), but useful for debugging
%angles = zeros(size(f));
%for i = 1:size(f,1)
%    angles(i,1) = atan2(norm(cross(p1(i,:),q1(i,:))),dot(p1(i,:),q1(i,:)));
%     angles(i,2) = atan2(norm(cross(p2(i,:),q2(i,:))),dot(p2(i,:),q2(i,:)));
%     angles(i,3) = atan2(norm(cross(p3(i,:),q3(i,:))),dot(p3(i,:),q3(i,:)));
%     assert(abs(sum(angles(i,:))-pi) < 0.00003);
%     assert(sum(double(angles(i,:) > pi)) == 0,'Angle larger pi: %f \n',(angles(i,angles(i,:) < pi)));
%     assert(sum(double(angles(i,:) > 0)) == 3);
% end

% ugly copy and paste but faster than creating cell and looping over pairs
angles = zeros(size(f));
cross1 = cross(p1,q1);
angles(:,1) = atan2(sqrt(sum(cross1.^2,2)),sum(p1.*q1,2));

cross2 = cross(p2,q2);
angles(:,2) = atan2(sqrt(sum(cross2.^2,2)),sum(p2.*q2,2));

cross3 = cross(p3,q3);
angles(:,3) = atan2(sqrt(sum(cross3.^2,2)),sum(p3.*q3,2));