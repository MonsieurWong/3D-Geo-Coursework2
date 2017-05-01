function L = build_LB_operator_facebased(v,f)

%load mesh data
NumF = size(f,1);
NumV = size(v,1);
angles = get_angles(v,f);

%allocate space
L = speye(NumV);
area_face = zeros(NumF,1);
area_vertex = zeros(NumV,1);


for i = 1:NumF
    %find 3 vertices and its index
    v1 = v(f(i,1),:);
    v2 = v(f(i,2),:);
    v3 = v(f(i,3),:);
    %compute the area
    p = v1 - v2;
    q = v3 - v2;
    area_face(i) = norm(cross(p,q))/2;
    %assign the area to each vertex
    area_vertex(f(i,1)) = area_vertex(f(i,1)) + area_face(i)/3;
    area_vertex(f(i,2)) = area_vertex(f(i,2)) + area_face(i)/3;
    area_vertex(f(i,3)) = area_vertex(f(i,3)) + area_face(i)/3;
    %assign the cotan value to the entries of L
    for j = 1:3
        n1_idx = f(i,j);
        n2_idx = f(i,rem(j,3)+1);
        n3_idx = f(i,rem(rem(j,3)+1,3)+1);
        L(n1_idx,n2_idx) = cot(angles(n3_idx))/2;
        L(n2_idx,n1_idx) = cot(angles(n3_idx))/2;
%         sum_cot(n1_idx) = sum_cot(n1_idx) + cot(angles(n3_idx))/2;
%         sum_cot(n2_idx) = sum_cot(n2_idx) + cot(angles(n3_idx))/2;
    end
end
% L = full(L);
%for diagonal entries, take the negtive sum of the other entries
for i = 1:NumV  
    L(i,:) = L(i,:) / area_vertex(i);
    L(i,i) = 0;
    L(i,i) = -sum(L(i,:));
end

end