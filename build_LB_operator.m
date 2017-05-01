function L = build_LB_operator(vNum,v,f,e)

L = -speye(vNum);
L = full(L);
H = zeros(vNum,1);
area = zeros(vNum,1);
for i = 1:vNum
    % for each vertex, find its 1-ring neighbors.
    p = v(i,:);   
    [rowf colf ~] = find(f == i);
    [row col ~] = find(e == i);
    nvalence = size(row,1);
    neighbor1 = [];
    face = [];
    edge = [];
    for j = 1:nvalence
        face(j,:) = f(rowf(j),:);
        edge(j,:) = e(row(j),:);
        neighbor1(j) = e(row(j),-col(j)+3);
    end
    
    %for each neighbor, calculate the angle and cot value
    v_last = 0;
    cij = zeros(nvalence,1);
    sum_cot = 0;
    for n = 1:nvalence
        vx = neighbor1(n);
        v_x = v(vx,:);
        [vx_f, ~, ~] = find(face == vx);
        nearest = [face(vx_f(1),:) face(vx_f(2),:)];
        nearest(find(nearest == vx)) = [];
        nearest(find(nearest == i)) = [];
        v_last = v(nearest(1),:);
        v_next = v(nearest(2),:);
        p1 = v_x - v_last;
        q1 = p - v_last;
        p2 = v_x - v_next;
        q2 = p - v_next;
        %calculate the angle and barycentric area
        angle = atan2(norm(cross(p1,p2')),dot(p1,p2'));
        
        alpha = atan2(norm(cross(q1,p1')),dot(q1,p1'));
        beta = atan2(norm(cross(q2,p2')),dot(q2,p2'));
        theta = atan2(norm(cross(q1,q2')),dot(q1,q2'));
        phi = atan2(norm(cross(p1,p2')),dot(p1,p2'));
        area(i) = area(i) + sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sin(theta) / 6 + (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)) * sin(phi) / 6;
        
        cij(n) =  cot(alpha) + cot(beta);
        sum_cot = sum_cot + cij(n);
    end
    area(i) = area(i)*2;
    %construct the Laplacian-Beltrami operator
    for n = 1:nvalence
            L(i,neighbor1(n)) = cij(n)/area(i);
    end
    L(i,i) = L(i,i) * sum_cot / area(i);

end
L = L;