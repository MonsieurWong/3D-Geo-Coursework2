%%============================Initialization============================%%
clear all;
close all;
% [v f] = read_ply('dragon_vrip_res3.ply');
% [v f] = read_off('decimated-knight.off');
[v f] = read_off('bunny.off');
v = v'; 
f = f';

e = compute_edges(f)';
n = compute_normal(v',f')';
vNum = size(v,1);
fNum = size(f,1);

%% =====================Uniform Mean Curvature======================%% 
L1 = -speye(vNum);
%build the Laplacian matrix, for each vertex, find the nearest 1-ring
%neighbors and assign weight to each neighbor.
for i = 1:vNum
    [row col ~] = find(e == i);
    nvalence = size(row,1);
    for j = 1:nvalence
        neighbor(j) = e(row(j),-col(j)+3);
        L1(i,neighbor(j)) = 1 / nvalence;
    end
end
%apply the matrix to vetices matrix
L1 = sparse(L1);
deltaV = L1 * v;
% compute the mean curvature at each vertex
for i = 1:vNum
    H(i) = deltaV(i,:) * n(i,:)' * -0.5;
end
H = H';
hN = H;

% draw the mesh
view(2);
lighting phong;
material shiny;
axis off;
colormap jet(256);
light( 'position',[10,10,10]);

h = patch('vertices',v,'faces',f,'FaceVertexCData',hN ,'VertexNormal',n,'FaceColor','flat','LineStyle','none');
pause;

%% ======================Gaussian Curvature=========================%%
%intialization
H_G = zeros(vNum,1);
area = zeros(vNum,1);
theta = zeros(vNum,1);
phi = zeros(vNum,1);
sum_theta =zeros(vNum,1);
%for each vertex, find nearest 1-ring neighbor and calculate the
%corresponding angle and are
for i = 1:vNum
    %find 1-ring neighbors
    p = v(i,:);
    face = [];
    edge = [];
    [rowf colf ~] = find(f == i);
    [row col ~] = find(e == i);
    nvalence = size(row,1);
    neighbor1 = [];
    
    for j = 1:nvalence
        face(j,:) = f(rowf(j),:);
        edge(j,:) = e(row(j),:);
        neighbor1(j) = e(row(j),-col(j)+3);
    end
    
    %for each neighbor, calculate the area and angle
    for n = 1:nvalence
        %find the neighboring edges
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
        %calculate the angles and area
        theta(i) = acos((q1) * (q2)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2)));
        phi(i) = acos((p1) * (p2)' / (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
        area(i) = area(i) + sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sin(theta(i)) / 6 + (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)) * sin(phi(i)) / 2  * 1/3;
        sum_theta(i) = sum_theta(i) + theta(i)/2;
    end
    %calculate the Gaussian curvature
    H_G(i) = (2*pi - sum_theta(i)) / area(i);
    H_G(i) = real(H_G(i));

end

%draw the mesh
figure;
view(2);
lighting phong;
material shiny
axis off;
colormap jet(256);
light( 'position',[10,5,10]);

h = patch('vertices',v,'faces',f,'FaceVertexCData',H_G ,'FaceColor','flat','LineStyle','none');
pause;


%% ================Discrete L-B operator Mean Curvature=============%%
% L1 = -speye(vNum);
L1 = eye(vNum);
H = zeros(vNum,1);
area = zeros(vNum,1);
%for each vertex, find nearest 1-ring neighbor and calculate the
%corresponding angle and cot value
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
        alpha = acos((q1) * (p1)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2)));
        beta = acos((q2) * (p2)' / (sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
        theta = acos((q1) * (q2)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2)));
        phi = acos((p1) * (p2)' / (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
        area(i) = area(i) + sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sin(theta) / 6 + (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)) * sin(phi) / 6;
        cij(n) =  cot(alpha) + cot(beta);
        sum_cot = sum_cot + cij(n); 
    end
    %construct the Laplacian-Beltrami operator
    for n = 1:nvalence
            L(i,neighbor1(n)) = cij(n)/area(i);
    end
    L(i,i) = L(i,i) * sum_cot / area(i);

end
%apply the L-B operator to vertices matrix
deltaX = L * v / 2;
for i = 1:vNum
    H(i) = sqrt(deltaX(i,1)^2 + deltaX(i,2)^2 + deltaX(i,3)^2);
end

%draw the mesh
figure;
view(2);
lighting phong;
material shiny;
axis off;
colormap jet(256);
light( 'position',[10,10,10]);

h = patch('vertices',v,'faces',f,'FaceVertexCData',H ,'FaceColor','flat','LineStyle','none');
pause;


%% =========================Explicit Smoothing=====================%%
ve = v;
for I = 1:15
%     L = -eye(vNum);
%     area = zeros(vNum,1);
%     for i = 1:vNum
%         p = ve(i,:);   
%         [rowf colf ~] = find(f == i);
%         [row col ~] = find(e == i);
%         nvalence = size(row,1);
%         neighbor1 = [];
%         face = [];
%         edge = [];
%         for j = 1:nvalence
%             face(j,:) = f(rowf(j),:);
%             edge(j,:) = e(row(j),:);
%             neighbor1(j) = e(row(j),-col(j)+3);
%         end
% 
%         v_last = 0;
%         cij = zeros(nvalence,1);
%         sum_cot = 0;
%         for n = 1:nvalence
%             vx = neighbor1(n);
%             v_x = ve(vx,:);
%             [vx_f, ~, ~] = find(face == vx);
%             nearest = [face(vx_f(1),:) face(vx_f(2),:)];
%             nearest(find(nearest == vx)) = [];
%             nearest(find(nearest == i)) = [];
%             v_last = ve(nearest(1),:);
%             v_next = ve(nearest(2),:);
%             p1 = v_x - v_last;
%             q1 = p - v_last;
%             p2 = v_x - v_next;
%             q2 = p - v_next;
%             alpha = acos((q1) * (p1)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2)));
%             beta = acos((q2) * (p2)' / (sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
%             theta = acos((q1) * (q2)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2)));
%             phi = acos((p1) * (p2)' / (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
%             area(i) = area(i) + sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sin(theta) / 6 + (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)) * sin(phi) / 6;
%             cij(n) =  cot(alpha) + cot(beta);
%             sum_cot = sum_cot + cij(n); %* (v_current - p));
%         end
%         for n = 1:nvalence
%                 L(i,neighbor1(n)) = cij(n)/area(i);
%         end
%         L(i,i) = L(i,i) * sum_cot / area(i);
    L = build_LB_operator(vNum,v,f,e,n);
    
%     L = sparse(L);
    deltaX = L * ve / 2;
    for i = 1:vNum
        H(i) = sqrt(deltaX(i,1)^2 + deltaX(i,2)^2 + deltaX(i,3)^2);
    end
    lamda = 5e-7; % step for the bunny mesh
%     lamda = 1e-2; % step for the bumpy cube mesh
    IS =eye(vNum) + lamda * L;
    ve = IS * ve;

    figure;
    view(2);
    lighting phong;
    material shiny
    axis off;
    colormap jet(256);
    light( 'position',[10,10,10]);

    h = patch('vertices',ve,'faces',f,'FaceVertexCData',H ,'FaceColor','flat','LineStyle','none');
    pause(0.2);
end
pause;
close all;

%% ======================Implicit Smoothing========================%%
vi = v;
for I = 1:15
%     H = zeros(vNum,1);
%     L = -eye(vNum);
%     area = zeros(vNum,1);
%     for i = 1:vNum
%         p = vi(i,:);   
%         [rowf colf ~] = find(f == i);
%         [row col ~] = find(e == i);
%         nvalence = size(row,1);
%         neighbor1 = [];
%         face = [];
%         edge = [];
%         for j = 1:nvalence
%             face(j,:) = f(rowf(j),:);
%             edge(j,:) = e(row(j),:);
%             neighbor1(j) = e(row(j),-col(j)+3);
%         end
% 
%         v_last = 0;
%         cij = zeros(nvalence,1);
%         sum_cot = 0;
%         for n = 1:nvalence
%             vx = neighbor1(n);
%             v_x = vi(vx,:);
%             [vx_f, ~, ~] = find(face == vx);
%             nearest = [face(vx_f(1),:) face(vx_f(2),:)];
%             nearest(find(nearest == vx)) = [];
%             nearest(find(nearest == i)) = [];
%             v_last = vi(nearest(1),:);
%             v_next = vi(nearest(2),:);
%             p1 = v_x - v_last;
%             q1 = p - v_last;
%             p2 = v_x - v_next;
%             q2 = p - v_next;
%             alpha = acos((q1) * (p1)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2)));
%             beta = acos((q2) * (p2)' / (sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
%             theta = acos((q1) * (q2)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2)));
%             phi = acos((p1) * (p2)' / (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
%             area(i) = area(i) + sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sin(theta) / 6 + (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)) * sin(phi) / 6;
%             cij(n) =  cot(alpha) + cot(beta);
%             sum_cot = sum_cot + cij(n); %* (v_current - p));
%         end
%         for n = 1:nvalence
%                 L(i,neighbor1(n)) = cij(n)/area(i);
%         end
%         L(i,i) = L(i,i) * sum_cot / area(i);
% 
%     end
%     L = sparse(L);
    L = build_LB_operator(vNum,vi,f,e);
%       L = build_LB_operator_facebased(v,f);
    deltaX = L * vi / 2;
    for i = 1:vNum
        H(i) = sqrt(deltaX(i,1)^2 + deltaX(i,2)^2 + deltaX(i,3)^2);
    end
    
    % set step length and update the smoothed mesh coordinates
    lamda = 2e-6; % step for the bunny mesh
%     lamda = 1e-2; % step for the bumpy cube mesh
    IS =eye(vNum) - lamda * L;
    vi = IS \ vi;

%     figure;
%     view(2);
%     lighting phong;
%     material shiny
%     axis off;
%     colormap jet(256);
%     light( 'position',[10,10,10]);
% 
%     h = patch('vertices',vi,'faces',f,'FaceVertexCData',H' ,'FaceColor','flat','LineStyle','none');
    h = drawMesh(vi, f, H);
    pause(0.1);
end
pause;
close all;




%% ============================Denoising============================%%
%draw the original mesh
figure;
view(2);
lighting phong;
material shiny
axis off;
colormap jet(256);
light( 'position',[10,10,10]);

h = patch('vertices',v,'faces',f,'FaceVertexCData',1 ,'FaceColor','flat','LineStyle','none');
pause(0.2);

%generate a white Gaussian noise and add it to the original mesh
sigma = 0.20 * mean(v(:));
[length dim] = size(v);
noise = sigma*randn(length,1);
noise = [noise noise noise];

v_noised = v + noise;
figure;
view(2);
lighting phong;
material shiny
axis off;
colormap jet(256);
light( 'position',[10,10,10]);

h = patch('vertices',v_noised,'faces',f,'FaceVertexCData',1 ,'FaceColor','flat','LineStyle','none');
pause(0.2);
difference = abs(v - v_noised);
current_error = sum(difference(:));
error_im_smooth = current_error;
last_error = Inf;

% denoising the mesh until the error stops decreasing
while current_error < last_error
%     L = -eye(vNum);
%     H = zeros(vNum,1);
%     area = zeros(vNum,1);
%     for i = 1:vNum
%         p = v_noised(i,:);   
%         [rowf colf ~] = find(f == i);
%         [row col ~] = find(e == i);
%         nvalence = size(row,1);
%         neighbor1 = [];
%         face = [];
%         edge = [];
%         for j = 1:nvalence
%             face(j,:) = f(rowf(j),:);
%             edge(j,:) = e(row(j),:);
%             neighbor1(j) = e(row(j),-col(j)+3);
%         end
% 
%         v_last = 0;
%         cij = zeros(nvalence,1);
%         sum_cot = 0;
%         for n = 1:nvalence
%             vx = neighbor1(n);
%             v_x = v_noised(vx,:);
%             [vx_f, ~, ~] = find(face == vx);
%             nearest = [face(vx_f(1),:) face(vx_f(2),:)];
%             nearest(find(nearest == vx)) = [];
%             nearest(find(nearest == i)) = [];
%             v_last = v_noised(nearest(1),:);
%             v_next = v_noised(nearest(2),:);
%             p1 = v_x - v_last;
%             q1 = p - v_last;
%             p2 = v_x - v_next;
%             q2 = p - v_next;
%             alpha = acos((q1) * (p1)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2)));
%             beta = acos((q2) * (p2)' / (sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
%             theta = acos((q1) * (q2)' / (sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2)));
%             phi = acos((p1) * (p2)' / (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)));
%             area(i) = area(i) + sqrt(q1(1)^2 + q1(2)^2 + q1(3)^2) * sqrt(q2(1)^2 + q2(2)^2 + q2(3)^2) * sin(theta) / 6 + (sqrt(p1(1)^2 + p1(2)^2 + p1(3)^2) * sqrt(p2(1)^2 + p2(2)^2 + p2(3)^2)) * sin(phi) / 6;
%             cij(n) =  cot(alpha) + cot(beta);
%             sum_cot = sum_cot + cij(n); 
%         end
%         for n = 1:nvalence
%                 L(i,neighbor1(n)) = cij(n)/area(i);
%         end
%         L(i,i) = L(i,i) * sum_cot / area(i);
% 
%     end
%     L = sparse(L);
    L = build_LB_operator(vNum,v_noised,f,e,n);
    deltaX = L * v_noised / 2;
    for i = 1:vNum
        H(i) = sqrt(deltaX(i,1)^2 + deltaX(i,2)^2 + deltaX(i,3)^2);
    end
%     lamda = 1e-6;
    lamda = 2e-6;
    IS =speye(vNum) - lamda * L;
    v_noised = IS \ v_noised;
    figure;
    view(2);
    lighting phong;
    material shiny
    axis off;
    colormap jet(256);
    light( 'position',[10,10,10]);

    h = patch('vertices',v_noised,'faces',f,'FaceVertexCData',1 ,'FaceColor','flat','LineStyle','none');
    pause(0.2);
    
    % update the error
    diffrence = abs(v - v_noised);
    last_error = current_error;
    current_error = sum(diffrence(:));
    error_im_smooth = [error_im_smooth current_error];
end
% count the iteration and draw the error decreasing curve
iteration_im_smooth = size(error_im_smooth,2);
figure;
plot(0:iteration_im_smooth - 2, error_im_smooth(1:size(error_im_smooth,2) - 1));
