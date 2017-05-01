%%============================edge collapse=============================%%% function [v_o f_o] = remesh()


%%============================Initialization============================%%
clear all;
close all;
% [v f] = read_ply('dragon_vrip_res3.ply');
[v f] = read_off('bumpy.off');
v = v'; 
f = f';

e = compute_edges(f)';
n = compute_normal(v',f')';
vNum = size(v,1);
fNum = size(f,1);


%%============================edge collapse============================%%


edge_length1 = edgeLength(v,e);
for i = 1:10
    vNum = size(v,1);
    L = build_LB_operator(vNum,v,f,e);
2; % step for the bumpy cube mesheLength(vi,e);
    % ratio = abs(edge_length2 - edge_length1)./edge_length1;
    a = edge_length2;
    a(a>0.2*mean(edge_length1)) = 0;
    [collapsedEdge, ~] = find(a > 0);

    vi(e(collapsedEdge,1),:) = (vi(e(collapsedEdge,1),:) + vi(e(collapsedEdge,2),:))/2;

    for i = 1:size(collapsedEdge,1)
        [collapsedFace1 ~] = find(f == e(collapsedEdge(i),1));
        [collapsedFace2 ~] = find(f == e(collapsedEdge(i),2));
        collapsedFacetem = intersect(collapsedFace1,collapsedFace2);
        collapsedFace = [collapsedFace; collapsedFacetem];
    end
    f(collapedFace,:) = [];
%     vi(e(collapsedEdge,2),:) = [];
    e = compute_edges(f)';
    
%     write_off('file.off',vi,f);
    clear v;
    v = vi;
end

%%=============================face split===============================%%
angle_thresh = 110*pi/180;
angles = get_angles(v,f);

for i = 1:size(f,1)
    bool triangle_split = false;
    split_angle = false(3,1);
    for j = 1:3
        split_angle(j) = (angles(i,j)>=angle_thresh) ;
        triangle_split = triangle_split || split_angle(j);
    end
    
    %split triangles
    if triangle_split
%         delete_idx = [delete_idx; i]; 
        col = find(split_angle>0);
        idx_v1 = col;
        idx_v2 = rem(col,3)+1;
        idx_v3 = rem(idx_v2,3)+1;
        v1 = v(f(i,idx_v1),:);
        v2 = v(f(i,idx_v2),:);
        v3 = v(f(i,idx_v3),:);
        p = (v3 - v2)/norm(v3 - v2);
        q = (v1 - v2)/norm(v1 - v2);
        tao = cross(p,cross(p,q));
        d = v1 - tao * dot(v1,(v1 - v2));
        idx_d = size(v) + 1;
        f(i,:) = [];
        f = [f ; idx_v1 idx_v2 idx_d ; idx_v1 idx_v3 idx_d];
    end
   
    
end


    


