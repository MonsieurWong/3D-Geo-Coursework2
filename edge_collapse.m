function [v, f] = edge_collapse(v,f)
    e = compute_edges(f)';
    edge_length1 = edgeLength(v,e);

 % step for the bumpy cube mesheLength(vi,e);
    % ratio = abs(edge_length2 - edge_length1)./edge_length1;
    a = edge_length1;
%     a(a>0.2*mean(edge_length1)) = 0;
    a(a>0.15) = 0;
    [collapsedEdge, ~] = find(a > 0);

    v(e(collapsedEdge,1),:) = (v(e(collapsedEdge,1),:) + v(e(collapsedEdge,2),:))/2;
    collapsedFace = [];
    for i = 1:size(collapsedEdge,1)
        [collapsedFace1 ~] = find(f == e(collapsedEdge(i),1));
        [collapsedFace2 ~] = find(f == e(collapsedEdge(i),2));
        collapsedFacetem = intersect(collapsedFace1,collapsedFace2);
        collapsedFace = [collapsedFace; collapsedFacetem];
    end
    f(collapedFace,:) = [];
%     vi(e(collapsedEdge,2),:) = [];
    e = compute_edges(f)';
    
    clear v;
    v = v;
end