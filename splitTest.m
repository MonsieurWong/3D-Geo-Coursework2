
clear all;
close all;
% [v f] = read_ply('dragon_vrip_res3.ply');
[v f] = read_off('bunny.off');
v = v'; 
f = f';

e = compute_edges(f)';
n = compute_normal(v',f')';
vNum = size(v,1);
fNum = size(f,1);
h = drawMesh(v,f,0.7);


L = build_LB_operator(vNum,v,f,e);
deltaX = L * v/ 2;
lamda = 2e-5;
IS =speye(vNum) - lamda * L;
vi = IS \ v;
h = drawMesh(v,f,0.7);


[vo,fo] = triangle_split(vi,f);
h = drawMesh(vo,fo,0.7);