function h = drawMesh(v,f,color)

    figure;
    view(2);
    lighting phong;
    material shiny
    axis off;
    colormap jet(256);
    light( 'position',[10,10,10]);
    if size(color,2) > 1
        color = color';
    end
    h = patch('vertices',v,'faces',f,'FaceVertexCData',color ,'FaceColor','flat','LineStyle','-');
end