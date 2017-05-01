function edge_length = edgeLength(v,e)

    e1 = e(:,1);
    e2 = e(:,2);
    fun = @(e1,e2) norm(v(e1,:) - v(e2,:));
    edge_length = arrayfun(fun,e1,e2);
    
end