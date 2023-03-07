function d = geographic_distance(i,j,XY)
    p1 = XY(i,:);
    p2 = XY(j,:);
    d = norm(p1-p2,2);
end
