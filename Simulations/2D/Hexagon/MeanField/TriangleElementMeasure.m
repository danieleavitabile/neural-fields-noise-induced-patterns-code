function S = TriangleElementMeasure(x,y,z)

    S = ...
    sqrt(...
       ( x(2)*y(1) - x(3)*y(1) - x(1)*y(2) + x(3)*y(2) + x(1)*y(3) - x(2)*y(3) )^2 ...
     + ( x(2)*z(1) - x(3)*z(1) - x(1)*z(2) + x(3)*z(2) + x(1)*z(3) - x(2)*z(3) )^2 ...
     + ( y(2)*z(1) - y(3)*z(1) - y(1)*z(2) + y(3)*z(2) + y(1)*z(3) - y(2)*z(3) )^2 ...
    )/2;

end