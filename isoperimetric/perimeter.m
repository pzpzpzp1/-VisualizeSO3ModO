function p = perimeter(x)

x = x(:);
x = reshape(x,2,[])';
dx = (x - circshift(x,1,1));
p = sum(sqrt(sum(dx.*dx,2)));

end


