function value = gridAle(grid1,grid2,x,w)

xi_iter=1;
if x(1,2)==grid2(2)
    xi_iter=2;
end

value=interpAle(grid1,w(:,xi_iter),x(:,1))';

end