function value = gridAle2(grid1,grid2,grid3,x,w)

xi_iter=1;
if x(1,3)==grid3(2)
    xi_iter=2;
end

value=interpAle2(grid1,grid2,w(:,:,xi_iter),x(:,1),x(:,2))';

end