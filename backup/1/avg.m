function B=avg(A,idim)
if(idim==1)
    B= [ A(2:end,:) + A(1:end-1,:) ]/2;
    elseif(idim==2)
    B=[ A(:,2:end) + A(:,1:end-1) ]/2;
else
    error('avg(A,idim):idim must be 1 or 2');
end
