function D = duplication_matrix(n)
    m=1/2*n*(n+1);
    nsq=n^2;
    Lind=tril(true(n));
    Lind=find(Lind(:));
    Lind=Lind(:);
    Uind=rem(Lind-1,n)*n+ceil(Lind/n);
    i=(1:m)';
    a=(i-1)*nsq+Lind;
    b=(i-1)*nsq+Uind;
    c=union(a,b);
    [I,J]=ind2sub([nsq,m],c);
    D=sparse(I,J,1,nsq,m);
end