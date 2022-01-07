def assemblingAF(c,b,a,f,xl,x,d,w):
    E,D=interpolatingMat(x,d)
    ET=numpy.transpose(E)
    DT=numpy.transpose(D)
    m=numpy.size(x)-1
    n=numpy.size(xl)-1
    n_total=n*m+1;
    A=sparse.coo_matrix((n_total,n_total))
    F=sparse.coo_matrix((n_total,1))
    i_loc=numpy.arange(m + 1)
    i_loc_sp_mat = numpy.repeat(i_loc, m + 1)
    j_loc_sp_mat = numpy.tile(i_loc, m + 1)
    for l in range(n):
        h=xl[l+1]-xl[l] # size of element l
        xd=xl[l]+h*d # quadrature nodes on element l
        c_loc=numpy.diag((w/h)*c(xd))
        b_loc=numpy.diag(w*b(xd))
        a_loc=numpy.diag((h*w)*a(xd))
        Al=(DT @ c_loc+ET @ b_loc) @ D+ET @ a_loc @ E
        row=l*m+i_loc_sp_mat
        col=l*m+j_loc_sp_mat
        A=A+sparse.coo_matrix((Al.flatten(),(row, col)), \ 
            shape=(n_total,n_total))
        Fl=ET @ (h*w*f(xd)).reshape(-1,1)
        F=F+sparse.coo_matrix((Fl.flatten(), \ 
            (l*m+i_loc, numpy.zeros(m+1))),shape=(n_total,1))
    return A,F