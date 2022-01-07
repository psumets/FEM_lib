def uD(x1,x2,sdl): 
    f=numpy.zeros(x1.shape)
    # segments in local numeration:
    # 1-st segment 
    I=sdl==0
    f[I]=numpy.square(x1[I])+numpy.square(x2[I])
    # 2-nd segment 
    I=sdl==1
    f[I]=numpy.square(x1[I])+numpy.square(x2[I])
    return f