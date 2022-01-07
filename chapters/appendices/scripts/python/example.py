import pygmsh
import numpy
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve 
#geometry definition
def disk(h):
    geom = pygmsh.built_in.Geometry()
    # cirlcle 
    cirlcle_p0 = geom.add_point([0.0, 0.0, 0.0], h)#center
    cirlcle_p1 = geom.add_point([-1.0, 0.0, 0.0], h)
    cirlcle_p2 = geom.add_point([0.0, -1.0, 0.0], h)
    cirlcle_p3 = geom.add_point([1.0, 0.0, 0.0], h)
    cirlcle_p4 = geom.add_point([0.0, 1.0, 0.0], h)
    
    arc_1=geom.add_circle_arc(cirlcle_p1, cirlcle_p0, cirlcle_p2)
    arc_2=geom.add_circle_arc(cirlcle_p2, cirlcle_p0, cirlcle_p3)
    arc_3=geom.add_circle_arc(cirlcle_p3, cirlcle_p0, cirlcle_p4)
    arc_4=geom.add_circle_arc(cirlcle_p4, cirlcle_p0, cirlcle_p1)
    
    circle = geom.add_line_loop([arc_1,arc_2,arc_3,arc_4])
    disk = geom.add_plane_surface(circle)
    #make boundary and area physical
    geom_labels={}
    group_id=geom._TAKEN_PHYSICALGROUP_IDS
    ph_disk = geom.add_physical_surface(disk,label="1")
    geom_labels[group_id[-1]]="1"
    ph_arc_1 = geom.add_physical_line(arc_1,label="1:L1:R0")
    geom_labels[group_id[-1]]="1:L1:R0"
    ph_arc_2 = geom.add_physical_line(arc_2,label="2:L1:R0")
    geom_labels[group_id[-1]]="2:L1:R0"
    ph_arc_3 = geom.add_physical_line(arc_3,label="3:L1:R0")
    geom_labels[group_id[-1]]="3:L1:R0"
    ph_arc_4 = geom.add_physical_line(arc_4,label="4:L1:R0")
    geom_labels[group_id[-1]]="4:L1:R0"
    return geom,geom_labels
    
# boundary conditions
def uD(x1,x2,sdl): 
    f=numpy.zeros(x1.shape)
    I=sdl==0
    f[I]=numpy.square(x1[I])+numpy.square(x2[I])
    I=sdl==1
    f[I]=numpy.square(x1[I])+numpy.square(x2[I])
    return f

def sigma(x1,x2,sdl): 
    f=numpy.zeros(x1.shape)
    I=sdl==0
    f[I]=0
    I=sdl==1
    f[I]=2
    return f

def mu(x1,x2,sdl): 
    f=numpy.zeros(x1.shape)
    I=sdl==0
    f[I]=2*(2+x1[I]+x2[I])*(numpy.square(x1[I])+ \ 
                            numpy.square(x2[I])) 
    I=sdl==1 
    f[I]=2*(2+x1[I]+x2[I]+1)*(numpy.square(x1[I])+ \ 
                            numpy.square(x2[I])) 
    return f
    
#main function
def testPDE():
    np=[] ; nt=[] ; errh2=[] ; err=[] ;
    for h in [0.5, 0.1, 0.05, 0.02]:
        geom,geom_labels=mg.disk(h)
        p,t,e = generateMesh(geom,geom_labels)
        p=numpy.transpose(p)
        t=numpy.transpose(t)
        e=numpy.transpose(e)
        x1=p[0,:]
        x2=p[1,:] 
        def exact(x1,x2,sdl): 
            return numpy.square(x1)+numpy.square(x2)
        def c(x1,x2,sdl): return 2+x1+x2 
        def a(x1,x2,sdl): return x1+x2  
        def f(x1,x2,sdl): 
            return -8-6*(x1+x2)+(2+x1+x2)* \ 
                    (numpy.square(x1)+numpy.square(x2)) 
        def b1(x1,x2,sdl): return x1 
        def b2(x1,x2,sdl): return x2 
        u_exact=numpy.transpose(exact(x1,x2,1))
        bc={}
        bc["bsD"]=[1, 3] 
        bc["uD"]=uD 
        bc["bsR"]=[4, 2] 
        bc["sg"]=sigma 
        bc["mu"]=mu 
        A,F=assemblingAF(p,t,c,a,b1,b2,f) 
        N,S,UD,M=assemblingBC(bc,p,e,4) 
        u=assemblingPDE(bc,p,e,t,c,a,b1,b2,f,1).toarray()
        u_exact=u_exact.reshape(u.shape)
        np.append(p.shape[1]) 
        nt.append(t.shape[1]) 
        h2=h*h
        err.append(numpy.linalg.norm(u_exact-u,numpy.inf)) 
        errh2.append(numpy.linalg.norm(u_exact-u,numpy.inf)/h2)
    print('{:<10s}{:<10s}{:^10s}{:^10s}'.format(
            "np","nt","err","errh2"))
    print()
    for i in range(len(err)):
        print('{:<10d}{:<10d}{:^10.5f}{:^10.5f}'.format(
                np[i],nt[i],err[i],errh2[i]))
