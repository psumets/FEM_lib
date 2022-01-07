import pygmsh
import numpy
def generateMesh(geom,geom_labels):
    # generate mesh
    points, cells, point_data, cell_data, field_data= \
                                pygmsh.generate_mesh(geom) 
    #extract p,t,e matrices
    p=points[:,0:2]
    t=numpy.zeros((cells["triangle"].shape[0],4))
    e=numpy.zeros((cells["line"].shape[0],7))
    t[:,0:3]=cells["triangle"]
    e[:,0:2]=cells["line"]
    for label in field_data:
        elemnt_index=field_data[label][0]
        group_index=field_data[label][1]
        names=list(cell_data.keys())
        if names[group_index-1]=="triangle":# domain labels
            t[cell_data["triangle"] \
                ["gmsh:physical"]==elemnt_index,3]= \ 
                int(geom_labels[elemnt_index])
        elif names[group_index-1]=="line":# boundary labels
            label=None
            left_label=None
            right_label=None
            for item in geom_labels[elemnt_index].split(":"):
                if item[0]=="L": left_label=int(item[1:])
                elif item[0]=="R": right_label=int(item[1:])
                else: label=int(item)
            if label==None or \ 
                left_label==None or \ 
                right_label==None:
                print('error, wrong geometry boundary labelling') 
            index=cell_data["line"]["gmsh:physical"]==elemnt_index   
            e[index,4]=label
            e[index,5]=left_label
            e[index,6]=right_label
    return p,t.astype(int),e.astype(int)