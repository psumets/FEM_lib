import pygmsh
import numpy
def geometry(h):
    geom = pygmsh.built_in.Geometry() # initialize geometry object
    #rectangular #1
    rect1_p1=geom.add_point([0.0, -1.0, 0.0], h)
    rect1_p2=geom.add_point([0.0, 0.0, 0.0], h)
    rect1_p3=geom.add_point([2.0, 0.0, 0.0], h)
    rect1_p4=geom.add_point([2.0, -1.0, 0.0], h)
    rect1_l1=geom.add_line(rect1_p1, rect1_p2)
    rect1_l2=geom.add_line(rect1_p2, rect1_p3)
    rect1_l3=geom.add_line(rect1_p3, rect1_p4)
    rect1_l4=geom.add_line(rect1_p4, rect1_p1)
    # circular hole 
    cirlcle_p0=geom.add_point([1.0, -0.5, 0.0], h) #center
    cirlcle_p1=geom.add_point([0.75, -0.5, 0.0], h)
    cirlcle_p2=geom.add_point([1.0, -0.75, 0.0], h)
    cirlcle_p3=geom.add_point([1.25, -0.5, 0.0], h)
    cirlcle_p4=geom.add_point([1.0, -0.25, 0.0], h)
    arc_1=geom.add_circle_arc(cirlcle_p1, cirlcle_p0, cirlcle_p2)
    arc_2=geom.add_circle_arc(cirlcle_p2, cirlcle_p0, cirlcle_p3)
    arc_3=geom.add_circle_arc(cirlcle_p3, cirlcle_p0, cirlcle_p4)
    arc_4=geom.add_circle_arc(cirlcle_p4, cirlcle_p0, cirlcle_p1)
    rect1_line=geom.add_line_loop([rect1_l1,rect1_l2, \
                                rect1_l3,rect1_l4, \
                                arc_1,arc_2,arc_3,arc_4])
    #rectangle #2
    rect2_p1=geom.add_point([0.75, 0.25, 0.0], h)
    rect2_p2=geom.add_point([0.75, 0.75, 0.0], h)
    rect2_p3=geom.add_point([1.25, 0.75, 0.0], h)
    rect2_p4=geom.add_point([1.25, 0.25, 0.0], h)
    rect2_l1=geom.add_line(rect2_p1, rect2_p2)
    rect2_l2=geom.add_line(rect2_p2, rect2_p3)
    rect2_l3=geom.add_line(rect2_p3, rect2_p4)
    rect2_l4=geom.add_line(rect2_p4, rect2_p1)
    rect2_line=geom.add_line_loop([rect2_l1,rect2_l2, \
                                rect2_l3,rect2_l4])
    #rectangle #3
    rect3_p1=geom.add_point([0.0, 0.0, 0.0], h)
    rect3_p2=geom.add_point([0.0, 1.0, 0.0], h)
    rect3_p3=geom.add_point([2.0, 1.0, 0.0], h)
    rect3_p4=geom.add_point([2.0, 0.0, 0.0], h)
    rect3_l1=geom.add_line(rect3_p1, rect3_p2)
    rect3_l2=geom.add_line(rect3_p2, rect3_p3)
    rect3_l3=geom.add_line(rect3_p3, rect3_p4)
    rect3_l4=geom.add_line(rect3_p4, rect3_p1)
    rect3_line=geom.add_line_loop([rect3_l1,rect3_l2, \
                            rect3_l3,rect3_l4,rect2_l1, \
                            rect2_l2,rect2_l3,rect2_l4])
    #create surfaces
    rect1_surf=geom.add_plane_surface(rect1_line)
    rect3_surf=geom.add_plane_surface(rect3_line)
    rect2_surf=geom.add_plane_surface(rect2_line)
    # define labels for geometry elements:
    geom_labels={}
    group_id=geom._TAKEN_PHYSICALGROUP_IDS
    # domain labels
    ph_rect1_surf=geom.add_physical_surface(rect1_surf,label="1")
    geom_labels[group_id[-1]]="1"
    ph_rect2_surf=geom.add_physical_surface(rect2_surf,label="2")
    geom_labels[group_id[-1]]="2"
    ph_rect3_surf=geom.add_physical_surface(rect3_surf,label="3")
    geom_labels[group_id[-1]]="3"
    # boundary labels
    ph_rect1_l1=geom.add_physical_line(rect1_l1,label="7:L0:R1")
    geom_labels[group_id[-1]]="7:L0:R1"
    ph_rect1_l2=geom.add_physical_line(rect1_l2,label="11:L3:R1")
    geom_labels[group_id[-1]]="11:L3:R1"
    ph_rect1_l3=geom.add_physical_line(rect1_l3,label="4:L0:R1")
    geom_labels[group_id[-1]]="4:L0:R1"
    ph_rect1_l4=geom.add_physical_line(rect1_l4,label="1:L0:R1")
    geom_labels[group_id[-1]]="1:L0:R1"
    ph_rect2_l1=geom.add_physical_line(rect2_l1,label="10:L3:R2")
    geom_labels[group_id[-1]]="10:L3:R2"
    ph_rect2_l2=geom.add_physical_line(rect2_l2,label="6:L3:R2")
    geom_labels[group_id[-1]]="6:L3:R2"
    ph_rect2_l3=geom.add_physical_line(rect2_l3,label="2:L3:R2")
    geom_labels[group_id[-1]]="2:L3:R2"
    ph_rect2_l4=geom.add_physical_line(rect2_l4,label="3:L3:R2")
    geom_labels[group_id[-1]]="3:L3:R2"
    ph_rect3_l1=geom.add_physical_line(rect3_l1,label="8:L0:R3")
    geom_labels[group_id[-1]]="8:L0:R3"
    ph_rect3_l2=geom.add_physical_line(rect3_l2,label="9:L0:R3")
    geom_labels[group_id[-1]]="9:L0:R3"
    ph_rect3_l3=geom.add_physical_line(rect3_l3,label="5:L0:R3")
    geom_labels[group_id[-1]]="5:L0:R3"
    ph_arc_1=geom.add_physical_line(arc_1,label="12:L0:R1")
    geom_labels[group_id[-1]]="12:L0:R1"
    ph_arc_2=geom.add_physical_line(arc_2,label="13:L0:R1")
    geom_labels[group_id[-1]]="13:L0:R1"
    ph_arc_3=geom.add_physical_line(arc_3,label="14:L0:R1")
    geom_labels[group_id[-1]]="14:L0:R1"
    ph_arc_4=geom.add_physical_line(arc_4,label="15:L0:R1")
    geom_labels[group_id[-1]]="15:L0:R1"
    return geom,geom_labels