reinitialize
load data/1XQF.pdb, prot
as cartoon, prot
bg_color white
set two_sided_lighting, on
set antialias, 2
set cgo_lighting, 1
set transparency_mode, 2
group vis
pseudoatom com, pos=[43.632,16.764,14.631]
show spheres, com
set sphere_scale, 0.8, com
color red, com
group vis, com
pseudoatom pin_0, pos=[73.210,-13.236,9.400]
color blue, pin_0
group vis, pin_0
pseudoatom pin_1, pos=[73.210,46.764,9.400]
color blue, pin_1
group vis, pin_1
pseudoatom pin_2, pos=[13.679,46.764,16.885]
color blue, pin_2
group vis, pin_2
pseudoatom pin_3, pos=[13.679,-13.236,16.885]
color blue, pin_3
group vis, pin_3
pseudoatom pout_0, pos=[73.585,-13.236,12.377]
color red, pout_0
group vis, pout_0
pseudoatom pout_1, pos=[73.585,46.764,12.377]
color red, pout_1
group vis, pout_1
pseudoatom pout_2, pos=[14.053,46.764,19.862]
color red, pout_2
group vis, pout_2
pseudoatom pout_3, pos=[14.053,-13.236,19.862]
color red, pout_3
group vis, pout_3
python
from pymol import cmd
from pymol.cgo import *
com = [43.632,16.764,14.631]
tip = [46.751,16.764,39.436]
tri_in  = [[73.210,-13.236,9.400],[73.210,46.764,9.400],[13.679,46.764,16.885],[73.210,-13.236,9.400],[13.679,46.764,16.885],[13.679,-13.236,16.885]]
tri_out = [[73.585,-13.236,12.377],[73.585,46.764,12.377],[14.053,46.764,19.862],[73.585,-13.236,12.377],[14.053,46.764,19.862],[14.053,-13.236,19.862]]
quad_in  = [[73.210,-13.236,9.400],[73.210,46.764,9.400],[13.679,46.764,16.885],[13.679,-13.236,16.885]]
quad_out = [[73.585,-13.236,12.377],[73.585,46.764,12.377],[14.053,46.764,19.862],[14.053,-13.236,19.862]]
normal_vec = [BEGIN, LINES, COLOR, 0.10,0.45,0.85, VERTEX, 43.632, 16.764, 14.631, VERTEX, 46.751, 16.764, 39.436, END,
 CONE, 46.751, 16.764, 39.436, 47.375, 16.764, 44.397, 1.8, 0.0, 0.10,0.45,0.85, 0.10,0.45,0.85, 1.0, 1.0]
cmd.load_cgo(normal_vec, 'normal_vec')
def _emit_tri(name, T, rgb):
    r,g,b = rgb
    obj=[BEGIN, TRIANGLES, COLOR, r,g,b, ALPHA, 0.28]
    for x,y,z in T: obj += [VERTEX, float(x),float(y),float(z)]
    obj += [END]
    cmd.load_cgo(obj, name)
def _emit_outline(name, Q, rgb):
    r,g,b = rgb
    order=[0,1,2,3,0]
    obj=[BEGIN, LINES, COLOR, r,g,b]
    for i in range(4):
        a,b = order[i], order[i+1]
        x1,y1,z1 = Q[a]; x2,y2,z2 = Q[b]
        obj += [VERTEX, float(x1),float(y1),float(z1), VERTEX, float(x2),float(y2),float(z2)]
    obj += [END]
    cmd.load_cgo(obj, name)
_emit_tri('plane_in',  tri_in,  (0.20,0.60,1.00))
_emit_tri('plane_out', tri_out, (1.00,0.30,0.30))
_emit_outline('plane_in_outline',  quad_in,  (0.20,0.60,1.00))
_emit_outline('plane_out_outline', quad_out, (1.00,0.30,0.30))
python end
enable vis
enable plane_in
enable plane_out
enable plane_in_outline
enable plane_out_outline
enable normal_vec
zoom vis, buffer=10
