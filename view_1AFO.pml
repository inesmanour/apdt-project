reinitialize
load data/1AFO.pdb, prot
as cartoon, prot
bg_color white
set two_sided_lighting, on
set antialias, 2
set cgo_lighting, 1
set transparency_mode, 2
group vis
pseudoatom com, pos=[-1.998,3.699,4.685]
show spheres, com
set sphere_scale, 0.8, com
color red, com
group vis, com
pseudoatom pin_0, pos=[-7.681,30.265,-29.402]
color blue, pin_0
group vis, pin_0
pseudoatom pin_1, pos=[-28.493,30.265,26.873]
color blue, pin_1
group vis, pin_1
pseudoatom pin_2, pos=[-14.483,-27.846,32.054]
color blue, pin_2
group vis, pin_2
pseudoatom pin_3, pos=[6.328,-27.846,-24.221]
color blue, pin_3
group vis, pin_3
pseudoatom pout_0, pos=[10.486,35.244,-22.683]
color red, pout_0
group vis, pout_0
pseudoatom pout_1, pos=[-10.325,35.244,33.592]
color red, pout_1
group vis, pout_1
pseudoatom pout_2, pos=[3.685,-22.867,38.773]
color red, pout_2
group vis, pout_2
pseudoatom pout_3, pos=[24.496,-22.867,-17.502]
color red, pout_3
group vis, pout_3
python
from pymol import cmd
from pymol.cgo import *
com = [-1.998,3.699,4.685]
tip = [20.711,9.923,13.084]
tri_in  = [[-7.681,30.265,-29.402],[-28.493,30.265,26.873],[-14.483,-27.846,32.054],[-7.681,30.265,-29.402],[-14.483,-27.846,32.054],[6.328,-27.846,-24.221]]
tri_out = [[10.486,35.244,-22.683],[-10.325,35.244,33.592],[3.685,-22.867,38.773],[10.486,35.244,-22.683],[3.685,-22.867,38.773],[24.496,-22.867,-17.502]]
quad_in  = [[-7.681,30.265,-29.402],[-28.493,30.265,26.873],[-14.483,-27.846,32.054],[6.328,-27.846,-24.221]]
quad_out = [[10.486,35.244,-22.683],[-10.325,35.244,33.592],[3.685,-22.867,38.773],[24.496,-22.867,-17.502]]
normal_vec = [BEGIN, LINES, COLOR, 0.10,0.45,0.85, VERTEX, -1.998, 3.699, 4.685, VERTEX, 20.711, 9.923, 13.084, END,
 CONE, 20.711, 9.923, 13.084, 25.253, 11.168, 14.764, 1.8, 0.0, 0.10,0.45,0.85, 0.10,0.45,0.85, 1.0, 1.0]
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
