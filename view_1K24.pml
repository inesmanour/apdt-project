reinitialize
load data/1K24.pdb, prot
as cartoon, prot
bg_color white
set two_sided_lighting, on
set antialias, 2
set cgo_lighting, 1
set transparency_mode, 2
group vis
pseudoatom com, pos=[30.403,-2.402,9.879]
show spheres, com
set sphere_scale, 0.8, com
color red, com
group vis, com
pseudoatom pin_0, pos=[65.630,-31.601,20.200]
color blue, pin_0
group vis, pin_0
pseudoatom pin_1, pos=[65.630,18.235,-13.212]
color blue, pin_1
group vis, pin_1
pseudoatom pin_2, pos=[9.759,6.055,-31.380]
color blue, pin_2
group vis, pin_2
pseudoatom pin_3, pos=[9.759,-43.781,2.032]
color blue, pin_3
group vis, pin_3
pseudoatom pout_0, pos=[51.047,-10.859,51.138]
color red, pout_0
group vis, pout_0
pseudoatom pout_1, pos=[51.047,38.977,17.726]
color red, pout_1
group vis, pout_1
pseudoatom pout_2, pos=[-4.824,26.796,-0.443]
color red, pout_2
group vis, pout_2
pseudoatom pout_3, pos=[-4.824,-23.040,32.970]
color red, pout_3
group vis, pout_3
python
from pymol import cmd
from pymol.cgo import *
com = [30.403,-2.402,9.879]
tip = [21.289,10.562,29.215]
tri_in  = [[65.630,-31.601,20.200],[65.630,18.235,-13.212],[9.759,6.055,-31.380],[65.630,-31.601,20.200],[9.759,6.055,-31.380],[9.759,-43.781,2.032]]
tri_out = [[51.047,-10.859,51.138],[51.047,38.977,17.726],[-4.824,26.796,-0.443],[51.047,-10.859,51.138],[-4.824,26.796,-0.443],[-4.824,-23.040,32.970]]
quad_in  = [[65.630,-31.601,20.200],[65.630,18.235,-13.212],[9.759,6.055,-31.380],[9.759,-43.781,2.032]]
quad_out = [[51.047,-10.859,51.138],[51.047,38.977,17.726],[-4.824,26.796,-0.443],[-4.824,-23.040,32.970]]
normal_vec = [BEGIN, LINES, COLOR, 0.10,0.45,0.85, VERTEX, 30.403, -2.402, 9.879, VERTEX, 21.289, 10.562, 29.215, END,
 CONE, 21.289, 10.562, 29.215, 19.466, 13.154, 33.082, 1.8, 0.0, 0.10,0.45,0.85, 0.10,0.45,0.85, 1.0, 1.0]
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
