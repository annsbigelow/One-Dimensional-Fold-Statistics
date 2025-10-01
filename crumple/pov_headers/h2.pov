// W=1600 H=1600
// Chris's example
#version 3.6;

global_settings {
	max_trace_level 64
}

camera {
    orthographic
	location <0,0,500>
	sky y
	right -135*x*image_width/image_height
	up 135*y
	look_at <0,0,0>
}

background{rgb 0}

light_source{<-200,240,500> color rgb <0.77,0.75,0.75>}
light_source{<400,-100,250> color rgb <0.38,0.40,0.40>}

#declare t_mesh=texture{
	pigment{gradient z
		pigment_map {
            [0 rgb <0.15,0.3,1>]
            [0.3 rgb <0.22,0.39,0.82>]
            [0.4 rgb <0.5,0.82,0.71>]
            [0.5 rgb <0.92,0.88,0.43>]
            [0.6 rgb <0.90,0.59,0.51>]
            [0.7 rgb <0.90,0.20,0.86>]
            [1 rgb <1,0.08,0.98>]
		}
		scale 24
		translate <0,0,-12>
	}
	finish{specular 0.4 ambient 0.3}
}

#include "msh.pov"

CYL:#declare s=0.0025;
CYL:union {
CYL:#include "cyl.pov"
CYL:	texture {
CYL:		pigment{rgb <0.2,0.3,0.8>}
CYL:		finish{specular 0.4 ambient 0.4}
CYL:	}
CYL:}
