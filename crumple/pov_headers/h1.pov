// W=1024 H=768
// Chris's example
#version 3.6;

global_settings {
	max_trace_level 64
}

camera {
	location <5,-27,13>
	sky z
	right -0.072*x*image_width/image_height
	up 0.072*z
	look_at <0,0,0>
}

background{rgb 0}

light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}

#declare t_mesh=texture{
	pigment{rgbft <1,0.8,0.6,0,0.4>}
	finish{specular 0.4 ambient 0.4}
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
