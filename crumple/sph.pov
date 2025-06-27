#version 3.6;

// Right-handed coordinate system in which the z-axis points upwards
camera {
	location <30,-50,25>
	sky z
	right -0.05*x*image_width/image_height
	up 0.05*z
	look_at <0,0,0>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}

// Radius of the cylinders
#declare s=0.002;

// Particles
union{
#include "sph_1000_c.pov"
	pigment{rgb <0.1,0.1,0.5>} finish{specular 0.4 ambient 0.42}
}

// Voronoi cells
#declare t_mesh=texture{
	pigment{rgb <0.2,0.7,0.85>} finish{specular 0.8 ambient 0.42}
}
#include "sph_1000_m.pov"
