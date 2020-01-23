#declare a=0.94;

camera {
  orthographic
  location <0,0,-200>
  look_at  <0,0,0>
  right 10*x
  up 11*y
}

background{color rgb 1}

light_source {<-28,35,-45> color rgb <0.76,0.73,0.73>}
light_source { <30,-10,-25> color rgb <0.33,0.36,0.36>}

#declare f1=finish{ambient 0.4 diffuse 1 specular .2 phong .15 reflection 0.08}

#declare t0=texture{pigment{rgb<0.35,0.6,0.8>} finish{f1}}
#declare t1=texture{pigment{rgb<0.8,0.6,0.35>} finish{f1}}
#declare t2=texture{pigment{rgb<0.8,0.6,0.35>} finish{f1}}
#declare t3=texture{pigment{rgb<0.8,0.2,0.2>} finish{f1}}
#declare t4=texture{pigment{rgb<0.4,0.15,0.5>} finish{f1}}
#declare t5=texture{pigment{rgb<0.65,0.75,0.65>} finish{f1}}

#include "sph.pov"
