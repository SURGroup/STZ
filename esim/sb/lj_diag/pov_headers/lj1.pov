#declare a=0.94;

camera {
  orthographic
  location <0,0,-200>
  look_at  <0,0,0>
  right 11*image_width/image_height*x
  up 11*y
}

background{color rgb 1}

light_source {<-25,35,-25> color rgb <1,0.3,1>}
light_source { <10,30,10> color rgb <0.5,1,0.5>}

#declare t1=texture{
pigment{rgbft<1,.6,.3,0,.3>}     //was 0.6
finish{ambient 0.6 diffuse 1 specular .2 phong .2 phong_size 10 reflection 0.1}
}
#declare t2=texture{
pigment{rgbft<1,.6,.3,0,0.2>}    //was 0.5
finish{ambient 0.6 diffuse 1 specular .2 phong .2 phong_size 10 reflection 0.1}
}
#declare t3=texture{
pigment{rgb<0.8,0,0>}
finish{diffuse 1.0 ambient 0.5 specular 0.3}
}
#declare t4=texture{
pigment{rgb<0,0.8,0>}
finish{diffuse 1.0 ambient 0.5 specular 0.2}
}

union{
#include "sph.pov"
texture{t1}
}
