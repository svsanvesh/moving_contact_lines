{

#line 1575 "/home/fpl/softwares/basilisk/src/draw.h"

      if (val(f,0,0,0) != nodata) {
 glPushMatrix();
 char s[80];
 sprintf (s, "%g", val(f,0,0,0));
 float scale = 0.8*Delta_x/(strlen(s)*width);
 glTranslatef (x - 0.4*Delta_x, y - scale*height/3., 0.);
 glScalef (scale, scale, 1.);
 gl_StrokeString (s);
 glPopMatrix();
      }