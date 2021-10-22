
  coord v = {(xp - x)/Delta, (yp - y)/Delta}, np;
  {
#line 732

    np.x = - interp (point, v, n.x);
#line 732

    np.y = - interp (point, v, n.y);}
  vertex_buffer_glNormal3d (np.x, np.y, 1.);
  glvertex3d (view, xp, yp, zp);
