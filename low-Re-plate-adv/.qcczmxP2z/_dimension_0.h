 {
     memcpy (b, &val(v.x,0,0,0), sizeof(double)*_attribute[v.x.i].block);
     b += _attribute[v.x.i].block;
     if (allocated(1,0,0))
       memcpy (b, &val(v.x,1,0,0), sizeof(double)*_attribute[v.x.i].block);
     else
       *b = nodata;
     b += _attribute[v.x.i].block;
   }