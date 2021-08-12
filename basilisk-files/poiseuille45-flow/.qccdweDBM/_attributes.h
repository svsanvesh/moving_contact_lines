typedef struct {

#line 945 "/opt/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  bool face, nodump, freed;
  int block;

#line 17 "/opt/basilisk/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

#line 8 "/opt/basilisk/src/grid/tree-common.h"

  void (* refine) (Point, scalar);

#line 94 "/opt/basilisk/src/grid/tree-common.h"

  void (* coarsen) (Point, scalar);

#line 81 "/opt/basilisk/src/fractions.h"

  vector n;

#line 206 "/opt/basilisk/src/embed-tree.h"

  void (* embed_gradient) (Point, scalar, coord *);

#line 175 "/opt/basilisk/src/embed.h"

  bool third;

} _Attributes;
_Attributes * _attribute;
