size_t datasize = 12*sizeof (double);
static int metric (const int i, const double t, Event * _ev);
static int metric_expr0 (int * ip, double * tp, Event * _ev);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int cleanup (const int i, const double t, Event * _ev);
static int cleanup_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_1 (const int i, const double t, Event * _ev);
static int defaults_1_expr0 (int * ip, double * tp, Event * _ev);
static int default_display (const int i, const double t, Event * _ev);
static int default_display_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int stability (const int i, const double t, Event * _ev);
static int stability_expr0 (int * ip, double * tp, Event * _ev);
static int vof (const int i, const double t, Event * _ev);
static int vof_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection (const int i, const double t, Event * _ev);
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion (const int i, const double t, Event * _ev);
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev);
static int properties (const int i, const double t, Event * _ev);
static int properties_expr0 (int * ip, double * tp, Event * _ev);
static int advection_term (const int i, const double t, Event * _ev);
static int advection_term_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int projection (const int i, const double t, Event * _ev);
static int projection_expr0 (int * ip, double * tp, Event * _ev);
static int end_timestep (const int i, const double t, Event * _ev);
static int end_timestep_expr0 (int * ip, double * tp, Event * _ev);
static int adapt (const int i, const double t, Event * _ev);
static int adapt_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int logfile (const int i, const double t, Event * _ev);
static int logfile_expr0 (int * ip, double * tp, Event * _ev);
static int logfile_expr1 (int * ip, double * tp, Event * _ev);
static int profile (const int i, const double t, Event * _ev);
static int profile_expr0 (int * ip, double * tp, Event * _ev);
static int movies (const int i, const double t, Event * _ev);
static int movies_expr0 (int * ip, double * tp, Event * _ev);
static int movies_expr1 (int * ip, double * tp, Event * _ev);
static void _set_boundary0 (void);
static void _set_boundary1 (void);
static void _set_boundary2 (void);
static void _set_boundary3 (void);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, metric, {metric_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/embed.h", 891, "metric"});
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/embed.h", 929, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/run.h", 42, "defaults"});
  event_register ((Event){ 0, 1, defaults_1, {defaults_1_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 126, "defaults"});
  event_register ((Event){ 0, 1, default_display, {default_display_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 179, "default_display"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 188, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((int *)0), ((double *)0),
    "poiseuille45.c", 43, "init"});
  event_register ((Event){ 0, 2, logfile, {logfile_expr0, logfile_expr1}, ((int *)0), ((double *)0),
    "poiseuille45.c", 68, "logfile"});
  event_register ((Event){ 0, 1, profile, {profile_expr0}, ((int *)0), ((double *)0),
    "poiseuille45.c", 75, "profile"});
  event_register ((Event){ 0, 2, movies, {movies_expr0, movies_expr1}, ((int *)0), ((double *)0),
    "poiseuille45.c", 99, "movies"});
  event_register ((Event){ 0, 1, cleanup, {cleanup_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/run.h", 50, "cleanup"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 216, "set_dtmax"});
  event_register ((Event){ 0, 1, stability, {stability_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 218, "stability"});
  event_register ((Event){ 0, 1, vof, {vof_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 228, "vof"});
  event_register ((Event){ 0, 1, tracer_advection, {tracer_advection_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 229, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_diffusion, {tracer_diffusion_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 230, "tracer_diffusion"});
  event_register ((Event){ 0, 1, properties, {properties_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 237, "properties"});
  event_register ((Event){ 0, 1, advection_term, {advection_term_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 314, "advection_term"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 344, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 380, "acceleration"});
  event_register ((Event){ 0, 1, projection, {projection_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 423, "projection"});
  event_register ((Event){ 0, 1, end_timestep, {end_timestep_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 438, "end_timestep"});
  event_register ((Event){ 0, 1, adapt, {adapt_expr0}, ((int *)0), ((double *)0),
    "/opt/basilisk/src/navier-stokes/centered.h", 448, "adapt"});
  _attribute = (_Attributes *) pcalloc (datasize/sizeof(double), sizeof (_Attributes), __func__, __FILE__, __LINE__);
  all = (scalar *) pmalloc (sizeof (scalar)*13,__func__, __FILE__, __LINE__);
  for (int i = 0; i < 12; i++)
    all[i].i = i;
  all[12].i = -1;
  set_fpe();
  quadtree_methods();
  init_scalar ((scalar){11}, "un");
  init_face_vector ((vector){{9},{10}}, "uf");
  init_scalar ((scalar){8}, "pf");
  init_vector ((vector){{6},{7}}, "g");
  init_vector ((vector){{4},{5}}, "u");
  init_scalar ((scalar){3}, "p");
  embed = new_bid();
  init_face_vector ((vector){{1},{2}}, "fs");
  init_scalar ((scalar){0}, "cs");
  init_const_scalar ((scalar){_NVARMAX+5}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+4}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+2},{_NVARMAX+3}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1}}, "zerof", (double []) {0.,0.,0.});
  _set_boundary0();
  _set_boundary1();
  _set_boundary2();
  _set_boundary3();
}
