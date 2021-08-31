#line 91 "/opt/basilisk/src/navier-stokes/centered.h"
static double _boundary0 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 90 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); } return 0.; } static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 90 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 91
return  neumann_homogeneous(); } return 0.; }
#line 92 "/opt/basilisk/src/navier-stokes/centered.h"
static double _boundary1 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 91 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); } return 0.; } static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 91 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 92
return  neumann_homogeneous(); } return 0.; }
#line 101 "/opt/basilisk/src/navier-stokes/centered.h"
static double _boundary2 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 100 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0)/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),0))); } return 0.; } static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 100 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 101
return  neumann_homogeneous(); } return 0.; }
#line 102 "/opt/basilisk/src/navier-stokes/centered.h"
static double _boundary3 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 101 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); } return 0.; } static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 101 "/opt/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 102
return  neumann_homogeneous(); } return 0.; }
#line 44 "test.c"
static double _boundary4 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 43 "test.c"
return  dirichlet(1.0); return 0.; } static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 43 "test.c"
return  dirichlet_homogeneous(); return 0.; }
#line 45 "test.c"
static double _boundary5 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 44 "test.c"
return  dirichlet(0.0); return 0.; } static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 44 "test.c"
return  dirichlet_homogeneous(); return 0.; }
#line 47 "test.c"
static double _boundary6 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 46 "test.c"
return  dirichlet(0.0); return 0.; } static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 46 "test.c"
return  dirichlet_homogeneous(); return 0.; }
#line 48 "test.c"
static double _boundary7 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 47 "test.c"
return  dirichlet(0.0); return 0.; } static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 47 "test.c"
return  dirichlet_homogeneous(); return 0.; }
#line 130 "test.c"
static double _boundary8 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 129 "test.c"
return  dirichlet(1.); return 0.; } static double _boundary8_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 129 "test.c"
return  dirichlet_homogeneous(); return 0.; }
#line 131 "test.c"
static double _boundary9 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 130 "test.c"
return  dirichlet(0.); return 0.; } static double _boundary9_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); POINT_VARIABLES; 
#line 130 "test.c"
return  dirichlet_homogeneous(); return 0.; }
