/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__it2
#define _nrn_initial _nrn_initial__it2
#define nrn_cur _nrn_cur__it2
#define _nrn_current _nrn_current__it2
#define nrn_jacob _nrn_jacob__it2
#define nrn_state _nrn_state__it2
#define _net_receive _net_receive__it2 
#define castate castate__it2 
#define evaluate_fct evaluate_fct__it2 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gcabar _p[0]
#define g _p[1]
#define m _p[2]
#define h _p[3]
#define Cai _p[4]
#define Cao _p[5]
#define Dm _p[6]
#define Dh _p[7]
#define iCa _p[8]
#define carev _p[9]
#define v _p[10]
#define _g _p[11]
#define _ion_Cai	*_ppvar[0]._pval
#define _ion_Cao	*_ppvar[1]._pval
#define _ion_iCa	*_ppvar[2]._pval
#define _ion_diCadv	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_efun(void);
 static void _hoc_evaluate_fct(void);
 static void _hoc_ghk(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_it2", _hoc_setdata,
 "efun_it2", _hoc_efun,
 "evaluate_fct_it2", _hoc_evaluate_fct,
 "ghk_it2", _hoc_ghk,
 0, 0
};
#define efun efun_it2
#define ghk ghk_it2
 extern double efun( _threadargsprotocomma_ double );
 extern double ghk( _threadargsprotocomma_ double , double , double );
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[6];
#define _gth 0
#define hinf1 hinf1_it2
 double hinf1 = 69.7;
#define hx hx_it2
 double hx = 1.5;
#define h_inf_it2 _thread1data[0]
#define h_inf _thread[_gth]._pval[0]
#define minf1 minf1_it2
 double minf1 = 46.2;
#define mx mx_it2
 double mx = 3;
#define m_inf_it2 _thread1data[1]
#define m_inf _thread[_gth]._pval[1]
#define phi_h_it2 _thread1data[2]
#define phi_h _thread[_gth]._pval[2]
#define phi_m_it2 _thread1data[3]
#define phi_m _thread[_gth]._pval[3]
#define shift1 shift1_it2
 double shift1 = -8;
#define sh sh_it2
 double sh = 6.4;
#define sm sm_it2
 double sm = 8.7;
#define shift2 shift2_it2
 double shift2 = 0;
#define tauh5 tauh5_it2
 double tauh5 = 46.9;
#define tauh4 tauh4_it2
 double tauh4 = 362.9;
#define tauh3 tauh3_it2
 double tauh3 = 9.1;
#define tauh2 tauh2_it2
 double tauh2 = 22.2;
#define tauh1 tauh1_it2
 double tauh1 = 21;
#define taum5 taum5_it2
 double taum5 = 13;
#define taum4 taum4_it2
 double taum4 = -0.54;
#define taum3 taum3_it2
 double taum3 = -19.7;
#define taum2 taum2_it2
 double taum2 = 125.7;
#define taum1 taum1_it2
 double taum1 = 5.4;
#define tau_h_it2 _thread1data[4]
#define tau_h _thread[_gth]._pval[4]
#define tau_m_it2 _thread1data[5]
#define tau_m _thread[_gth]._pval[5]
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "shift1_it2", "mV",
 "shift2_it2", "mV",
 "tau_m_it2", "ms",
 "tau_h_it2", "ms",
 "gcabar_it2", "mho/cm2",
 "g_it2", "mho/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "hx_it2", &hx_it2,
 "mx_it2", &mx_it2,
 "minf1_it2", &minf1_it2,
 "hinf1_it2", &hinf1_it2,
 "taum1_it2", &taum1_it2,
 "taum2_it2", &taum2_it2,
 "taum3_it2", &taum3_it2,
 "taum4_it2", &taum4_it2,
 "taum5_it2", &taum5_it2,
 "tauh1_it2", &tauh1_it2,
 "tauh2_it2", &tauh2_it2,
 "tauh3_it2", &tauh3_it2,
 "tauh4_it2", &tauh4_it2,
 "tauh5_it2", &tauh5_it2,
 "sm_it2", &sm_it2,
 "sh_it2", &sh_it2,
 "shift1_it2", &shift1_it2,
 "shift2_it2", &shift2_it2,
 "m_inf_it2", &m_inf_it2,
 "tau_m_it2", &tau_m_it2,
 "h_inf_it2", &h_inf_it2,
 "tau_h_it2", &tau_h_it2,
 "phi_m_it2", &phi_m_it2,
 "phi_h_it2", &phi_h_it2,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"it2",
 "gcabar_it2",
 0,
 "g_it2",
 0,
 "m_it2",
 "h_it2",
 0,
 0};
 static Symbol* _Ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	gcabar = 8.5e-06;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_Ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* Cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* Cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* iCa */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_diCadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _it2_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("Ca", 2.0);
 	_Ca_sym = hoc_lookup("Ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "Ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "Ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "Ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "Ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 it2 /Users/carol.h.wang/Documents/GitHub/RISE-Group-D/Biophysical_thalamocortical_system/Geir2011/x86_64/it2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double R = 8.3145;
static int _reset;
static char *modelname = "Low threshold calcium current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int evaluate_fct(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int castate(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   evaluate_fct ( _threadargscomma_ v ) ;
   Dm = ( m_inf - m ) / tau_m ;
   Dh = ( h_inf - h ) / tau_h ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 evaluate_fct ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_m )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_h )) ;
  return 0;
}
 /*END CVODE*/
 static int castate (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   evaluate_fct ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_m)))*(- ( ( ( m_inf ) ) / tau_m ) / ( ( ( ( - 1.0 ) ) ) / tau_m ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_h)))*(- ( ( ( h_inf ) ) / tau_h ) / ( ( ( ( - 1.0 ) ) ) / tau_h ) - h) ;
   }
  return 0;
}
 
static int  evaluate_fct ( _threadargsprotocomma_ double _lv ) {
   m_inf = 1.0 / ( 1.0 + exp ( - ( _lv + shift1 + minf1 ) / sm ) ) ;
   h_inf = 1.0 / ( 1.0 + exp ( ( _lv + shift2 + hinf1 ) / sh ) ) ;
   tau_m = ( taum1 + 1.0 / ( exp ( ( _lv + shift1 + taum2 ) / ( taum3 ) ) + exp ( ( _lv + shift1 + taum4 ) / taum5 ) ) ) / phi_m ;
   tau_h = ( tauh1 + 1.0 / ( exp ( ( _lv + shift2 + tauh2 ) / tauh3 ) + exp ( - ( _lv + shift2 + tauh4 ) / tauh5 ) ) ) / phi_h ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double ghk ( _threadargsprotocomma_ double _lv , double _lCi , double _lCo ) {
   double _lghk;
 double _lz , _leci , _leco ;
 _lz = ( 1e-3 ) * 2.0 * FARADAY * _lv / ( R * ( celsius + 273.15 ) ) ;
   _leco = _lCo * efun ( _threadargscomma_ _lz ) ;
   _leci = _lCi * efun ( _threadargscomma_ - _lz ) ;
   _lghk = ( .001 ) * 2.0 * FARADAY * ( _leci - _leco ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ghk ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double efun ( _threadargsprotocomma_ double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  efun ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  Cai = _ion_Cai;
  Cao = _ion_Cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  Cai = _ion_Cai;
  Cao = _ion_Cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(6, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(Datum* _thread) {
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   
/*VERBATIM*/
	Cai = _ion_Cai;
	Cao = _ion_Cao;
 phi_m = pow( mx , ( ( celsius - 23.5 ) / 10.0 ) ) ;
   phi_h = pow( hx , ( ( celsius - 23.5 ) / 10.0 ) ) ;
   evaluate_fct ( _threadargscomma_ v ) ;
   m = m_inf ;
   h = h_inf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  Cai = _ion_Cai;
  Cao = _ion_Cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = gcabar * m * m * h ;
   iCa = g * ghk ( _threadargscomma_ v , Cai , Cao ) ;
   }
 _current += iCa;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  Cai = _ion_Cai;
  Cao = _ion_Cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _diCa;
  _diCa = iCa;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_diCadv += (_diCa - iCa)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_iCa += iCa ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  Cai = _ion_Cai;
  Cao = _ion_Cao;
 {   castate(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/carol.h.wang/Documents/GitHub/RISE-Group-D/Biophysical_thalamocortical_system/Geir2011/it2.mod";
static const char* nmodl_file_text = 
  "TITLE Low threshold calcium current\n"
  ":\n"
  ":   Ca++ current responsible for low threshold spikes (LTS)\n"
  ":\n"
  ":   Written by Alain Destexhe, Salk Institute, Sept 18, 1992\n"
  ":   Modified by Geir Halnes, Norwegian University of Life Sciences, June 2011:\n"
  ":\n"
  ":     - Kinetics adapted to LGN interneuron data from Broicher et al.: Mol Cell Neurosci 36: 132-145, 2007.\n"
  ":         using Q10 values of 3 and 1.5 for activation/inactivation.\n"
  ":     - Activation variable shifted 8mV to account for dLGN interneuron data in Halnes et al. 2011\n"
  "\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "THREADSAFE\n"
  "	SUFFIX it2\n"
  "	USEION Ca READ Cai, Cao WRITE iCa VALENCE 2\n"
  "	RANGE gcabar, g\n"
  "	GLOBAL m_inf, tau_m, h_inf, tau_h, shift2, sm, sh, phi_m, phi_h, hx, mx, shift1\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mV) =	(millivolt)\n"
  "	(mA) =	(milliamp)\n"
  "	(mM) =	(millimolar)\n"
  "	FARADAY = (faraday) (coulomb)\n"
  "	R = (k-mole) (joule/degC)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER {\n"
  "	v		(mV)\n"
  "	celsius	= 36	(degC)\n"
  "	gcabar	= 8.5e-6	(mho/cm2)\n"
  "      hx      = 1.5\n"
  "      mx      = 3.0\n"
  "	Cai	= 5e-5 (mM) : Initial Ca concentration\n"
  "	Cao	= 2	(mM) : External Ca concentration\n"
  "\n"
  ": GH, parameters fitted to Broicher et al. 07 - data\n"
  "	minf1 = 46.2\n"
  "	hinf1 = 69.7\n"
  "	taum1 = 5.4\n"
  "	taum2 = 125.7\n"
  "	taum3 = -19.7\n"
  "	taum4 = -0.54\n"
  "	taum5 = 13\n"
  "	tauh1 = 21\n"
  "	tauh2 = 22.2\n"
  "	tauh3 = 9.1\n"
  "	tauh4 = 362.9\n"
  "	tauh5 = 46.9\n"
  "      sm = 8.7\n"
  "      sh = 6.4\n"
  "	shift1 = -8 	(mV) : Halnes et al. 2011\n"
  "      shift2  = 0    	(mV) : Halnes et al. 2011\n"
  "}\n"
  "\n"
  "\n"
  "STATE {\n"
  "	m h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	iCa	(mA/cm2)\n"
  "	g       (mho/cm2)\n"
  "	carev	(mV)\n"
  "	m_inf\n"
  "	tau_m	(ms)\n"
  "	h_inf\n"
  "	tau_h	(ms)\n"
  "	phi_m\n"
  "	phi_h\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE castate METHOD cnexp\n"
  "	g = gcabar * m*m*h\n"
  "	iCa = g * ghk(v, Cai, Cao)\n"
  "}\n"
  "\n"
  "DERIVATIVE castate {\n"
  "	evaluate_fct(v)\n"
  "	m' = (m_inf - m) / tau_m\n"
  "	h' = (h_inf - h) / tau_h\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "INITIAL {\n"
  "	VERBATIM\n"
  "	Cai = _ion_Cai;\n"
  "	Cao = _ion_Cao;\n"
  "	ENDVERBATIM\n"
  ":\n"
  "	phi_m = mx ^ ((celsius-23.5)/10)\n"
  "	phi_h = hx ^ ((celsius-23.5)/10)\n"
  "\n"
  "	evaluate_fct(v)\n"
  "	m = m_inf\n"
  "	h = h_inf\n"
  "}\n"
  "\n"
  "PROCEDURE evaluate_fct(v(mV)) { \n"
  "	m_inf = 1.0 / ( 1 + exp(-(v+shift1+minf1)/sm) )\n"
  "	h_inf = 1.0 / ( 1 + exp((v+shift2+hinf1)/sh) )\n"
  "	tau_m = (taum1+1.0/(exp((v+shift1+taum2)/(taum3))+exp((v+shift1+taum4)/taum5)))/ phi_m\n"
  "	tau_h = (tauh1+1/(exp((v+shift2+tauh2)/tauh3)+exp(-(v+shift2+tauh4)/tauh5)))/phi_h\n"
  "}\n"
  "\n"
  "FUNCTION ghk(v(mV), Ci(mM), Co(mM)) (.001 coul/cm3) {\n"
  "	LOCAL z, eci, eco\n"
  "	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))\n"
  "	eco = Co*efun(z)\n"
  "	eci = Ci*efun(-z)\n"
  "	:high Cao charge moves inward\n"
  "	:negative potential charge moves inward\n"
  "	ghk = (.001)*2*FARADAY*(eci - eco)\n"
  "}\n"
  "\n"
  "FUNCTION efun(z) {\n"
  "	if (fabs(z) < 1e-4) {\n"
  "		efun = 1 - z/2\n"
  "	}else{\n"
  "		efun = z/(exp(z) - 1)\n"
  "	}\n"
  "}\n"
  "UNITSON\n"
  ;
#endif
