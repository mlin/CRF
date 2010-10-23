#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>

#ifndef NAN
#error "wtf"
#endif

#define MYISNAN(x) (x!=x) /*isnan(x)*/

inline double elnsum(double eln_x,double eln_y) {
  if(MYISNAN(eln_x))
	return eln_y;
  if(MYISNAN(eln_y))
	return eln_x;

  if(eln_x >= eln_y)
	return (eln_x + log(1.0 + exp(eln_y - eln_x)));
  else
	return (eln_y + log(1.0 + exp(eln_x - eln_y)));
}

value caml_elnsum(value eln_x,value eln_y) {
  CAMLparam2(eln_x,eln_y);
  double ans = elnsum(Double_val(eln_x),Double_val(eln_y));
  CAMLreturn(caml_copy_double(ans));
}

inline double elnproduct(double eln_x,double eln_y) {
  /* The following happens automatically:
  if(MYISNAN(eln_x)||MYISNAN(eln_y))
	return NAN;
  */
  return (eln_x+eln_y);
}

value caml_elnproduct(value eln_x,value eln_y) {
  CAMLparam2(eln_x,eln_y);
  double ans = elnproduct(Double_val(eln_x),Double_val(eln_y));
  CAMLreturn(caml_copy_double(ans));
}

inline double elnproduct3(double eln_x,double eln_y,double eln_z) {
  return (eln_x+eln_y+eln_z);
}

value caml_elnproduct3(value eln_x,value eln_y,value eln_z) {
  CAMLparam3(eln_x,eln_y,eln_z);
  double ans = elnproduct3(Double_val(eln_x),Double_val(eln_y),Double_val(eln_z));
  CAMLreturn(caml_copy_double(ans));
}

inline double elnproduct4(double eln_w,double eln_x,double eln_y,double eln_z) {
  return (eln_w+eln_x+eln_y+eln_z);
}

value caml_elnproduct4(value eln_w,value eln_x,value eln_y,value eln_z) {
  CAMLparam4(eln_w,eln_x,eln_y,eln_z);
  double ans = elnproduct4(Double_val(eln_w),Double_val(eln_x),Double_val(eln_y),Double_val(eln_z));
  CAMLreturn(caml_copy_double(ans));
}

inline double elnsumprod(double eln_a,double eln_b,double eln_c) {
  return(elnsum(elnproduct(eln_a,eln_b),eln_c));
}

value caml_elnsumprod(value eln_a,value eln_b,value eln_c) {
  CAMLparam3(eln_a,eln_b,eln_c);
  double ans = elnsumprod(Double_val(eln_a),Double_val(eln_b),Double_val(eln_c));
  CAMLreturn(caml_copy_double(ans));
}

inline double elnsumprod3(double eln_a,double eln_b,double eln_c,double eln_d) {
  return(elnsum(elnproduct3(eln_a,eln_b,eln_c),eln_d));
}

value caml_elnsumprod3(value eln_a,value eln_b,value eln_c,value eln_d) {
  CAMLparam4(eln_a,eln_b,eln_c,eln_d);
  double ans = elnsumprod3(Double_val(eln_a),Double_val(eln_b),Double_val(eln_c),Double_val(eln_d));
  CAMLreturn(caml_copy_double(ans));
}
