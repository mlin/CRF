(** Floating-point operations intended for working with small numbers (e.g. probabilities) in log-space.

These operations provide a special representation for log(zero), since it is often necessary to deal with zero probabilities. They are inspired by Tobias P. Mann's tutorial, "Numerically Stable Hidden Markov Model Implementation". *)

type eln = float

(** Raised when any function encounters an [FP_subnormal] or [FP_infinite] as indicated by OCaml's [classify_float]. *)
exception Numerical_artifact of string*fpclass

(** The special value indicating log(zero). The implementation uses [nan] for this value meaning that you {b cannot} use [=] or [<>] to test for [log_zero]. *)
val log_zero : eln

(** Tests whether the value is log(zero) *)
val is_log_zero : eln -> bool

(** [eexp(x)] is [exp(x)] for real [x], and 0.0 for [x = log_zero]. *)
val eexp : float -> float

(** [eln(x)] is [ln(x)] for positive real [x], and [log_zero] for [x = 0]. *)
val eln : float -> eln

(** [elninv(eln(x))] is [1/ln(x)] for nonzero [x].
	@raise Numerical_artifact if [x = 0]. *)
val elninv : eln -> eln

(** [elnsum(eln(x),eln(y))] is [ln(x+y)] for nonnegative real [x] and [y]. *)
val elnsum : eln -> eln -> eln

(** [elndiff(eln(x),eln(y))] is [ln(x-y)] for nonnegative real [x] and [y], [x >= y].
    @raise Numerical_artifact if [x < y] *)
val elndiff : eln -> eln -> eln

(** [elnproduct(eln(x),eln(y)) = ln(xy)] for positive [x] and [y] ([log_zero] if either [x] or [y] is zero) *)
val elnproduct : eln -> eln -> eln

(** [elnproduct(eln(x),eln(y),eln(z)) = ln(xyz)] *)
val elnproduct3 : eln -> eln -> eln -> eln

val elnproduct4 : eln -> eln -> eln -> eln -> eln

(** [elnsumprod(eln(a),eln(b),eln(c)) = ln(ab+c)] *)
val elnsumprod : eln -> eln -> eln -> eln

(** [elnsumprod(eln(a),eln(b),eln(c),eln(d)) = ln(abc+d)] *)
val elnsumprod3 : eln -> eln -> eln -> eln -> eln

val elncompare : eln -> eln -> int
