type eln = float
exception Numerical_artifact of string*fpclass
let log_zero = nan

let eexp x =
	match classify_float x with
		| FP_normal | FP_zero -> exp x
		| FP_nan -> 0.
		| sux -> raise (Numerical_artifact ("eexp",sux))

let eln x = (* log x *)
	match classify_float x with
		| FP_zero -> log_zero
		| FP_normal when x > 0.0 -> log x
		| FP_normal -> raise (Numerical_artifact ("eln of a negative number",FP_normal))
		| sux -> raise (Numerical_artifact (Printf.sprintf "eln %F" x,sux))

let is_log_zero (x:float) = (x <> x)

let elninv eln_x =
	match classify_float eln_x with
		| FP_nan -> raise (Numerical_artifact ("elninv divide-by-zero",FP_nan)) 
		| FP_normal | FP_zero -> 0. -. eln_x
		| sux -> raise (Numerical_artifact ("elninv",sux))


external elnsum : float -> float -> float = "caml_elnsum" "elnsum" "float"
external elnproduct : float -> float -> float = "caml_elnproduct" "elnproduct" "float"
external elnproduct3 : float -> float -> float -> float = "caml_elnproduct3" "elnproduct3" "float"
external elnproduct4 : float -> float -> float -> float -> float = "caml_elnproduct4" "elnproduct4" "float"
external elnsumprod : float -> float -> float -> float = "caml_elnsumprod" "elnsumprod" "float"
external elnsumprod3 : float -> float -> float -> float -> float = "caml_elnsumprod3" "elnsumprod3" "float"

let elncompare eln_x eln_y =
	match (classify_float eln_x,classify_float eln_y) with
		| (c1,c2) when (c1 = FP_normal || c1 = FP_zero) && (c2 = FP_normal || c2 = FP_zero) ->
			compare eln_x eln_y
		| (c,FP_nan) when (c = FP_normal || c = FP_zero) -> 1
		| (FP_nan,c) when (c = FP_normal || c = FP_zero) -> (-1)
		| (FP_nan,FP_nan) -> 0
		| ((FP_subnormal as sux),_) | (_,(FP_subnormal as sux))
		| ((FP_infinite as sux),_) | (_,(FP_infinite as sux)) ->
			raise (Numerical_artifact ("elncompare",sux))

let elndiff eln_x eln_y =
	if is_log_zero eln_y then
		eln_x
	else if (elncompare eln_x eln_y) >= 0 then
		eln_x +. log (1. -. exp (eln_y -. eln_x))
	else
		raise (Numerical_artifact ("elndiff negative answer",FP_nan))

