(*pp ocaml+twt *)
open Printf

let print_matrix mat =
	let (m,n) = Gsl_matrix.dims mat
	for row = 0 to (m-1) do
		for col = 0 to (n-1) do
			eprintf " %8.5f" mat.{row,col}
		eprintf "\n"

let print_caml_matrix mat =
	let m = Array.length mat
	for row = 0 to (m-1) do
		let n = Array.length mat.(row)
		for col = 0 to (n-1) do
			eprintf " %8.5f" mat.(row).(col)
		eprintf "\n"


let print_vector vec =
	let n = Gsl_vector.length vec
	for i = 0 to n-1 do
		eprintf " %8.5f" vec.{i}
	eprintf "\n"

let print_caml_vector vec =
	let n = Array.length vec
	for i = 0 to n-1 do
		eprintf " %8.5f" vec.(i)
	eprintf "\n"

let print_intvector vec =
	let n = Bigarray.Array1.dim vec
	for i = 0 to n-1 do
		eprintf " %2d" vec.{i}
	eprintf "\n"

let matrix_set_all (m:float array array) x =
	for i = 0 to (Array.length m)-1 do
		for j = 0 to (Array.length m.(i))-1 do
			m.(i).(j) <- x

let array_index_of a x =
	let l = Array.length a
	let i = ref 0
	while !i < l && a.(!i) <> x do
		i := !i + 1
	if !i = l then raise Not_found
	!i

let rec range lower upper =
	if lower >= upper then
		[]
	else
		lower :: (range (lower+1) upper)

module Expando = struct
	type 'a t = { mutable a : 'a array; default : 'a }
	let make l default = { a = Array.make l default; default = default }
	let access e = e.a
	let ensure e l_new =
		let l_old = Array.length e.a
		if l_new > l_old then
			let a_new = Array.init (max (l_old*2) l_new) (fun i -> if i < l_old then e.a.(i) else e.default)
			e.a <- a_new
	let access_ensure e l_new =
		ensure e l_new
		access e

module Rotator = struct
	type 'a t = { mutable zero : int; items : 'a array }
	let make items = { zero = 0; items = items }
	let rotate r =
		let l = Array.length r.items
		if r.zero = 0 then
			r.zero <- (l-1)
		else
			r.zero <- r.zero-1
	let get r i =
		if i >= (Array.length r.items) then invalid_arg "CRF.SMCRF.Rotator.get"
		r.items.((r.zero + i) mod (Array.length r.items))

