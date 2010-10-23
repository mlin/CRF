(*pp ocaml+twt *)

module K = struct
	type t = int
	let compare = compare
module M = Map.Make(K)
		
type t = {
	length : int;
	map : float M.t;
}
	
let make l =
	if l < 0 then invalid_arg "SparseVector.make"
	{ length = l; map = M.empty }
let length { length = x } = x
let get v k =
	if k < 0 || k >= v.length then invalid_arg "SparseVector.get"
	try M.find k v.map with Not_found -> 0.
let set v k x =
	if k < 0 || k >= v.length then invalid_arg "SparseVector.set"
	{ length = v.length; map = M.add k x v.map }
let grow v l = if l > v.length then { length = l; map = v.map } else v
	
let clean ?(threshold=0.) v =
	let m = 
		M.fold
			fun k x m -> if abs_float x > threshold then M.add k x m else m
			v.map
			M.empty
	{ length = v.length; map  = m }

let sub v lo n =
	if lo < 0 || n < 0 then invalid_arg "SparseVector.sub"
	let hi = lo+n-1
	let m =
		M.fold
			fun k x m -> if k >= lo && k <= hi then M.add (k-lo) x m else m
			v.map
			M.empty
	{ length = n; map = m }

let to_array v =
	let a = Array.create v.length 0.
	M.iter (fun k x -> Array.unsafe_set a k x) v.map
	a
	
let of_array a =
	let l = Array.length a
	let m = ref M.empty
	for k = 0 to l-1 do
		let x = Array.unsafe_get a k
		if x <> 0. then m := M.add k x !m
	{ length = l; map = !m }

let count_set v = M.fold (fun _ _ z -> z+1) v.map 0
let iter_set f v = M.iter (fun _ x -> f x) v.map
let iteri_set f v = M.iter f v.map
let map_set f v = { length = v.length; map = M.map f v.map }
let mapi_set f v = { length = v.length; map = M.mapi f v.map }
let foldi_set f v z = M.fold f v.map z
let fold_set f v z = M.fold (fun _ x z -> f x z) v.map z

let norm ?(l=2) v =
	if l < 1 then invalid_arg "SparseArray.norm"
	let l = float l
	(M.fold (fun _ x z -> z +. x ** l) v.map 0.) ** (1. /. l)

let addf k x m = M.add k (try x +. (M.find k m) with Not_found -> x) m
let add v1 v2 =
	if v1.length <> v2.length then invalid_arg "SparseVector.add: length mismatch"
	{ length = v1.length; map = M.fold addf v1.map v2.map }
	
(* nts: carefully document semantics for dot and dot_array in case of any entries being inf or -inf (zero*infinity=zero) *)
let dot v1 v2 =
	if v1.length <> v2.length then invalid_arg "SparseVector.dot: length mismatch"
	M.fold (fun k x z -> try z +. (M.find k v2.map) *. x with Not_found -> z) v1.map 0.
	
let dot_array v a =
	if Array.length a <> v.length then invalid_arg "SparseVector.dot_array: length mismatch"
	M.fold (fun k x z -> z +. x *. (Array.unsafe_get a k)) v.map 0.

let add_to_array ?(ofs=0) v a =
	if Array.length a < ofs+v.length then invalid_arg "SparseVector.add_to_array: length mismatch"
	M.iter (fun k x -> a.(ofs+k) <- a.(ofs+k) +. x) v.map	
	
type double_bigarray = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
let dot_bigarray ?(ofs=0) v ba =
	if Bigarray.Array1.dim ba < ofs+v.length then invalid_arg "SparseVector.dot_bigarray: length mismatch"
	M.fold (fun k x z -> z +. x *. ba.{ofs+k}) v.map 0.

let add_to_bigarray ?(ofs=0) v ba =
	if Bigarray.Array1.dim ba < ofs+v.length then invalid_arg "SparseVector.add_to_bigarray: length mismatch"
	M.iter (fun k x -> ba.{ofs+k} <- ba.{ofs+k} +. x) v.map

let append v1 v2 =
	let m =
		M.fold
			fun k x m -> M.add (k+v1.length) x m
			v2.map
			v1.map
	{ length = v1.length+v2.length; map = m }

let concat vs =
	if vs <> [] then
		List.fold_left append (List.hd vs) (List.tl vs)
	else
		make 0
