(*pp ocaml+twt *)

(** SMCRF implements core algorithms for chain-structured, semi-Markov conditional random fields. *)

open Printf
open Elog
open Tools

(******************************************************************************)
(* types *)
(******************************************************************************)

(** labels are specified by number. label 0 is always start. The highest valid label is always end. *)
type label = int

(** "Settings" for each user-specified label. *)
type label_spec = {
	(** The maximum semi-Markov segment length for this label *)
	max_length : int
}

type segment = {
	seg_label : label;
	seg_end : int;
	seg_length : int
}

type segmentation = segment list

(** The type of feature functions. The algorithm can use the information that a feature function is independent of certain parameters to cache values, reducing the number of times it has to call that feature function. Feature functions may return any finite real value. *)
type 'obs_seq feature_function =
	| SemiMarkov_feature of		(prev_label:int ->	label:int ->	seq:'obs_seq ->	i:int ->	seg_length:int ->	chain_length:int ->	float)*[`NoCache]
	| Segment_feature of		(					label:int ->	seq:'obs_seq ->	i:int ->	seg_length:int ->	chain_length:int ->	float)*[`NoCache]
	| Markov_feature of			(prev_label:int ->	label:int ->	seq:'obs_seq -> i:int ->						chain_length:int ->	float)*[`NoCache|`DenseCache|`SparseCache]
	| Emission_feature of		(					label:int ->	seq:'obs_seq -> i:int -> 						chain_length:int -> float)*[`NoCache|`DenseCache|`SparseCache]
	| Length_feature of			(					label:int ->	seq:'obs_seq ->				seg_length:int ->	chain_length:int -> float)
	| Transition_feature of		(prev_label:int ->	label:int ->	seq:'obs_seq ->									chain_length:int -> float)

	| SemiMarkov_feature_vec of	(prev_label:int -> label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> values:(float array) -> unit)*[`NoCache|`DenseCache|`SparseCache]
	| Segment_feature_vec of	(label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> values:(float array) -> unit)*[`NoCache|`DenseCache|`SparseCache]

type 'obs_seq sparse_feature_set =
	| Emission_features of int*(label:int -> seq:'obs_seq -> i:int -> chain_length:int -> SparseVector.t)
	| Markov_features of int*(prev_label:int -> label:int -> seq:'obs_seq -> i:int -> chain_length:int -> SparseVector.t)
	| Segment_features_vec of int*(label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> SparseVector.t array)
	| SemiMarkov_features_vec of int*(prev_label:int -> label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> SparseVector.t array)

(* TODO caching of length restrictor *)
type 'obs_seq seg_length_restrictor = (prev_label:int -> label:int -> seq:'obs_seq -> i:int -> max_seg_length:int -> chain_length:int -> fwd:bool -> int array)

let d0 = [||]
let d1 = [| 1 |]
let default_seg_length_restrictor ~prev_label ~label ~seq ~i ~max_seg_length ~chain_length ~fwd =
	if (not fwd) && prev_label = 0 then
		if (i+1) <= max_seg_length then
			[| i+1 |]
		else
			d0
	else
		let maxd = if fwd then (min max_seg_length (chain_length-i-1)) else (min max_seg_length (i+1))
		if maxd = 1 then
			d1
		else
			Array.init maxd (fun d -> d+1)

(** The abstract type of a trained CRF. *)
type 'obs_seq t = {
	n_user_labels : int;
	label_specs : label_spec array;
	max_max_length : int;
	feature_functions : 'obs_seq feature_function array;
	sparse_feature_sets : 'obs_seq sparse_feature_set array;
	seg_length_restrictor : 'obs_seq seg_length_restrictor;
	lambda : Gsl_vector.vector
}

exception False

(* an internal data structure for caching and retrieving feature values during CRF dynamic programming algorithms*)
type 'obs_seq featcache = {
	crf : 'obs_seq t;
	chain_length : int;
	seq : 'obs_seq;
	cache : float array array array;
	vec_cache : float array array array array;
	sparse_cache : (int*int*int*int,float) Hashtbl.t;
	sparse_vec_cache : (int*int*int*int*bool,float array) Hashtbl.t;
	mutable cache_primed_forward : bool;
	mutable cache_primed_backward : bool;
	transitions : int array array;
	mutable transition_ct : int;
	transition_map : int array
}

type 'obs_seq training_sequence = {
	train_seq : 'obs_seq;
	train_chain_length : int;
	train_segs : segmentation;
	train_id : string
}

(******************************************************************************)
(* Feature cache implementation *)
(******************************************************************************)

(* feature memoization *)
let markov_lookup fc featnum f prev_label label i =
	if (Array.length fc.cache.(featnum)) = 0 then
		fc.cache.(featnum) <- Array.make_matrix (fc.transition_ct+fc.crf.n_user_labels+2) (fc.chain_length+1) infinity
	let trn = if prev_label=label then (fc.transition_ct+label) else fc.transition_map.((fc.crf.n_user_labels+2)*prev_label+label) (* hack needed because Markov_features with d>1 always use the self-loop transition *)
	assert (trn >= 0)
	let cacheline = fc.cache.(featnum).(trn)
	if fc.cache_primed_forward || fc.cache_primed_backward then
		cacheline.(i)
	else
		let v = f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length
		cacheline.(i) <- v
		v

let markov_sparse_lookup fc featnum f prev_label label i =
	if fc.cache_primed_forward || fc.cache_primed_backward then
		try
			Hashtbl.find fc.sparse_cache (i,featnum,prev_label,label)
		with
			| Not_found -> 0.
	else
		let v = f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length
		if v <> 0. then
			Hashtbl.replace fc.sparse_cache (i,featnum,prev_label,label) v
		v

let emission_lookup fc featnum f label i =
	if (Array.length fc.cache.(featnum)) = 0 then
		fc.cache.(featnum) <- Array.make_matrix (fc.crf.n_user_labels+2) (fc.chain_length+1) infinity
	let cacheline = fc.cache.(featnum).(label)
	if fc.cache_primed_forward || fc.cache_primed_backward then
		cacheline.(i)
	else
		let v = f ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length
		cacheline.(i) <- v
		v

let length_lookup fc featnum f label seg_length =
	assert (seg_length > 0 && seg_length <= fc.crf.max_max_length)
	if (Array.length fc.cache.(featnum)) = 0 then
		fc.cache.(featnum) <- Array.make_matrix 1 fc.crf.max_max_length infinity
	let cacheline = fc.cache.(featnum).(0)
	if fc.cache_primed_forward || fc.cache_primed_backward then
		cacheline.(seg_length-1)
	else
		let v = f ~label:label ~seq:fc.seq ~seg_length:seg_length ~chain_length:fc.chain_length
		cacheline.(seg_length-1) <- v
		v

let transition_lookup fc featnum f prev_label label =
	if (Array.length fc.cache.(featnum)) = 0 then
		fc.cache.(featnum) <- Array.make_matrix (fc.crf.n_user_labels+2) (fc.crf.n_user_labels+2) infinity
	let cacheline = fc.cache.(featnum).(prev_label)
	if fc.cache_primed_forward || fc.cache_primed_backward then
		cacheline.(label)
	else
		let v = f ~prev_label:prev_label ~label:label ~seq:fc.seq ~chain_length:fc.chain_length
		cacheline.(label) <- v
		v
			
let semimarkov_vec_lookup fc featnum f prev_label label i seg_lengths fwd values =
	let nds = Array.length seg_lengths
	if (Array.length fc.vec_cache.(featnum)) = 0 then
		fc.vec_cache.(featnum) <- Array.make_matrix (2*fc.transition_ct) (fc.chain_length+2) [||]
	let trn = fc.transition_map.((fc.crf.n_user_labels+2)*prev_label+label) + (if fwd then fc.transition_ct else 0)
	assert (trn >= 0)
	let cacheline = fc.vec_cache.(featnum).(trn)
	if ((not fwd) && fc.cache_primed_forward) || (fwd && fc.cache_primed_backward) then
		let a = cacheline.(i+1) (* +1 because on a backwards pass there will be an i=-1 *)
		assert (Array.length a = nds)
		Array.blit a 0 values 0 nds
	else
		let a = Array.make nds nan
		f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~seg_lengths:seg_lengths ~chain_length:fc.chain_length ~fwd:fwd ~values:a
		cacheline.(i+1) <- a
		Array.blit a 0 values 0 nds
		
let segment_vec_lookup fc featnum f label i seg_lengths fwd values =
	let nds = Array.length seg_lengths
	if (Array.length fc.vec_cache.(featnum)) = 0 then
		fc.vec_cache.(featnum) <- Array.make_matrix (2*(fc.crf.n_user_labels+2)) (fc.chain_length+2) [||]
	let cacheline = fc.vec_cache.(featnum).(label + (if fwd then fc.crf.n_user_labels+2 else 1))
	if ((not fwd) && fc.cache_primed_forward) || (fwd && fc.cache_primed_backward) then
		let a = cacheline.(i+1)
		assert (Array.length a = nds)
		Array.blit a 0 values 0 nds
	else
		let a = Array.make nds nan
		f ~label:label ~seq:fc.seq ~i:i ~seg_lengths:seg_lengths ~chain_length:fc.chain_length ~fwd:fwd ~values:a
		cacheline.(i+1) <- a
		Array.blit a 0 values 0 nds

let sparse_vec_lookup fc featnum f prev_label label i seg_lengths fwd values =
	let nds = Array.length seg_lengths
	if nds>0 then
		let id = (i,featnum,prev_label,label,fwd)
		if (fwd && fc.cache_primed_backward) || ((not fwd) && fc.cache_primed_forward) then
			let a = Hashtbl.find fc.sparse_vec_cache id
			(* this is dangerous because we're assuming seg_lengths will always be the same whenever the feature is evaluated here... *)
			assert (Array.length a = nds)
			Array.blit a 0 values 0 nds
		else
			let a = Array.make nds nan
			f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~seg_lengths:seg_lengths ~chain_length:fc.chain_length ~fwd:fwd ~values:a
			Hashtbl.replace fc.sparse_vec_cache id a
			Array.blit a 0 values 0 nds

let box = [| 0. |]

(* return the value of the feature function featnum applied to a segment ending
	at i, with the given label, previous label, and length=1 *)
let feature_value fc featnum prev_label label i =
	match fc.crf.feature_functions.(featnum) with
		| Length_feature f -> length_lookup fc featnum f label 1
		| Transition_feature f -> transition_lookup fc featnum f prev_label label
		| Markov_feature (f,`DenseCache) -> markov_lookup fc featnum f prev_label label i
		| Markov_feature (f,`SparseCache) -> markov_sparse_lookup fc featnum f prev_label label i
		| Markov_feature (f,_) -> f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length
		| Emission_feature (f,`DenseCache) -> emission_lookup fc featnum f label i
		| Emission_feature (f,`SparseCache) -> markov_sparse_lookup fc featnum (fun ~prev_label -> f) 0 label i
		| Emission_feature (f,_) -> f ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length
		| SemiMarkov_feature (f,_) -> f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~seg_length:1 ~chain_length:fc.chain_length
		| Segment_feature (f,_) -> f ~label:label ~seq:fc.seq ~i:i ~seg_length:1 ~chain_length:fc.chain_length
		(* no caching of vec features because we may only be looking at a subset of the lengths *)
		| SemiMarkov_feature_vec (f,_) ->
			f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~seg_lengths:d1 ~chain_length:fc.chain_length ~fwd:false ~values:box
			box.(0)
		| Segment_feature_vec (f,_) ->
			f ~label:label ~seq:fc.seq ~i:i ~seg_lengths:d1 ~chain_length:fc.chain_length ~fwd:false ~values:box
			box.(0)

let meta_segments_markov_sums reduce null fwd f i seg_lengths values =
	let nds = Array.length seg_lengths
	if nds>0 then
		let v = ref null
		let id = ref 0
		for d = 1 to seg_lengths.(nds-1) do
			assert (!id < nds && d <= seg_lengths.(!id))
			if d = seg_lengths.(!id) then
				(* special handling needed for the "first" position in the segment, so that the correct prev_label is passed to the feature function
                   when fwd=true, this is the "first" position in the segment only when d=1. 
				   when not fwd, we're always in the "first" position in the segment when d=seg_lengths.(!id) *)
				let first = (not fwd) || d=1
				values.(!id) <- reduce !v (f first (if fwd then (i + d) else (i - d + 1)))
				assert (!id = 0 || seg_lengths.(!id) > seg_lengths.(!id-1))
				incr id

			v := reduce !v (f (fwd && d=1) (if fwd then (i + d) else (i - d + 1)))
		assert (!id = (Array.length seg_lengths))
		
let segments_markov_sums = meta_segments_markov_sums (+.) 0.


let segments_semimarkov_values fwd f i seg_lengths values =
	Array.iteri
		fun id d ->
			values.(id) <- (if fwd then f (i+d) d else f i d)
		seg_lengths

(*
	evaluate the feature function at (y',y,i) for each segment length in
	seg_lengths. seg_lengths is an array of length < fc.crf.max_max_length
	containing unique natural numbers in increasing order. The results
	are placed in values which must be a float array of the same length.
	if fwd=false, then the segments evaluated are
		(i-seg_lengths[0],i], (i-seg_lengths[1],i], (i-seg_lengths[2],i], etc.
	if fwd=true, then the segments evaluated are
		(i,i+seg_lengths[0]], (i,i+seg_lengths[2]], (i,i+seg_lengths[3]], etc.
*)		
let segments_values ?(no_vec_cache=false) ?(fwd=false) fc featnum prev_label label i seg_lengths values =
	match fc.crf.feature_functions.(featnum) with

		| Markov_feature (f,`DenseCache) ->
			segments_markov_sums fwd (fun first i -> markov_lookup fc featnum f (if first then prev_label else label) label i) i seg_lengths values
		| Markov_feature (f,`SparseCache) ->
			segments_markov_sums fwd (fun first i -> markov_sparse_lookup fc featnum f (if first then prev_label else label) label i) i seg_lengths values
		| Markov_feature (f,_) ->
			segments_markov_sums fwd (fun first i -> f ~prev_label:(if first then prev_label else label) ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length) i seg_lengths values

		| Emission_feature (f,`DenseCache) ->
			segments_markov_sums fwd (fun _ i -> emission_lookup fc featnum f label i) i seg_lengths values
		| Emission_feature (f,`SparseCache) ->
			segments_markov_sums fwd (fun _ i -> markov_sparse_lookup fc featnum (fun ~prev_label -> f) 0 label i) i seg_lengths values
		| Emission_feature (f,_) ->
			segments_markov_sums fwd (fun _ i -> f ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length) i seg_lengths values

		| Length_feature f ->
			segments_semimarkov_values fwd (fun i d -> length_lookup fc featnum f label d) i seg_lengths values

		| Transition_feature f -> 
			let v = transition_lookup fc featnum f prev_label label
			for id = 0 to (Array.length seg_lengths)-1 do
				values.(id) <- v

		| SemiMarkov_feature (f,_) ->
			let pf = (fun i d -> f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~seg_length:d ~chain_length:fc.chain_length)
			segments_semimarkov_values fwd pf i seg_lengths values

		| Segment_feature (f,_) ->
			let pf = (fun i d -> f ~label:label ~seq:fc.seq ~i:i ~seg_length:d ~chain_length:fc.chain_length)
			segments_semimarkov_values fwd pf i seg_lengths values

		| SemiMarkov_feature_vec (f,`DenseCache) when not no_vec_cache -> semimarkov_vec_lookup fc featnum f prev_label label i seg_lengths fwd values
		| SemiMarkov_feature_vec (f,`SparseCache) when not no_vec_cache -> sparse_vec_lookup fc featnum f prev_label label i seg_lengths fwd values
		| SemiMarkov_feature_vec (f,_) -> f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~seg_lengths:seg_lengths ~chain_length:fc.chain_length ~fwd:fwd ~values:values

		| Segment_feature_vec (f,`DenseCache) when not no_vec_cache -> segment_vec_lookup fc featnum f label i seg_lengths fwd values
		| Segment_feature_vec (f,`SparseCache) when not no_vec_cache -> sparse_vec_lookup fc featnum (fun ~prev_label -> f) (-314159) label i seg_lengths fwd values
		| Segment_feature_vec (f,_) -> f ~label:label ~seq:fc.seq ~i:i ~seg_lengths:seg_lengths ~chain_length:fc.chain_length ~fwd:fwd ~values:values

let sparse_segments_values ?(fwd=false) fc sfeatnum prev_label label i seg_lengths =
	match fc.crf.sparse_feature_sets.(sfeatnum) with
		| Emission_features (n,f) ->
			let s0 = SparseVector.make n
			let sas = Array.make (Array.length seg_lengths) s0
			meta_segments_markov_sums SparseVector.add s0 fwd (fun _ i -> f ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length) i seg_lengths sas
			sas
		| Markov_features (n,f) ->
			let s0 = SparseVector.make n
			let sas = Array.make (Array.length seg_lengths) s0
			meta_segments_markov_sums SparseVector.add s0 fwd (fun first i -> f ~prev_label:(if first then prev_label else label) ~label:label ~seq:fc.seq ~i:i ~chain_length:fc.chain_length) i seg_lengths sas
			sas
		| Segment_features_vec (n,f) -> f ~label:label ~seq:fc.seq ~i:i ~seg_lengths:seg_lengths ~chain_length:fc.chain_length ~fwd:fwd
		| SemiMarkov_features_vec (n,f) -> f ~prev_label:prev_label ~label:label ~seq:fc.seq ~i:i ~seg_lengths:seg_lengths ~chain_length:fc.chain_length ~fwd:fwd

(* featcache setup *)
let setup_transitions fc =
	let n = (fc.crf.n_user_labels+2)
	let labels = range 0 n
	let ntr = ref 0
	for y = 0 to n-1 do
		let a =
			Array.of_list
				List.filter
					fun y' ->
						if y <> 0 && y' <> n-1 then
							try
								for i = 0 to (Array.length fc.crf.feature_functions) - 1 do
									match fc.crf.feature_functions.(i) with
										| Transition_feature f -> if is_log_zero (transition_lookup fc i f y' y) then raise False
										| _ -> ()
								let trnum = y'*n+y
								fc.transition_map.(trnum) <- !ntr
								incr ntr
								true
							with False -> false
						else
							false
					labels
		fc.transitions.(y) <- a
	fc.transition_ct <- !ntr
			

let make_featcache ~crf ~chain_length ~seq =
	let fc = {
		crf = crf;
		chain_length = chain_length;
		seq = seq;
		cache = Array.make (Array.length crf.feature_functions) (Array.make_matrix 0 0 infinity);
		vec_cache = Array.make (Array.length crf.feature_functions) (Array.init 0 (fun _ -> Array.make_matrix 0 0 infinity));
		sparse_cache = (Hashtbl.create 1024);
		sparse_vec_cache = (Hashtbl.create 1024);
		cache_primed_forward = false; cache_primed_backward = false;
		transitions = Array.make (crf.n_user_labels+2) (Array.make 0 0);
		transition_ct = 0;
		transition_map = Array.make ((crf.n_user_labels+2)*(crf.n_user_labels+2)) (-1)
	}
	setup_transitions fc
	fc

(******************************************************************************)
(* core algorithms *)
(******************************************************************************)

(* compute the potential matrix M for a given i
	if fwd = false, then
		M(y,y') = the segment potentials for segments ENDING at i, with label y,
					preceded by a segment with label y'.
	if fwd = true, then
		M(y,y') = the segment potentials for segments BEGINNING at i+1, with label y,
					preceded by a segment with label y'.
	segment potentials are specified as a tuple of two arrays. The first array
	contains the LENGTHS of all allowable segments ending/beginning at i/i+1,
	given y' and y. The second array contains the corresponding potentials.
	The "allowable" segments are specified by the segment length restrictor.
*)

type mi = (int array*float array) array array

let pot0 = [||]
let dpot0 = ([||],[||])
let compute_Mi ~fwd fc i =
	let nff = Array.length fc.crf.feature_functions
	let nsff = Array.length fc.crf.sparse_feature_sets
	let start_label = 0
	let end_label = fc.crf.n_user_labels+1
	let start_pos = 0
	let end_pos = fc.chain_length
	let lambda = fc.crf.lambda
	let segvals = Expando.make (min 32 fc.crf.max_max_length) log_zero

	let state_ok y' y =
		assert (y' <> end_label)
		assert (y <> start_label)
		if fwd then
			((i >= start_pos || y' = start_label) &&
			 (i < start_pos || y' <> start_label) &&
			 ((i+1) < end_pos || y = end_label) &&
			 ((i+1) = end_pos || y <> end_label) &&
			 (y' <> start_label || y <> end_label))
		else
			((i <> end_pos || y = end_label) &&
			 (i == end_pos || y <> end_label) &&
			 (i <> 0 || y' == start_label) &&
			 (y' <> start_label || y <> end_label) ) (* &&
			 (i == 0 || y' <> start_label)) wrong in case the first segment has length > 1 *)

	let mi =
		Array.init (fc.crf.n_user_labels+2)
			fun y ->
				Array.init (Array.length fc.transitions.(y))
					fun iy' ->			
						let y' = fc.transitions.(y).(iy')
						if not (state_ok y' y) then
							dpot0
						else
							let ds =
								if y <> end_label then
									(fc.crf.seg_length_restrictor
										~prev_label:y' ~label:y	~seq:fc.seq	~i:i
										~max_seg_length:fc.crf.label_specs.(y-1).max_length
										~chain_length:fc.chain_length ~fwd:fwd)
									(* TODO it would be nice to have some optional validation of ds *)
								else
									d1
							if (Array.length ds) = 0 then
								dpot0
							else
								if (not fwd) && y' = start_label then
									if (Array.length ds) > 1 || ds.(0) <> (i+1) then
										failwith (sprintf "CRF.SMCRF: segment length restrictor returned a invalid first segment(s) (y'=%d,y=%d,ds=[%s],i=%d)" y' y (String.concat " " (List.map string_of_int (Array.to_list ds))) i)
								let segvals = Expando.access_ensure segvals (Array.length ds)
								let pots = Array.make (Array.length ds) 0.
(*								if not ((Array.length ds) = 1 && ds.(0) = 1) then *)
								if true then
									for j = 0 to nff-1 do
										segments_values ~fwd:fwd fc j y' y i ds segvals
										for id = 0 to (Array.length pots)-1 do
											(* assert (if fwd then (i+ds.(id)) <= fc.chain_length else (i-ds.(id)+1) >= 0) *)
											pots.(id) <- pots.(id) +. lambda.{j} *. segvals.(id)
(*								else
									(* faster way for Markov segments (ds=[|1|]) *)
									pots.(0) <- lambda.{0} *. (feature_value fc 0 y' y (if fwd then i+1 else i))
									let j = ref 1
									while !j < nff && not (is_log_zero pots.(0)) do
										pots.(0) <- pots.(0) +. lambda.{!j} *. (feature_value fc !j y' y (if fwd then i+1 else i))
										incr j *)
									
								if nsff > 0 then
									let ssvs = Array.init nsff (fun sj -> sparse_segments_values ~fwd:fwd fc sj y' y i ds)
									for id = 0 to (Array.length pots)-1 do
										let ofs = ref nff
										for sj = 0 to nsff-1 do
											let sv = ssvs.(sj).(id)
											pots.(id) <- pots.(id) +. SparseVector.dot_bigarray ~ofs:!ofs sv lambda
											ofs := !ofs + SparseVector.length sv
										assert (!ofs = Gsl_vector.length lambda)

								(ds,pots)
	(mi:mi)

(* captures a common pattern of iterating through the sequence and performing a computation on the matrices M(i) at each position *)
let fold_Mis fc f init =
	let rslt = ref init

	(* iterate over positions *)		
	for i = 0 to fc.chain_length do
		rslt := f i !rslt (compute_Mi ~fwd:false fc i)
	fc.cache_primed_forward <- true
	!rslt

(* compute the alphas. note, a_0 is virtual (not actually represented in the matrix) *)
let logalpha_folder fc =
	let a = Gsl_matrix.create (fc.chain_length+1) (fc.crf.n_user_labels+2)
	let lafolder i a (mi:mi) =
		for y = 0 to (fc.crf.n_user_labels+1) do				(* to compute each a.{i,y} *)
			let miy = mi.(y)
			let tot = ref log_zero
			for iy' = 0 to (Array.length fc.transitions.(y)) - 1 do 			(* accumulate a.{i-d} dot M(i,y,d) *)
				let y' = fc.transitions.(y).(iy')
				let (seg_lengths,seg_potentials) = miy.(iy')

				for id = 0 to (Array.length seg_lengths) - 1 do
					let d = seg_lengths.(id)
					let logalpha_imd_y' =
						if (i-d) = -1 then
							if y' = 0 then 0.0 else log_zero
						else
							a.{i-d,y'}
					tot := elnsumprod logalpha_imd_y' seg_potentials.(id) !tot
			a.{i,y} <- !tot
		a
	(lafolder,a)

(* compute the potential for the given segmentation. the potential can be divided by the normalization factor Z(x) to obtain the probability of the segmentation *)

exception Zero_potential of string*segment*int*string
(* TODO compute only the needed Mis; big speedup with lengthy segments *)
let potential_folder tsid fc segs =
	let potfolder i ((sc,prev_label,segs) as sofar) (mi:mi) =
		assert (not (is_log_zero sc))
		match segs with
			| seg :: rest ->
				if i = seg.seg_end then
					let iy' =
						try
							array_index_of fc.transitions.(seg.seg_label) prev_label
						with Not_found -> raise (Zero_potential (tsid,seg,prev_label,"unallowable transition"))

					let miyy' = mi.(seg.seg_label).(iy')
					let id =
						try
							array_index_of (fst miyy') seg.seg_length
						with Not_found -> raise (Zero_potential (tsid,seg,prev_label,"unallowable segment length"))
					let segpot = (snd miyy').(id)
					if is_log_zero segpot then raise (Zero_potential (tsid,seg,prev_label,"feature function returned log_zero"))
					(elnproduct sc segpot,
						seg.seg_label,
						rest)
				else
					sofar
			| [] ->
				assert (i = fc.chain_length)
				let end_label = fc.crf.n_user_labels+1
				try
					let iy' = array_index_of fc.transitions.(end_label) prev_label
					(elnproduct sc (snd mi.(end_label).(iy')).(0),
						end_label,
						[])
				with Not_found -> raise (Zero_potential (tsid,{seg_label = end_label; seg_end = i; seg_length = 1 },prev_label,"unallowable transition into end"))
	(potfolder,(0.0,0,segs))

(* compute the normalization factor Z(x), given the alphas *)
let logz logalphas =
	let (m,n) = Gsl_matrix.dims logalphas
	let tot = ref log_zero
	for j = 0 to n-1 do
		tot := elnsum !tot logalphas.{m-1,j}
	!tot

(* compute the betas. note, beta_0 is not actually represented in the matrix,
	but returned in the vector b0 *)
let compute_logbetas fc =
	let b = Gsl_matrix.create ~init:log_zero (fc.chain_length+1) (fc.crf.n_user_labels+2)
	let b0 = Gsl_vector.create ~init:log_zero (fc.crf.n_user_labels+2)

	(* initialize the last column *)	
	for y = 0 to (fc.crf.n_user_labels) do
		b.{fc.chain_length,y} <- log_zero
	b.{fc.chain_length,fc.crf.n_user_labels+1} <- 0.0 (* log 1.0 *)

	for i' = fc.chain_length-1 downto -1 do
		let (mi':mi) = compute_Mi ~fwd:true fc i'
		for y = 0 to (fc.crf.n_user_labels+1) do
			let mi'y = mi'.(y)
			for iy' = 0 to (Array.length fc.transitions.(y)) - 1 do
				let y' = fc.transitions.(y).(iy')
				let (seg_lengths,seg_potentials) = mi'y.(iy')

				for id = 0 to (Array.length seg_lengths)-1 do
					let d = seg_lengths.(id)
					if i' >= 0 then
						b.{i',y'} <- elnsumprod b.{i'+d,y} seg_potentials.(id) b.{i',y'}
					else
						assert (y' = 0) (* first position can only be start state *)
						b0.{y'} <- elnsumprod b.{i'+d,y} seg_potentials.(id) b0.{y'}
	fc.cache_primed_backward <- true


	(b0,b)

(* compute the total value of each feature for the training sequence (NOT weighted by
	lambda). This is used in SMCRF training *)
let compute_empirical_featvec fc segs = 
	let nff = (Array.length fc.crf.feature_functions)
	let nsff = (Array.length fc.crf.sparse_feature_sets)
	let v = Gsl_vector.create ~init:0.0 (Gsl_vector.length fc.crf.lambda)
	let dbuf = [| 0 |]
	let buf = [| 0. |]
	let rec iter y' = function
		| seg :: rest ->
			let i = seg.seg_end
			let y = seg.seg_label
			let d = seg.seg_length
			assert (i < fc.chain_length)
			assert (d <= fc.crf.label_specs.(y-1).max_length)
			dbuf.(0) <- d

			for j = 0 to (nff-1) do
				segments_values ~no_vec_cache:true fc j y' y i dbuf buf
				v.{j} <- v.{j} +. buf.(0)
				
				
			if nsff > 0 then
				let ssvs = Array.init nsff (fun sj -> sparse_segments_values fc sj y' y i dbuf)
				let ofs = ref nff
				for sj = 0 to nsff-1 do
					let sv = ssvs.(sj).(0)
					SparseVector.add_to_bigarray ~ofs:!ofs sv v
					ofs := !ofs + (SparseVector.length sv)
				assert (!ofs = Gsl_vector.length fc.crf.lambda)

			iter y rest
		| [] ->
			(* add potential for the end position *)
			let i = fc.chain_length
			let y = fc.crf.n_user_labels+1
			
			for j = 0 to (nff-1) do
				segments_values ~no_vec_cache:true fc j y' y i [|1|] buf
				v.{j} <- v.{j} +. buf.(0)
				
			if nsff > 0 then
				let ssvs = Array.init nsff (fun sj -> sparse_segments_values fc sj y' y i [|1|])
				let ofs = ref nff
				for sj = 0 to nsff-1 do
					let sv = ssvs.(sj).(0)
					SparseVector.add_to_bigarray ~ofs:!ofs sv v
					ofs := !ofs + (SparseVector.length sv)
				assert (!ofs = Gsl_vector.length fc.crf.lambda)
				
			v
	iter 0 segs

(* compute the expected value of the feature vector under the current lambda, given X. this
	is the gradient of the partition function.
 *)
let expected_featvec_folder fc logalphas logbetas logZx =
	let nff = (Array.length fc.crf.feature_functions)
	let nsff = (Array.length fc.crf.sparse_feature_sets)
	let inv_logZx = elninv logZx
	let fvs = Expando.make 1 0.
	let exfv_folder i v (mi:mi) =
		let abz = ref log_zero
		let totp = ref log_zero
		for y = 0 to fc.crf.n_user_labels+1 do

			let miy = mi.(y)

			
			assert (abz := elnsumprod logalphas.{i,y} logbetas.{i,y} !abz; true)
			

			for iy' = 0 to (Array.length fc.transitions.(y)) - 1 do
				let y' = fc.transitions.(y).(iy')
				let seg_lengths, seg_potentials = miy.(iy')
				if (Array.length seg_lengths) > 0 then
					let seg_ps =
						Array.init (Array.length seg_lengths)
							fun id ->
								let d = seg_lengths.(id)
								let logalpha_imd_y' =
									if (i-d) = -1 then
										if y' = 0 then 0.0 else log_zero
									else							
										logalphas.{i-d,y'}
								let logp = (elnproduct4 logalpha_imd_y' seg_potentials.(id) logbetas.{i,y} inv_logZx)
								assert (totp := elnsum !totp logp; true)
								eexp logp
					let fvs = Expando.access_ensure fvs (Array.length seg_lengths)
					for j = 0 to nff-1 do
						segments_values fc j y' y i seg_lengths fvs
						for id = 0 to (Array.length seg_lengths)-1 do
							if seg_ps.(id) <> 0. then v.{j} <- v.{j} +. (seg_ps.(id) *. fvs.(id))
					if nsff > 0 then
						let ssvs = Array.init nsff (fun sj -> sparse_segments_values fc sj y' y i seg_lengths)
						for id = 0 to (Array.length seg_lengths)-1 do
							if seg_ps.(id) <> 0. then
								let ofs = ref nff
								for sj = 0 to nsff-1 do
									let sv = ssvs.(sj).(id)
									let psv = SparseVector.map_set (( *.) seg_ps.(id)) sv
									SparseVector.add_to_bigarray ~ofs:!ofs psv v
									ofs := !ofs + (SparseVector.length sv)
								assert (!ofs = Gsl_vector.length fc.crf.lambda)
		if fc.crf.max_max_length = 1 then (* these conditions don't necessarily hold in the semi-Markov case *)
			assert
				let zdiff = logZx -. !abz
				if zdiff *. zdiff < 0.1 then
					true
				else
					eprintf "%d %f %f\n" i logZx !abz
					false
			assert
				let pdiff = 0.0 -. !totp
				if pdiff *. pdiff < 0.1 then
					true
				else
					eprintf "%d %f\n" i !totp
					false
		v
	(exfv_folder,(Gsl_vector.create ~init:0.0 (Gsl_vector.length fc.crf.lambda)))

(******************************************************************************)
(* likelihood functions *)
(******************************************************************************)

(* compute the log-likelihood *)
let log_likelihood_f ?(verbose=false) crf ts fc =
	let (laf,logalpha_zero) = logalpha_folder fc
	let (pf,potential_zero) = potential_folder ts.train_id fc ts.train_segs
	let llfolder i (laf_sofar,pf_sofar) ms = (laf i laf_sofar ms,pf i pf_sofar ms)
	let (logalphas,(potential,_,_)) =
		try 
			fold_Mis fc llfolder (logalpha_zero,potential_zero)
		with
			| Zero_potential (tsid,seg,prev_label,reason) ->
				eprintf  "%s (i=%d,y'=%d,y=%d,d=%d)\n" tsid seg.seg_end prev_label seg.seg_label seg.seg_length
				eprintf "%s\n" reason
				failwith "Training sequence has zero potential due to the above segment; check model design."

	(*
		formally,
			log Z(x) = log(Sum logalpha_L(s))
							s
		we've designed our model with an end state, so we don't have to perform the sum
	*)

	let logZx = logalphas.{ts.train_chain_length,crf.n_user_labels+1}

	if is_log_zero logZx then failwith "SMCRF: partition function evaluated to zero; check model design"

	(* now, log p(y|x) = (log Z(x) + log p(y|x)) - log Z(x) *)
	let lpy = potential -. logZx

	if verbose then
		eprintf "log Psi(Y|X) = %f\n" potential
		eprintf "log Z(X) = %f\n" logZx
		eprintf "log P(Y|X) = %f\n" lpy
		flush stderr

	lpy



(* compute both the log-likelihood and its gradient. (whilst computing the gradient, also computing the log-
	likelihood is essentially free) *)
let log_likelihood_fdf ?(verbose=false) crf ts fc df pv =
	(* we need two passes:
		1. backwards pass to compute the betas and Z(x)
		2. forwards pass to compute the alphas, the likelihood, and the expected feature
			vector (the crucial piece needed for computing the gradient)

		the dynamic programming step for computing E[f_k(X)] only depends on the alphas at previous
		positions. so, we get to compute both the alphas and the expensive part of the gradient in
		one pass. *)
	(* first pass *)
	let (logb0,logbetas) = compute_logbetas fc
	let logZx = logb0.{0}

	if is_log_zero logZx then failwith "SMCRF: partition function evaluated to zero; check model design"

	let (laf,logalpha_zero) = logalpha_folder fc
	let (pf,potential_zero) = potential_folder ts.train_id fc ts.train_segs
	(* we're being tricky here - efvf uses the logalphas before they are fully computed. again,
		at each position, it only needs the alphas at previous positions. *)
	let (efvf,efv_zero) = expected_featvec_folder fc logalpha_zero logbetas logZx

	let metafolder i (laf_sofar,pf_sofar,efvf_sofar) ms =
		let laf_next = laf i laf_sofar ms
		(laf_next,pf i pf_sofar ms,efvf i efvf_sofar ms)

	(* second pass *)
	let (logalphas,(potential,_,_),expected_featvec) =
		try
			fold_Mis fc metafolder (logalpha_zero,potential_zero,efv_zero)
		with
			| Zero_potential (tsid,seg,prev_label,reason) ->
				eprintf  "%s (i=%d,y'=%d,y=%d,d=%d)\n" tsid seg.seg_end prev_label seg.seg_label seg.seg_length
				eprintf "%s\n" reason
				failwith "Training sequence has zero potential due to the above segment; check model design."

	(* sanity check: Z(x) = alpha[L,end] = beta[1,start] *)
	let logZx_diff = abs_float (logZx -. logalphas.{fc.chain_length,fc.crf.n_user_labels+1})
	if logZx_diff /. logZx > 1e-6 then
		failwith (sprintf "CRF.SMCRF: partition function evaluated to different values during forward and backward passes; check forward/backward equivalence of vectorized feature functions. (forward = %e backward = %e)" logalphas.{fc.chain_length,fc.crf.n_user_labels+1} logZx)
		
	let lpy = potential -. logZx
	let empirical_featvec = compute_empirical_featvec fc ts.train_segs
	assert (Gsl_vector.length empirical_featvec = Gsl_vector.length expected_featvec)
	assert (Gsl_vector.length empirical_featvec = Gsl_vector.length fc.crf.lambda)
	for j = 0 to Gsl_vector.length empirical_featvec - 1 do
		df.{j} <- empirical_featvec.{j} -. expected_featvec.{j}
		pv.{j} <- empirical_featvec.{j}

	if verbose then
		eprintf "D(Lambda) = "
		print_vector df
		eprintf "FV(Y|X) = "
		print_vector empirical_featvec
		eprintf "E[FV(Y|X)] = "
		print_vector expected_featvec
		eprintf "log Psi(Y|X) = %f\n" potential
		eprintf "log Z(X) = %f\n" logZx
		eprintf "log P(Y|X) = %f\n" lpy
		flush stderr

	lpy

(******************************************************************************)
(* API *)
(******************************************************************************)

let sparse_feature_set_size = function
	| Emission_features (n,_)
	| Markov_features (n,_)
	| Segment_features_vec (n,_)
	| SemiMarkov_features_vec (n,_) -> n

let make ?lambda ?(seg_length_restrictor=default_seg_length_restrictor) ?(sparse_feature_sets=[||]) ~labels ~feature_functions () =
	let np = (Array.length feature_functions) + (Array.fold_left (+) 0 (Array.map sparse_feature_set_size sparse_feature_sets))

	let lambda =
		match lambda with
			| Some v ->
				if (Array.length v) <> np then
					invalid_arg "SMCRF.make: specified lambda vector must be the same length as the number of feature functions"
				Gsl_vector.of_array v
			| None -> Gsl_vector.create ~init:1.0 np
	let slr = seg_length_restrictor (*match seg_length_restrictor with Some f -> f | None -> default_seg_length_restrictor*)
	let max_max_length = ref 0
	Array.iter (fun l -> max_max_length := max !max_max_length l.max_length) labels

	(* TODO validate input *)

	{ n_user_labels = (Array.length labels);
		label_specs = labels;
		max_max_length = !max_max_length;
		feature_functions = Array.copy feature_functions;
		sparse_feature_sets = Array.copy sparse_feature_sets;
		seg_length_restrictor = slr;
		lambda = lambda }

let get_lambda crf = Gsl_vector.to_array crf.lambda

let set_lambda crf lambda =
	let lambda = Gsl_vector.of_array lambda
	if (Gsl_vector.length lambda) <> (Gsl_vector.length crf.lambda) then invalid_arg "CRF.SMCRF.set_lambda: wrong dimension"
	Gsl_vector.memcpy ~src:lambda ~dst:crf.lambda

type training_configuration = {
	uselambda : bool;
	tol : float;
	epsabs : float;
	random_step : bool;
	restart_b : int;
	restart_m : int;
	maxiter : int;
	regularization_sigma : float;
	verbose : bool
}
let default_training_configuration = {
	uselambda = false;
	tol = 0.1;
	epsabs = 0.1;
	random_step = false;
	restart_b = 0;
	restart_m = 0;
	maxiter = 25;
	regularization_sigma = 0.;
	verbose = false
}

type training_result = {
	converged : bool;			(** true if the numerical step of the training optimization converged in at most [maxiter] steps, false otherwise. If false, the CRF might still be useful: maybe your convergence criteria were too strict. *)
	iterations : int;			(** the number of training optimization iterations performed. *)
	result_f : float;			(** final log-likelihood of the trained model. *)
	result_df : float array		(** final magnitude of the gradient *)
}

let validate_training_sequence crf ts =
	let rec f pos = function
		| seg :: rest ->
			if pos + seg.seg_length <> seg.seg_end then
				invalid_arg (Printf.sprintf "CRF.SMCRF.train: discontiguous segmentation in training sequence %s after position %d" ts.train_id pos)
			if seg.seg_label > 0 && seg.seg_label <= crf.n_user_labels then
				let maxlen = crf.label_specs.(seg.seg_label-1).max_length
				if seg.seg_length > maxlen then
					invalid_arg (Printf.sprintf "CRF.SMCRF.train: training segment length %d exceeds maximum segment length %d for label %d ending at position %d in sequence %s" seg.seg_length maxlen seg.seg_label (pos+seg.seg_length) ts.train_id)
			f seg.seg_end rest
		| [] ->
			if pos+1 <> ts.train_chain_length then
				invalid_arg (Printf.sprintf "CRF.SMCRF.train: training sequence %s has length %d segmentation ends before position %d" ts.train_id ts.train_chain_length (pos+1))
	f (-1) ts.train_segs

type training_sequence_f = (verbose:bool -> lambda:(float array) -> float)*(verbose:bool -> lambda:(float array) -> float*(float array)*(float array))
let prepare_training_sequence ~crf ~training_sequence =
	let ts = training_sequence
	validate_training_sequence crf ts
	let fc = make_featcache crf ts.train_chain_length ts.train_seq
	let f ~verbose ~lambda =
		set_lambda crf lambda
		log_likelihood_f ~verbose:verbose crf ts fc
	let fdf ~verbose ~lambda =
		set_lambda crf lambda
		let df = Gsl_vector.create (Array.length lambda)	
		Gsl_vector.set_zero df
		let pv = Gsl_vector.create (Array.length lambda)	
		Gsl_vector.set_zero pv
		let ll = log_likelihood_fdf ~verbose:verbose crf ts fc df pv
		ll, (Gsl_vector.to_array df), (Gsl_vector.to_array pv)
	(f,fdf)

type progress_function = i:int -> lambda:(float array) -> f:float -> df:(float array) -> pv:(float array) -> [ `Stop | `Go]

let train_advanced ?(configuration=default_training_configuration) ?progress ~crf ~training_sequence_f () =
	
	let np = Gsl_vector.length crf.lambda
	let exn = ref None

(*
	if configuration.verbose then
		eprintf "Allowed transition table:\n"
		Array.iteri
			fun label prev_labels ->
				eprintf "%d <- " label
				Array.iter (eprintf " %d") prev_labels
				eprintf "\n"
				flush stderr
			fcs.(0).transitions
		flush stderr
*)
	
	let total_fdf lambda dodf =
		set_lambda crf lambda

		if configuration.verbose then
			eprintf "-- Evaluate %s at Lambda = " (if dodf then "FDF" else "F")
			print_caml_vector lambda
			flush stderr

		if dodf then
			let f, df, pv = (snd training_sequence_f) ~verbose:configuration.verbose ~lambda:lambda
			if configuration.regularization_sigma <> 0. then
				let sigsq = configuration.regularization_sigma *. configuration.regularization_sigma
				let penalty = (Array.fold_left (fun tot x -> tot +. x *. x) 0. lambda) /. (2. *. sigsq)
				let dpenalty = Array.map (fun l_k -> l_k /. sigsq) lambda
				f -. penalty, (Array.init (Array.length lambda) (fun k -> df.(k) -. dpenalty.(k))), pv
			else
				f, df, pv
		else
			let f = (fst training_sequence_f) ~verbose:configuration.verbose ~lambda:lambda
			if configuration.regularization_sigma <> 0.0 then	
				let sigsq = configuration.regularization_sigma *. configuration.regularization_sigma
				let penalty = (Array.fold_left (fun tot x -> tot +. x *. x) 0. lambda) /. (2. *. sigsq)
				f -. penalty, [||], [||]
			else
				f, [||], [||]

	let minfdf x dodf =
		try
			let f, df, pv = total_fdf x dodf
			match classify_float f with
				| FP_normal | FP_zero -> 0. -. f, (Array.map (fun df_k -> 0. -. df_k) df), pv (* here is where we flip f & df since GSL provides a minimizer *)
				| FP_nan -> failwith "SMCRF.train: training sequence has zero probability; check well-formedness under specified model"
				| sux -> raise (Elog.Numerical_artifact ("minfdf",sux))
		with
			| ex ->
				if configuration.verbose then
					eprintf "%s\n" (Printexc.to_string ex)
					flush stderr
				exn := Some ex
				nan, [||], [||]

	(* for some reason GSL's optimizer sometimes repeatedly evaluates the same point *)
	let memotbl = Hashtbl.create 16
	let memo_minfdf x g dodf =
		let xa = (Gsl_vector.to_array x)
		try
			let (mf,maybe_mg) = Hashtbl.find memotbl xa
			if dodf then
				match maybe_mg with
					| Some (mdf,_) ->
						if configuration.verbose then
							eprintf "-- Recall memoized FDF at Lambda = "
							print_vector x
							eprintf "F = %f DF = " mf
							print_caml_vector mdf
							flush stderr
						Gsl_vector.set_zero g; Gsl_vector.add g (Gsl_vector.of_array mdf)
						mf
					| None ->
						let f, df, pv = minfdf xa true
						Hashtbl.replace memotbl xa (f,Some (df,pv))
						assert
							let fdiff = mf -. f
							fdiff *. fdiff < 0.1
						Gsl_vector.set_zero g; Gsl_vector.add g (Gsl_vector.of_array df)
						f
			else
				if configuration.verbose then
					eprintf "-- Recall memoized F at Lambda = "
					print_vector x
					eprintf "F = %f\n" mf
					flush stderr
				mf
		with
			| Not_found ->
				let f, df, pv = minfdf xa dodf
				Hashtbl.add memotbl xa (f,(if dodf then Some (df,pv) else None))
				if dodf then
					Gsl_vector.set_zero g; Gsl_vector.add g (Gsl_vector.of_array df)
				f

	let fdf = 
			{
				Gsl_fun.multim_f = (fun ~x -> memo_minfdf x (Gsl_vector.create np) false);
				Gsl_fun.multim_df = (fun ~x ~g -> ignore (memo_minfdf x g true));
				Gsl_fun.multim_fdf = (fun ~x ~g -> memo_minfdf x g true)
			}

	let init_lambda = if configuration.uselambda then crf.lambda else Gsl_vector.of_array (Array.init np (fun _ -> Random.float 1.0))
	let make_minimizer lambda = Gsl_multimin.Deriv.make Gsl_multimin.Deriv.VECTOR_BFGS2 np fdf ~x:lambda ~step:(if configuration.random_step then Random.float 1. else 1.) ~tol:configuration.tol
	let minimizer = ref (make_minimizer init_lambda)
	let iteri = ref 1
	let pf_go = ref true
	while (not (Gsl_multimin.Deriv.test_gradient !minimizer configuration.epsabs)) && !iteri <= configuration.maxiter && !exn = None && !pf_go do
		match progress with
			| Some pf ->
				let x = Gsl_vector.create (Gsl_vector.length crf.lambda)
				let f = Gsl_multimin.Deriv.minimum ~x:x !minimizer
				let xa = Gsl_vector.to_array x
				match Hashtbl.find memotbl xa with
					| (_,Some (df,pv)) -> if pf ~i:!iteri ~lambda:xa ~f:(0.0 -. f) ~df:(Array.map ((-.) 0.) df) ~pv:pv = `Stop then pf_go := false
					| _ -> assert false
			| None -> ()
		if (configuration.restart_b > 1 && (!iteri = configuration.restart_b)) || (configuration.restart_m > 1 && !iteri > configuration.restart_b && ((!iteri - configuration.restart_b) mod configuration.restart_m) = 0) then
			if configuration.verbose then
				eprintf "-- Optimizer restart at iteration %d\n" !iteri
				flush stderr
			ignore (Gsl_multimin.Deriv.minimum ~x:crf.lambda !minimizer)
			Gsl_multimin.Deriv.restart !minimizer
			minimizer := make_minimizer crf.lambda
		incr iteri
		Gsl_multimin.Deriv.iterate !minimizer
	match !exn with Some ex -> raise ex | None -> ()
	let g = Gsl_vector.copy crf.lambda
	let ll = Gsl_multimin.Deriv.minimum ~x:crf.lambda !minimizer ~g:g
	Gsl_vector.scale g (-1.0)
	{ converged = (!iteri <= configuration.maxiter);
		iterations = (!iteri - 1);
		result_f = 0.0 -. ll;
		result_df = (Gsl_vector.to_array g) }

let train ?(configuration=default_training_configuration) ?progress ~crf ~training_sequences () =
	let fs = Array.map (fun ts -> prepare_training_sequence ~crf:crf ~training_sequence:ts) training_sequences
	let meta_f ~verbose ~lambda = Array.fold_left (fun ll (f, _) -> ll +. (f ~verbose:verbose ~lambda:lambda)) 0. fs
	let meta_df ~verbose ~lambda =
		let (_, fdf_0) = fs.(0)
		let ll_0, df, pv = fdf_0 ~verbose:verbose ~lambda:lambda
		let ll = ref ll_0
		for i = 1 to (Array.length fs)-1 do
			let (_, fdf) = fs.(i)
			let ll_i, df_i, pv_i = fdf ~verbose:verbose ~lambda:lambda
			ll := !ll +. ll_i
			Array.iteri (fun k df_k -> df.(k) <- df.(k) +. df_k) df_i
			Array.iteri (fun k pv_k -> pv.(k) <- pv.(k) +. pv_k) pv_i
		!ll, df, pv
	let tsf = (meta_f,meta_df)
	match progress with
		| Some pf -> train_advanced ~configuration:configuration ~progress:pf ~crf:crf ~training_sequence_f:tsf ()
		| None -> train_advanced ~configuration:configuration ~crf:crf ~training_sequence_f:tsf ()

let recover_featvec fc y' y i d =
	let nff = Array.length fc.crf.feature_functions
	let nsff = Array.length fc.crf.sparse_feature_sets
	let ds = [| d |]
	let fvs =
		Array.init nff
			fun j ->
				segments_values ~fwd:false fc j y' y i ds box
				box.(0)
	let sfvs = (Array.init nsff (fun sj -> (sparse_segments_values ~fwd:false fc sj y' y i ds).(0)))
	(fvs,sfvs)

let decode_internal score crf seq chain_length =
	let tby' = Bigarray.Array2.create Bigarray.int Bigarray.c_layout (chain_length+1) (crf.n_user_labels+2)
	let tbd = Bigarray.Array2.create Bigarray.int Bigarray.c_layout (chain_length+1) (crf.n_user_labels+2)
	let v = Rotator.make (Array.make_matrix (crf.max_max_length+1) (crf.n_user_labels+2) log_zero)
	let fc = make_featcache crf chain_length seq

	let viterbi_folder i () (mi:mi) =
		for y = 0 to (crf.n_user_labels+1) do
			let miy = mi.(y)
			let cur = ref (log_zero,(-1),(-1))
			assert (Array.length miy = Array.length fc.transitions.(y))

			for iy' = 0 to (Array.length fc.transitions.(y)) - 1 do
				let y' = fc.transitions.(y).(iy')
				let seg_lengths, seg_potentials = miy.(iy')
				
				for id = 0 to (Array.length seg_lengths)-1 do
					let d = seg_lengths.(id)
					assert (d > 0 && (y < 1 || y > crf.n_user_labels || d <= crf.label_specs.(y-1).max_length))
					let (maxscore,_,_) = !cur
					let vdy' = 
						if i-d >= 0 then
							(Rotator.get v d).(y')
						else if y' = 0 then
							0.0 (* log 1.0 *)
						else
							log_zero
					let sc = elnproduct vdy' seg_potentials.(id)
					if (elncompare sc maxscore) > 0 then
						cur := (sc,d,y')

			let (sc,d,y') = !cur

			(Rotator.get v 0).(y) <- sc
			tby'.{i,y} <- y'
			tbd.{i,y} <- d
	
		Rotator.rotate v
	
	fold_Mis fc viterbi_folder ()

	for ylast = 0 to crf.n_user_labels do
		assert (is_log_zero (Rotator.get v 1).(ylast))

	let nff = Array.length fc.crf.feature_functions
	let nsff = Array.length fc.crf.sparse_feature_sets
	let fc = if score then make_featcache crf chain_length seq else fc (* for stupid reasons, we need to restart the feature cache if we're going to request feature values for non-exhausitve segment lengths from it *)
	let rslt = ref []
	let y = ref (crf.n_user_labels+1)
	let i = ref chain_length
	while !i >= 0 do
		assert (!y >= 0)
		let d = tbd.{!i,!y}
		assert (d > 0)
		let y' = tby'.{!i,!y}
		if (!i < chain_length) then
			let sc = ref nan
			if score then
				sc := 0.
				let (dfvs,sfvs) = recover_featvec fc y' !y !i d
				for j = 0 to nff-1 do sc := !sc +. (dfvs.(j) *. fc.crf.lambda.{j})
				let ofs = ref nff
				for sj = 0 to nsff-1 do
					let sv = sfvs.(sj)
					sc := !sc +. (SparseVector.dot_bigarray ~ofs:!ofs sv fc.crf.lambda)
					ofs := !ofs + SparseVector.length sv
				assert (!ofs = Gsl_vector.length fc.crf.lambda)
			rslt := ({ seg_label = !y; seg_end = !i; seg_length = d },!sc) :: !rslt
		y := y'
		i := !i - d
	assert (!i = (-1))
	!rslt, (Rotator.get v 1).(crf.n_user_labels+1)

let decode ~crf ~seq ~chain_length () =
	let lst, psi = decode_internal false crf seq chain_length
	(fst (List.split lst)), psi
let decode_score ~crf ~seq ~chain_length () =
	decode_internal true crf seq chain_length
	

let feature_vectors ~crf ~seq ~chain_length ~segs () =
	let fc = make_featcache crf chain_length seq
	Array.map 
		fun (y',{seg_label = y; seg_end = i; seg_length = d}) ->
			let dfvs, sfvs = recover_featvec fc y' y i d
			Array.concat (dfvs :: (List.map SparseVector.to_array (Array.to_list sfvs)))
		segs

let posterior_decode ~crf ~seq ~chain_length () =
	let fc = make_featcache crf chain_length seq
	let (b0,logbetas) = compute_logbetas fc
	let logZx = b0.{0}
	let (laf,logalpha_zero) = logalpha_folder fc
	let logalphas = fold_Mis fc laf logalpha_zero
	let pd = Gsl_matrix.create (chain_length+1) (crf.n_user_labels+2)
	for i = 0 to chain_length do
		for y = 0 to crf.n_user_labels+1 do
			pd.{i,y} <- eexp (elnproduct3 logalphas.{i,y} logbetas.{i,y} (0.0 -. logZx))
	pd,logZx

(*
	bit of code that may be useful for posterior decoding in the semi-Markov case:

	let logn_folder fc la lb =
		let n =
			if fc.crf.max_max_length > 1 then
				Gsl_matrix.create ~init:log_zero (fc.chain_length+1) (fc.crf.n_user_labels+2)
			else
				Gsl_matrix.create 0 0
		let folder im1 n (mim1:mi) =
			if fc.crf.max_max_length > 1 then
				let i = im1+1
				if i >= 2 then
					let mfim2 = compute_Mi ~fwd:true fc (i-2)
					for y = 0 to (fc.crf.n_user_labels+1) do
						if fc.crf.label_specs.(y).max_length >= 2 then
							n.{i,y} <- n.{i-1,y}
							for iy' = 0 to (Array.length fc.transitions.(y))-1 do
								let y' = fc.transitions.(y).(iy')
							
								let myla = la.{i-2,y'}
								let (fwd_lens,fwd_pots) = mfim2.(y).(iy')
								for id = 0 to (Array.length fwd_lens)-1 do
									let d = fwd_lens.(id)
									if d >= 2 then
										n.{i,y} <- elnsum n.{i,y} (elnproduct3 myla fwd_pots.(id) lb.{i+d-2,y})
							
								let mylb = lb.{i-1,y}
								let (bwd_lens,bwd_pots) = mim1.(y).(iy')
								for id = 0 to (Array.length bwd_lens)-1 do
									let d = bwd_lens.(id)
									if d >= 2 then
										n.{i,y} <- elndiff n.{i,y} (elnproduct3 la.{i-d-1,y'} bwd_pots.(id) mylb)
			n
		(folder,n)
*)

module MESS = struct
	let (+=.) rf x = rf := !rf +. x
	type 'obs_seq score_function =
		| Markov_score of (seq:'obs_seq training_sequence -> prev_label:int -> label:int -> i:int -> chain_length:int -> float)
		| SemiMarkov_score of (seq:'obs_seq training_sequence -> prev_label:int -> label:int -> i:int -> seg_length:int -> chain_length:int -> float)
		| SemiMarkov_score_vec of (seq:'obs_seq training_sequence -> prev_label:int -> label:int -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> values:(float array) -> unit)

	let segments_scores fc ts fwd y' y i lens scores = function
		| SemiMarkov_score_vec f -> f ~seq:ts ~prev_label:y' ~label:y ~i:i ~seg_lengths:lens ~chain_length:(fc.chain_length) ~fwd:fwd ~values:scores
		| SemiMarkov_score f -> segments_semimarkov_values fwd (fun i d -> f ~seq:ts ~prev_label:y' ~label:y ~i:i ~seg_length:d ~chain_length:fc.chain_length) i lens scores
		| Markov_score f -> segments_markov_sums fwd (fun first i -> f ~seq:ts ~prev_label:(if first then y' else y) ~label:y ~i:i ~chain_length:fc.chain_length) i lens scores

	let epsilon = 1e-3

	let gamma_folder fc ts la b0 lb inv_logZx s =
		let scs = Expando.make 1 0.
		let folder i g (mi:mi) =
			for y = 0 to (fc.crf.n_user_labels+1) do
				for iy' = 0 to (Array.length fc.transitions.(y))-1 do
					let y' = fc.transitions.(y).(iy')
					let lens, pots = mi.(y).(iy')
					let scs = Expando.access_ensure scs (Array.length lens)
					segments_scores fc ts false y' y i lens scs s

					for id = 0 to (Array.length lens)-1 do
						let d = lens.(id)
						let pot = pots.(id)
						let myla = if i-d>=0 then la.{i-d,y'} else if y'=0 then 0. else log_zero
						let lpr = elnproduct4 myla pot lb.{i,y} inv_logZx
						if not (is_log_zero lpr) then
							assert (lpr <= epsilon)
							let mylb = if i-d>=0 then lb.{i-d,y'} else b0.{y'}
							let clpr = elnproduct lpr (elninv (elnproduct3 myla mylb inv_logZx))
							assert ((not (is_log_zero clpr)) && clpr <= epsilon && lpr -. clpr <= epsilon)
							g.{i,y} <- g.{i,y} +. (eexp lpr) *. scs.(id) +. (eexp clpr) *. (if i-d>=0 then g.{i-d,y'} else 0.)
			g
		folder, (Gsl_matrix.create ~init:0. (fc.chain_length+1) (fc.crf.n_user_labels+2))

	let compute_delta_i' fc i' mi' ts la lb inv_logZx s drot =
		let di' = Rotator.get drot 0
		let scs = Expando.make 1 0.
		for y = 0 to (fc.crf.n_user_labels+1) do di'.(y) <- 0.
		for y = 0 to (fc.crf.n_user_labels+1) do
			for iy' = 0 to (Array.length fc.transitions.(y)) - 1 do
				let y' = fc.transitions.(y).(iy')
				let lens, pots = mi'.(y).(iy')
				let scs = Expando.access_ensure scs (Array.length lens)
				segments_scores fc ts true y' y i' lens scs s

				for id = 0 to (Array.length lens)-1 do
					let d = lens.(id)
					let myla = if i'>=0 then la.{i',y'} else if y'=0 then 0. else log_zero
					let lpr = elnproduct4 myla pots.(id) lb.{i'+d,y} inv_logZx
					if not (is_log_zero lpr) then
						assert (lpr <= epsilon)
						let clpr = elnproduct lpr (elninv (elnproduct3 la.{i'+d,y} lb.{i'+d,y} inv_logZx))
						assert ((not (is_log_zero clpr)) && clpr <= epsilon && lpr -. clpr <= epsilon)
						di'.(y') <- di'.(y') +. (eexp lpr) *. scs.(id) +. (eexp clpr) *. (Rotator.get drot d).(y)

	let esc_folder fc ts la lb inv_logZx s =
		let scs = Expando.make 1 0.
		let folder i sofar (mi:mi) =
			let sci = ref 0.
			for y = 0 to fc.crf.n_user_labels+1 do	
				for iy' = 0 to (Array.length fc.transitions.(y)) - 1 do
					let y' = fc.transitions.(y).(iy')
					let lens, pots = mi.(y).(iy')
					let scs = Expando.access_ensure scs (Array.length lens)
					segments_scores fc ts false y' y i lens scs s
					for id = 0 to (Array.length lens)-1 do
						let d = lens.(id)
						let myla = if i-d>=0 then la.{i-d,y'} else if y'=0 then 0. else log_zero
						let lpr = elnproduct4 myla pots.(id) lb.{i,y} inv_logZx
						if not (is_log_zero lpr) then
							assert (lpr <= epsilon)
							sci +=. (eexp lpr) *. scs.(id)
			sofar +. !sci
		folder, 0.

	let compute_expected_scorefeatvec fc ts la b0 lb inv_logZx g s =
		let nff = Array.length fc.crf.feature_functions
		let nsff = (Array.length fc.crf.sparse_feature_sets)
		let esc = ref 0.
		let efv = Gsl_vector.create ~init:0. (Gsl_vector.length fc.crf.lambda)
		let escfv = Gsl_vector.create ~init:0. (Gsl_vector.length fc.crf.lambda)
		let drot = Rotator.make (Array.make_matrix (fc.crf.max_max_length+2) (fc.crf.n_user_labels+2) 0.)
		let scs = Expando.make 1 0.
		let lprs = Expando.make 1 0.
		let esparts = Expando.make 1 0.
		let fvs = Expando.make 1 0.

		for i' = fc.chain_length-1 downto -1 do
			let (mi':mi) = compute_Mi ~fwd:true fc i'
			compute_delta_i' fc i' mi' ts la lb inv_logZx s drot
			
			for y = 0 to fc.crf.n_user_labels+1 do	
				for iy' = 0 to (Array.length fc.transitions.(y)) - 1 do
					let y' = fc.transitions.(y).(iy')
					let lens, pots = mi'.(y).(iy')
					let nds = Array.length lens

					let prs = Expando.access_ensure lprs nds
					let esparts = Expando.access_ensure esparts nds
					let scs = Expando.access_ensure scs nds
					let any = ref false

					segments_scores fc ts true y' y i' lens scs s

					for id = 0 to nds-1 do
						let d = lens.(id)
						let myla = if i'>=0 then la.{i',y'} else if y'=0 then 0. else log_zero
						let lpr = elnproduct4 myla pots.(id) lb.{i'+d,y} inv_logZx
						
						if is_log_zero lpr then
							prs.(id) <- 0.
						else
							assert (lpr <= epsilon)
							any := true
							let pr = eexp lpr
							prs.(id) <- pr
							let sc = scs.(id)
							esc +=. pr *. sc

							let gammapart =
								if i'<0 then
									0.
								else
									let lpr_cy' = elnproduct lpr (elninv (elnproduct3 myla lb.{i',y'} inv_logZx))
									assert ((not (is_log_zero lpr_cy')) && lpr_cy' <= epsilon && lpr -. lpr_cy' <= epsilon)
									(eexp lpr_cy') *. g.{i',y'}

							let lpr_cy = elnproduct lpr (elninv (elnproduct3 la.{i'+d,y} lb.{i'+d,y} inv_logZx))
							assert ((not (is_log_zero lpr_cy)) && lpr_cy <= epsilon && lpr -. lpr_cy <= epsilon)
							let deltapart =	(eexp lpr_cy) *. (Rotator.get drot d).(y)

							esparts.(id) <- gammapart +. deltapart +. pr *. sc

					if !any then
						let fvs = Expando.access_ensure fvs nds
						for k = 0 to nff-1 do
							segments_values ~fwd:true fc k y' y i' lens fvs
							for id = 0 to nds-1 do
								if prs.(id) > 0. then
									efv.{k} <- efv.{k} +. prs.(id) *. fvs.(id)
									escfv.{k} <- escfv.{k} +. fvs.(id) *. esparts.(id)
									
						if nsff > 0 then
							let ssvs = Array.init nsff (fun sj -> sparse_segments_values ~fwd:true fc sj y' y i' lens)
							for id = 0 to nds-1 do
								if prs.(id) > 0. then
									let ofs = ref nff
									for sj = 0 to nsff-1 do
										let sv = ssvs.(sj).(id)
										let psv = SparseVector.map_set (( *.) prs.(id)) sv
										SparseVector.add_to_bigarray ~ofs:!ofs psv efv
										let escfsv = SparseVector.map_set (( *.) esparts.(id)) sv
										SparseVector.add_to_bigarray ~ofs:!ofs escfsv escfv
										ofs := !ofs + (SparseVector.length sv)
									assert (!ofs = Gsl_vector.length fc.crf.lambda)

			Rotator.rotate drot

		!esc, efv, escfv, (Rotator.get drot 1).(0)

	let objective_function_f ?(verbose=false) fc ts s =
		let (b0,lb) = compute_logbetas fc
		let inv_logZx = elninv b0.{0}

		let (laf,la) = logalpha_folder fc
		let (escf,esc_zero) = esc_folder fc ts la lb inv_logZx s

		let metafolder i (la_sofar,esc_sofar) mi =
			let la_next = laf i la_sofar mi
			(la_next,escf i esc_sofar mi)

		let (la,esc) = fold_Mis fc metafolder (la,esc_zero)

		if verbose then
			eprintf "log Z(X) = %f\n" b0.{0}
			eprintf "E[S(Y,Y*,X)] = %f\n" esc
			flush stderr

		esc

	let objective_function_fdf ?(verbose=false) fc ts s =
		(* 1. backward pass computes betas *)
		let (b0,lb) = compute_logbetas fc
		let inv_logZx = elninv b0.{0}

		(* 2. forward pass computes alphas and gammas *)
		let (laf,la) = logalpha_folder fc
		let (gf,g) = gamma_folder fc ts la b0 lb inv_logZx s

		let metafolder i (la_sofar,g_sofar) mi =
			let la_next = laf i la_sofar mi
			(la_next,gf i g_sofar mi)

		let (la,g) = fold_Mis fc metafolder (la,g)
		
		(* 3. backward pass to compute E[FV*SC], E[SC], and E[FV] (and deltas on-the-fly) *)

		let esc, efv, escfv, d0start = compute_expected_scorefeatvec fc ts la b0 lb inv_logZx g s
		
		let df = Array.init (Gsl_vector.length fc.crf.lambda) (fun k -> escfv.{k} -. esc *. efv.{k})
		let pv = Gsl_vector.to_array escfv

		if verbose then
			eprintf "log Z(X) = %f\n" b0.{0}
			eprintf "gamma(L,End) = %f\n" g.{fc.chain_length,fc.crf.n_user_labels+1}
			eprintf "delta(-1,Start) = %f\n" d0start
			eprintf "E[S(Y,Y*,X)FV(Y|X)] = "
			print_vector escfv
			eprintf "E[FV(Y|X)] = "
			print_vector efv
			eprintf "E[S(Y,Y*,X)] = %f\n" esc
			eprintf "DF = "
			print_caml_vector df
			flush stderr

		esc, df, pv

	let prepare_training_sequence ~crf ~training_sequence ~score_function =
		let ts = training_sequence
		validate_training_sequence crf ts
		let fc = make_featcache crf ts.train_chain_length ts.train_seq
		let f ~verbose ~lambda =
			set_lambda crf lambda
			objective_function_f ~verbose:verbose fc ts score_function
		let fdf ~verbose ~lambda =
			set_lambda crf lambda
			objective_function_fdf ~verbose:verbose fc ts score_function
		(f,fdf)
