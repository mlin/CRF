(*pp ocaml+twt *)

(** Chain-structured, semi-Markov conditional random fields. *)

(** {1 Creating a CRF} *)

(** labels (hidden rv assignments) are specified by number. label 0 is always start. The highest valid label is always end. The start and end labels are inherent to the model: if you specify to create a CRF with two labels, the CRF will actually have four: 0=start, 1=yourstate1, 2=yourstate2, and 3=end. *)
type label = int

(** "Settings" for each user-specified label. The model can be made purely Markov by setting the [max_length] of all labels to 1. *)
type label_spec = {
	(** The maximum semi-Markov segment length for this label. Must be at least 1. *)
	max_length : int
}

type segment = {
	seg_label : label;
	seg_end : int;
	seg_length : int
}

(** A segmentation is a list of segments where [s(i).seg_end = s(i+1).seg_end - s(i+1).seg_length] and [seg_length] is at least 1. The segmentation does not include the start or end positions. *)
type segmentation = segment list

(** The type of feature functions. The algorithm can use the information that a feature function is independent of certain parameters to cache values, reducing the number of times it has to call that feature function. Feature functions may return any finite real value, or [nan] to indicate that the given combination of inputs is invalid (leads to a zero-probability path). *)
type 'obs_seq feature_function =
	| SemiMarkov_feature of		(prev_label:int ->	label:int ->	seq:'obs_seq ->	i:int ->	seg_length:int ->	chain_length:int ->	float)*[`NoCache]
	| Segment_feature of		(					label:int ->	seq:'obs_seq ->	i:int ->	seg_length:int ->	chain_length:int ->	float)*[`NoCache]
	| Markov_feature of			(prev_label:int ->	label:int ->	seq:'obs_seq -> i:int ->						chain_length:int ->	float)*[`NoCache|`DenseCache|`SparseCache]
	| Emission_feature of		(					label:int ->	seq:'obs_seq -> i:int -> 						chain_length:int -> float)*[`NoCache|`DenseCache|`SparseCache]
	| Length_feature of			(					label:int ->	seq:'obs_seq ->				seg_length:int	->	chain_length:int -> float)
	| Transition_feature of		(prev_label:int ->	label:int ->	seq:'obs_seq ->									chain_length:int -> float)

	| SemiMarkov_feature_vec of	(prev_label:int -> label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> values:(float array) -> unit)*[`NoCache|`DenseCache|`SparseCache]
	| Segment_feature_vec of	(label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> values:(float array) -> unit)*[`NoCache|`DenseCache|`SparseCache]

type 'obs_seq sparse_feature_set =
	| Emission_features of int*(label:int -> seq:'obs_seq -> i:int -> chain_length:int -> SparseVector.t)
	| Markov_features of int*(prev_label:int -> label:int -> seq:'obs_seq -> i:int -> chain_length:int -> SparseVector.t)
	| Segment_features_vec of int*(label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> SparseVector.t array)
	| SemiMarkov_features_vec of int*(prev_label:int -> label:int -> seq:'obs_seq -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> SparseVector.t array)

(** A SMCRF can optionally have a "segment length restrictor" that causes the algorithms to only consider certain segment lengths at any point (that is, where they would normally iterate over segment lengths from one to the maximum, they can be restricted to look at a subset of those lengths). This can tremendously speed it up.

If [fwd=false], the restrictor returns an array of all allowable lengths for a segment ending at position [i] inclusive, with label [label] and previous segment label [prev_label].

If [fwd=true], the restrictor returns an array of all allowable lengths for a segment beginning at position [i+1] inclusive, with label [label] and previous segment label [prev_label].

In either case, the returned array must be sorted in increasing order and each element must be between 1 and [max_seg_length], inclusive.

The segment length restrictor must compute identical sets of segments on forward and backward passes. That is, the union of all segments computed when [fwd=false] (at all combinations of [i], [label], and [prev_label]) must equal the corresponding union of all segments computed when [fwd=true], even though it must return different subsets at each position depending on [fwd]. The MCL training routine will raise an exception if it detects that this forward-backward equivalence is violated.
 *)
type 'obs_seq seg_length_restrictor = (prev_label:int -> label:int -> seq:'obs_seq -> i:int -> max_seg_length:int -> chain_length:int -> fwd:bool -> int array)

(** The abstract type of a CRF. *)
type 'obs_seq t

(** Create a CRF. The weight vector is initialized to [1.0, 1.0, 1.0, ...] if you do not specify one. The default segment length restrictor considers all possible segment lengths (up to the specified maximums). *)
val make : ?lambda:(float array) -> ?seg_length_restrictor:('obs_seq seg_length_restrictor) -> ?sparse_feature_sets:('obs_seq sparse_feature_set array) -> labels:(label_spec array) -> feature_functions:('obs_seq feature_function array) -> unit -> 'obs_seq t

(** Retrieve a copy of the weight vector. The lambda vector can be serialized, and a trained CRF can be resurrected by providing it to [make] or [set_lambda]. *)
val get_lambda : ('obs_seq t) -> float array

(** Replace the weight vector of a CRF. *)
val set_lambda : ('obs_seq t) -> float array -> unit

(** {1 Training your CRF (simple interface for conditional maximum likelihood training)} *)

(** Settings controlling various detailed aspects of training. *)
type training_configuration = {
	uselambda : bool;	(** (default false) If true, use the lambda vector already in the CRF as the starting point. Otherwise, use a random starting point. *)
	tol : float;		(** (default 0.1) Tolerance parameter affecting line search accuracy at each step. *)
	epsabs : float;		(** (default 0.1) Convergence criterion. The numerical optimization stops when the magnitude of the gradient at the present maximum is less than [epsabs]. *)
	random_step : bool; (** (default false) If true, use a random initial step size. *)
	restart_b : int;	(** (default 0) If greater than 1, clear the BFGS state once after [restart_b] iterations *)
	restart_m : int;	(** (default 0) If greater than 1, clear the BFGS state every [restart_m] iterations (after [restart_b] iterations) *)
	maxiter : int;		(** (default 25) Maximum numerical optimization iterations to perform. *)
	regularization_sigma : float; (** (default 0.0) If nonzero, applies regularization to control overfitting. *)
	verbose : bool;		(** (default false) print verbose debugging information on each iteration *)
}

val default_training_configuration : training_configuration

(** The type of a training sequence. *)
type 'obs_seq training_sequence = {
	train_seq : 'obs_seq;			(** observation sequence *)
	train_chain_length : int;		(** length of the CRF chain (usually the same as the length of the observation sequence) *)
	train_segs : segmentation;		(** the known segmentation for the chain. *)
	train_id : string				(** identifying string by which this sequence will be referred in error messages *)
}

type training_result = {
	converged : bool;			(** true if the numerical step of the training optimization converged in at most [maxiter] steps, false otherwise. If false, the CRF might still be useful: maybe your convergence criteria were too strict. *)
	iterations : int;			(** the number of training optimization iterations performed. *)
	result_f : float;			(** final value of the objective function (normally log-likelihood) *)
	result_df : float array		(** gradient of the objection function at the final lambda *)
}

(** Zero_potential may be raised during training if the model assigns zero potential to a training sequence. In this case, the model design should be double-checked. *)
exception Zero_potential of string*segment*int*string

type progress_function = i:int -> lambda:(float array) -> f:float -> df:(float array) -> pv:(float array) -> [ `Stop | `Go ]

(** Train a CRF using the given training data.

@return true if the optimization converged to the specified criteria in at most [maxiter] steps, false otherwise. *)
val train : ?configuration:training_configuration -> ?progress:progress_function -> crf:('obs_seq t) -> training_sequences:('obs_seq training_sequence array) -> unit -> training_result

(** {1 Training (advanced interface for parallelization and different objective functions)} *)

(** The type of an objective function *)
type training_sequence_f = (verbose:bool -> lambda:(float array) -> float)*(verbose:bool -> lambda:(float array) -> float*(float array)*(float array))

(** Create the CML objective function for the given training sequence *)
val prepare_training_sequence : crf:('obs_seq t) -> training_sequence:('obs_seq training_sequence) -> training_sequence_f

(** Train using the given objective function *)
val train_advanced : ?configuration:training_configuration -> ?progress:progress_function -> crf:('obs_seq t) -> training_sequence_f:training_sequence_f -> unit -> training_result

(** {1 Inference} *)

(** Find the highest-potential (most probable) segmentation.
	@return the segmentation and its total log-potential *)
val decode : crf:('obs_seq t) -> seq:'obs_seq -> chain_length:int -> unit -> segmentation*float

(** Decode, and also recover the log-potential of each individual segment in the resulting segmentation. This requires a little bit of extra overhead. *)
val decode_score : crf:('obs_seq t) -> seq:'obs_seq -> chain_length:int -> unit -> ((segment*float) list)*float


(** Find the posterior decoding for a chain with the specified length, given an input sequence.
	@return a matrix P such that [P(i,y)] = the posterior probability that the chain has label [y] at position [i].
*)
val posterior_decode : crf:('obs_seq t) -> seq:'obs_seq -> chain_length:int -> unit -> Gsl_matrix.matrix*float

(** Recover the feature vectors for the given segments. For each segment, the label of the preceding segment must also be specified.
*)
val feature_vectors : crf:('obs_seq t) -> seq:'obs_seq -> chain_length:int -> segs:((int*segment) array) -> unit -> float array array

module MESS : sig
	type 'obs_seq score_function =
		| Markov_score of (seq:'obs_seq training_sequence -> prev_label:int -> label:int -> i:int -> chain_length:int -> float)
		| SemiMarkov_score of (seq:'obs_seq training_sequence -> prev_label:int -> label:int -> i:int -> seg_length:int -> chain_length:int -> float)
		| SemiMarkov_score_vec of (seq:'obs_seq training_sequence -> prev_label:int -> label:int -> i:int -> seg_lengths:(int array) -> chain_length:int -> fwd:bool -> values:(float array) -> unit)

	val prepare_training_sequence : crf:('obs_seq t) -> training_sequence:('obs_seq training_sequence) -> score_function:('obs_seq score_function) -> training_sequence_f

