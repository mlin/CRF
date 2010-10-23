(*pp ocaml+twt *)

(** purely-applicative representation of sparse numerical vectors *)

type t

(** [make len] creates a sparse vector of length [len] with all elements initialized to zero, in [O(1)] time *)
val make : int -> t

(** return [L], the length of the vector, in O(1) time *)
val length : t -> int

(** retrieve an element from the vector in [O(lg N)] time where [N] is the number of explicitly-set elements in the sparse vector *)
val get : t -> int -> float

(** return a vector with an element replaced, in amortized [O(lg N)] time *)
val set : t -> int -> float -> t

(** grow a vector in [O(1)] time. no effect if the vector is already at least as long as requested *)
val grow : t -> int -> t

(** garbage-collect any elements that have been explicitly set to zero, in [O(N lg N)] time. they'll still be zero, but won't take up memory.
@param threshold any element whose absolute value is less than or equal to [threshold] is reset to zero, and memory requirements correspondingly reduced. default 0
*)
val clean : ?threshold:float -> t -> t

(** [sub v lo hi] extracts part of the vector in [O(N)] time *)
val sub : t -> int -> int -> t

(** create a sparse vector from a float array. any zeroes in the array will not take up memory in the sparse vector.*)
val of_array : float array -> t

(** create a float array from a sparse vector, in [O(L+N)] time *)
val to_array : t -> float array

(** return [N] in [O(N)] time *)
val count_set : t -> int

(** iterate only over the elements in the vector that have been explicitly set. (note: elements can be explicitly set to zero, so this is not necessarily only nonzero elements) *)
val iter_set : (float -> unit) -> t -> unit
val map_set : (float -> float) -> t -> t
val fold_set : (float -> 'a -> 'a) -> t -> 'a -> 'a


val iteri_set : (int -> float -> unit) -> t -> unit
val mapi_set : (int -> float -> float) -> t -> t
val foldi_set : (int -> float -> 'a -> 'a) -> t -> 'a -> 'a

(** calculate the norm of the vector
@param l if [l=2], calculates the l2-norm; if [l=3], calculates the l3-norm; etc. default 2 *)
val norm : ?l:int -> t -> float

(** add two vectors of the same length in [O(N lg N)] time (where [N] is the larger of the two) *)
val add : t -> t -> t

(** calculate the dot product of two vectors of the same length, in [O(N lg N)] time. The behavior is undefined if any element of either vector is not finite! *)
val dot : t -> t -> float

(** calculate the dot product of a vector and a float array of the same length, in [O(N)] time. The behavior is undefined if any element of either vector is not finite! *)
val dot_array : t -> float array -> float

(** add the elements of the sparse vector to the corresponding elements in a float array of the same length or longer.
@param ofs if specified, entry [k] of the sparse vector is added to entry [k+ofs] of the float array. The float array must have length at least [L+ofs].
*)
val add_to_array : ?ofs:int -> t -> float array -> unit

type double_bigarray = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
val dot_bigarray : ?ofs:int -> t -> double_bigarray -> float
val add_to_bigarray : ?ofs:int -> t -> double_bigarray -> unit

(** append two sparse vectors in [O(N2 lg (N1+N2))] time *)
val append : t -> t -> t
val concat : t list -> t
