module vqubit_native

// Complex numbers cross the Python boundary as a flat, interleaved array of
// `f64`: `[re0, im0, re1, im1, ...]`. V has a complex type in the standard
// library, but the marshalling layer knows real numbers, so a dense buffer of
// them is what both sides agree on and what keeps every call copy-free enough
// to be fast.
//
// A matrix of dimension `d` is row-major and `2 * d * d` floats long. A state
// vector of `n` qubits has `2**n` amplitudes and is `2 * 2**n` floats long.
// Qubits are ordered big-endian, as in Cirq: qubit 0 is the most significant
// bit of a basis-state index.

// matmul returns the product `a @ b` of two `dim`-by-`dim` complex matrices.
@[vc_fn]
pub fn matmul(a []f64, b []f64, dim int) []f64 {
	if dim <= 0 {
		return []f64{}
	}
	mut out := []f64{len: 2 * dim * dim, init: 0.0}
	for i in 0 .. dim {
		for k in 0 .. dim {
			ar := a[2 * (i * dim + k)]
			ai := a[2 * (i * dim + k) + 1]
			if ar == 0.0 && ai == 0.0 {
				continue
			}
			for j in 0 .. dim {
				br := b[2 * (k * dim + j)]
				bi := b[2 * (k * dim + j) + 1]
				out[2 * (i * dim + j)] += ar * br - ai * bi
				out[2 * (i * dim + j) + 1] += ar * bi + ai * br
			}
		}
	}
	return out
}

// kron returns the Kronecker (tensor) product of two complex matrices.
//
// `a` has dimension `a_dim` and `b` has dimension `b_dim`; the result has
// dimension `a_dim * b_dim`.
@[vc_fn]
pub fn kron(a []f64, a_dim int, b []f64, b_dim int) []f64 {
	if a_dim <= 0 || b_dim <= 0 {
		return []f64{}
	}
	d := a_dim * b_dim
	mut out := []f64{len: 2 * d * d, init: 0.0}
	for ia in 0 .. a_dim {
		for ja in 0 .. a_dim {
			ar := a[2 * (ia * a_dim + ja)]
			ai := a[2 * (ia * a_dim + ja) + 1]
			for ib in 0 .. b_dim {
				for jb in 0 .. b_dim {
					i := ia * b_dim + ib
					j := ja * b_dim + jb
					br := b[2 * (ib * b_dim + jb)]
					bi := b[2 * (ib * b_dim + jb) + 1]
					out[2 * (i * d + j)] = ar * br - ai * bi
					out[2 * (i * d + j) + 1] = ar * bi + ai * br
				}
			}
		}
	}
	return out
}

// apply_gate applies a complex `gate` to a state vector of `num_qubits` qubits.
//
// The gate has dimension `2**targets.len` and acts on the qubits listed in
// `targets`, with the first target the most significant. `state` is a flat,
// interleaved array of `2**num_qubits` amplitudes and the result has the same
// shape.
@[vc_fn]
pub fn apply_gate(state []f64, gate []f64, num_qubits int, targets []int) []f64 {
	if num_qubits <= 0 || targets.len == 0 {
		return state.clone()
	}
	n := 1 << num_qubits
	g := 1 << targets.len
	mut out := []f64{len: 2 * n, init: 0.0}
	mut tmask := 0
	for t in targets {
		tmask |= 1 << (num_qubits - 1 - t)
	}
	for i in 0 .. n {
		mut sub := 0
		for t in targets {
			bit := (i >> (num_qubits - 1 - t)) & 1
			sub = (sub << 1) | bit
		}
		base := i & ~tmask
		sr := state[2 * i]
		si := state[2 * i + 1]
		for row in 0 .. g {
			mut new_i := base
			for j, t in targets {
				bit := (row >> (targets.len - 1 - j)) & 1
				new_i |= bit << (num_qubits - 1 - t)
			}
			gr := gate[2 * (row * g + sub)]
			gi := gate[2 * (row * g + sub) + 1]
			out[2 * new_i] += gr * sr - gi * si
			out[2 * new_i + 1] += gr * si + gi * sr
		}
	}
	return out
}

// dagger returns the conjugate transpose of a `dim`-by-`dim` complex matrix.
@[vc_fn]
pub fn dagger(a []f64, dim int) []f64 {
	if dim <= 0 {
		return []f64{}
	}
	mut out := []f64{len: 2 * dim * dim, init: 0.0}
	for i in 0 .. dim {
		for j in 0 .. dim {
			out[2 * (i * dim + j)] = a[2 * (j * dim + i)]
			out[2 * (i * dim + j) + 1] = -a[2 * (j * dim + i) + 1]
		}
	}
	return out
}

// add_scaled returns `a + scale * b` for two complex matrices of dimension
// `dim`, where `scale` is a real (not complex) multiplier. It is the building
// block for sums of unitaries that share a shape.
@[vc_fn]
pub fn add_scaled(a []f64, b []f64, dim int, scale f64) []f64 {
	if dim <= 0 {
		return []f64{}
	}
	mut out := []f64{len: 2 * dim * dim, init: 0.0}
	for k in 0 .. 2 * dim * dim {
		out[k] = a[k] + scale * b[k]
	}
	return out
}
