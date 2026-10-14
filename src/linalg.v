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

import math

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

// matvec returns the product `matrix @ vector` for a `dim`-by-`dim` complex
// matrix and a complex vector of length `dim`.
@[vc_fn]
pub fn matvec(matrix []f64, vector []f64, dim int) []f64 {
	if dim <= 0 {
		return []f64{}
	}
	mut out := []f64{len: 2 * dim, init: 0.0}
	for i in 0 .. dim {
		mut re := 0.0
		mut im := 0.0
		for j in 0 .. dim {
			ar := matrix[2 * (i * dim + j)]
			ai := matrix[2 * (i * dim + j) + 1]
			xr := vector[2 * j]
			xi := vector[2 * j + 1]
			re += ar * xr - ai * xi
			im += ar * xi + ai * xr
		}
		out[2 * i] = re
		out[2 * i + 1] = im
	}
	return out
}

// sandwich returns `u @ a @ u^dagger`, the update of a density matrix by a
// unitary, for `dim`-by-`dim` matrices.
@[vc_fn]
pub fn sandwich(a []f64, u []f64, dim int) []f64 {
	if dim <= 0 {
		return []f64{}
	}
	size := 2 * dim * dim
	mut tmp := []f64{len: size, init: 0.0}
	for i in 0 .. dim {
		for k in 0 .. dim {
			ur := u[2 * (i * dim + k)]
			ui := u[2 * (i * dim + k) + 1]
			if ur == 0.0 && ui == 0.0 {
				continue
			}
			for j in 0 .. dim {
				ar := a[2 * (k * dim + j)]
				ai := a[2 * (k * dim + j) + 1]
				tmp[2 * (i * dim + j)] += ur * ar - ui * ai
				tmp[2 * (i * dim + j) + 1] += ur * ai + ui * ar
			}
		}
	}
	mut out := []f64{len: size, init: 0.0}
	for i in 0 .. dim {
		for k in 0 .. dim {
			tr := tmp[2 * (i * dim + k)]
			ti := tmp[2 * (i * dim + k) + 1]
			if tr == 0.0 && ti == 0.0 {
				continue
			}
			for j in 0 .. dim {
				ur := u[2 * (j * dim + k)]
				ui := -u[2 * (j * dim + k) + 1]
				out[2 * (i * dim + j)] += tr * ur - ti * ui
				out[2 * (i * dim + j) + 1] += tr * ui + ti * ur
			}
		}
	}
	return out
}

// channel returns `sum_i K_i @ rho @ K_i^dagger` for a list of `num_ops`
// `dim`-by-`dim` Kraus operators concatenated in `ops`.
@[vc_fn]
pub fn channel(rho []f64, ops []f64, num_ops int, dim int) []f64 {
	size := 2 * dim * dim
	mut out := []f64{len: size, init: 0.0}
	for o in 0 .. num_ops {
		k := ops[o * size..(o + 1) * size]
		term := sandwich(rho, k, dim)
		for idx in 0 .. size {
			out[idx] += term[idx]
		}
	}
	return out
}

// embed_operator places a `op_dim`-by-`op_dim` operator on `targets` of a
// `num_qubits`-qubit space, returning the full `2**num_qubits` square matrix.
@[vc_fn]
pub fn embed_operator(operator []f64, op_dim int, targets []int, num_qubits int) []f64 {
	dim := 1 << num_qubits
	mut out := []f64{len: 2 * dim * dim, init: 0.0}
	mut tmask := 0
	for t in targets {
		tmask |= 1 << (num_qubits - 1 - t)
	}
	for i in 0 .. dim {
		mut sub := 0
		for t in targets {
			sub = (sub << 1) | ((i >> (num_qubits - 1 - t)) & 1)
		}
		base := i & ~tmask
		for row in 0 .. op_dim {
			mut new_i := base
			for j, t in targets {
				new_i |= ((row >> (targets.len - 1 - j)) & 1) << (num_qubits - 1 - t)
			}
			out[2 * (new_i * dim + i)] = operator[2 * (row * op_dim + sub)]
			out[2 * (new_i * dim + i) + 1] = operator[2 * (row * op_dim + sub) + 1]
		}
	}
	return out
}

fn outcome_index(index int, targets []int, num_qubits int) int {
	mut sub := 0
	for t in targets {
		sub = (sub << 1) | ((index >> (num_qubits - 1 - t)) & 1)
	}
	return sub
}

// state_probabilities returns the probability of each outcome of measuring
// `targets` on a state vector, in big-endian outcome order.
@[vc_fn]
pub fn state_probabilities(state []f64, targets []int, num_qubits int) []f64 {
	mut probs := []f64{len: 1 << targets.len, init: 0.0}
	n := 1 << num_qubits
	for i in 0 .. n {
		sub := outcome_index(i, targets, num_qubits)
		re := state[2 * i]
		im := state[2 * i + 1]
		probs[sub] += re * re + im * im
	}
	return probs
}

// collapse_state projects a state vector onto ``outcome`` of ``targets`` and
// renormalizes it.
@[vc_fn]
pub fn collapse_state(state []f64, targets []int, outcome int, num_qubits int) []f64 {
	n := 1 << num_qubits
	mut norm := 0.0
	for i in 0 .. n {
		if outcome_index(i, targets, num_qubits) == outcome {
			re := state[2 * i]
			im := state[2 * i + 1]
			norm += re * re + im * im
		}
	}
	scale := if norm > 0.0 { 1.0 / math.sqrt(norm) } else { 0.0 }
	mut out := []f64{len: 2 * n, init: 0.0}
	for i in 0 .. n {
		if outcome_index(i, targets, num_qubits) == outcome {
			out[2 * i] = state[2 * i] * scale
			out[2 * i + 1] = state[2 * i + 1] * scale
		}
	}
	return out
}

// density_probabilities returns the probability of each outcome of measuring
// `targets` on a density matrix.
@[vc_fn]
pub fn density_probabilities(rho []f64, targets []int, num_qubits int) []f64 {
	mut probs := []f64{len: 1 << targets.len, init: 0.0}
	n := 1 << num_qubits
	for i in 0 .. n {
		probs[outcome_index(i, targets, num_qubits)] += rho[2 * (i * n + i)]
	}
	return probs
}

// density_collapse projects a density matrix onto ``outcome`` of ``targets``
// and renormalizes it.
@[vc_fn]
pub fn density_collapse(rho []f64, targets []int, outcome int, num_qubits int) []f64 {
	n := 1 << num_qubits
	mut norm := 0.0
	for i in 0 .. n {
		if outcome_index(i, targets, num_qubits) == outcome {
			norm += rho[2 * (i * n + i)]
		}
	}
	scale := if norm > 0.0 { 1.0 / norm } else { 0.0 }
	mut out := []f64{len: 2 * n * n, init: 0.0}
	for i in 0 .. n {
		if outcome_index(i, targets, num_qubits) != outcome {
			continue
		}
		for j in 0 .. n {
			if outcome_index(j, targets, num_qubits) == outcome {
				out[2 * (i * n + j)] = rho[2 * (i * n + j)] * scale
				out[2 * (i * n + j) + 1] = rho[2 * (i * n + j) + 1] * scale
			}
		}
	}
	return out
}

// pauli_entry returns the complex entry of a single-qubit Pauli matrix, with
// 0=I, 1=X, 2=Y, 3=Z.
@[inline]
fn pauli_entry(digit int, row int, col int) (f64, f64) {
	if digit == 0 {
		return if row == col { 1.0, 0.0 } else { 0.0, 0.0 }
	}
	if digit == 1 {
		return if row != col { 1.0, 0.0 } else { 0.0, 0.0 }
	}
	if digit == 2 {
		if row == 0 && col == 1 {
			return 0.0, 1.0
		}
		if row == 1 && col == 0 {
			return 0.0, -1.0
		}
		return 0.0, 0.0
	}
	if row == 0 && col == 0 {
		return 1.0, 0.0
	}
	if row == 1 && col == 1 {
		return -1.0, 0.0
	}
	return 0.0, 0.0
}

// pauli_expansion returns the coefficient (real, imaginary) of every Pauli
// string in the expansion of `matrix`, in the order of products of (I, X, Y, Z)
// with the first qubit most significant.
@[vc_fn]
pub fn pauli_expansion(matrix []f64, num_qubits int) []f64 {
	dim := 1 << num_qubits
	combos := 1 << (2 * num_qubits)
	mut out := []f64{len: 2 * combos, init: 0.0}
	for combo in 0 .. combos {
		mut tr_re := 0.0
		mut tr_im := 0.0
		for i in 0 .. dim {
			for j in 0 .. dim {
				mut pr := 1.0
				mut pi := 0.0
				for q in 0 .. num_qubits {
					digit := (combo >> (2 * (num_qubits - 1 - q))) & 3
					bi := (i >> (num_qubits - 1 - q)) & 1
					bj := (j >> (num_qubits - 1 - q)) & 1
					er, ei := pauli_entry(digit, bi, bj)
					nr := pr * er - pi * ei
					ni := pr * ei + pi * er
					pr = nr
					pi = ni
				}
				mr := matrix[2 * (j * dim + i)]
				mi := matrix[2 * (j * dim + i) + 1]
				tr_re += pr * mr - pi * mi
				tr_im += pr * mi + pi * mr
			}
		}
		out[2 * combo] = tr_re / f64(dim)
		out[2 * combo + 1] = tr_im / f64(dim)
	}
	return out
}
