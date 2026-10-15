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
		out_base := 2 * i * dim
		a_base := 2 * i * dim
		for k in 0 .. dim {
			ar := a[a_base + 2 * k]
			ai := a[a_base + 2 * k + 1]
			if ar == 0.0 && ai == 0.0 {
				continue
			}
			b_base := 2 * k * dim
			for j in 0 .. dim {
				br := b[b_base + 2 * j]
				bi := b[b_base + 2 * j + 1]
				out[out_base + 2 * j] += ar * br - ai * bi
				out[out_base + 2 * j + 1] += ar * bi + ai * br
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
// apply_gate_1q applies a 2x2 gate on one qubit with contiguous inner runs.
fn apply_gate_1q(mut state []f64, m []f64, num_qubits int, t int) {
	q := num_qubits - 1 - t
	stride := 1 << q
	block := stride << 1
	n := 1 << num_qubits
	m00r := m[0]
	m00i := m[1]
	m01r := m[2]
	m01i := m[3]
	m10r := m[4]
	m10i := m[5]
	m11r := m[6]
	m11i := m[7]
	unsafe {
		mut base := 0
		for base < n {
			for low in 0 .. stride {
				i0 := 2 * (base + low)
				i1 := i0 + 2 * stride
				a0r := state[i0]
				a0i := state[i0 + 1]
				a1r := state[i1]
				a1i := state[i1 + 1]
				state[i0] = m00r * a0r - m00i * a0i + m01r * a1r - m01i * a1i
				state[i0 + 1] = m00r * a0i + m00i * a0r + m01r * a1i + m01i * a1r
				state[i1] = m10r * a0r - m10i * a0i + m11r * a1r - m11i * a1i
				state[i1 + 1] = m10r * a0i + m10i * a0r + m11r * a1i + m11i * a1r
			}
			base += block
		}
	}
}

// apply_gate_2q_adjacent applies a 4x4 gate on two adjacent qubits with
// contiguous inner runs; `low` is the less significant of the two.
fn apply_gate_2q_adjacent(mut state []f64, m []f64, num_qubits int, low int) {
	q := num_qubits - 1 - low
	stride := 1 << q
	block := stride << 2
	n := 1 << num_qubits
	mut gr := []f64{len: 16, init: 0.0}
	mut gi := []f64{len: 16, init: 0.0}
	for idx in 0 .. 16 {
		gr[idx] = m[2 * idx]
		gi[idx] = m[2 * idx + 1]
	}
	// Permutation fast path (CNOT, CZ, SWAP, ...): gather then scatter.
	mut perm := []int{len: 4, init: -1}
	mut sign := []f64{len: 4, init: 1.0}
	mut is_perm := true
	for row in 0 .. 4 {
		mut count := 0
		mut col := 0
		mut sgn := 1.0
		for c in 0 .. 4 {
			if gr[row * 4 + c] * gr[row * 4 + c] + gi[row * 4 + c] * gi[row * 4 + c] > 1e-18 {
				count++
				col = c
				if math.abs(gi[row * 4 + c]) > 1e-15 {
					is_perm = false
				}
				sgn = if gr[row * 4 + c] < 0 { -1.0 } else { 1.0 }
			}
		}
		if count != 1 || !is_perm {
			is_perm = false
			break
		}
		perm[row] = col
		sign[row] = sgn
	}
	if is_perm {
		mut vals_r := []f64{len: 4, init: 0.0}
		mut vals_i := []f64{len: 4, init: 0.0}
		unsafe {
			mut base := 0
			for base < n {
				for l in 0 .. stride {
					i0 := 2 * (base + l)
					vals_r[0] = state[i0]
					vals_i[0] = state[i0 + 1]
					vals_r[1] = state[i0 + 2 * stride]
					vals_i[1] = state[i0 + 2 * stride + 1]
					vals_r[2] = state[i0 + 4 * stride]
					vals_i[2] = state[i0 + 4 * stride + 1]
					vals_r[3] = state[i0 + 6 * stride]
					vals_i[3] = state[i0 + 6 * stride + 1]
					state[i0] = sign[0] * vals_r[perm[0]]
					state[i0 + 1] = sign[0] * vals_i[perm[0]]
					state[i0 + 2 * stride] = sign[1] * vals_r[perm[1]]
					state[i0 + 2 * stride + 1] = sign[1] * vals_i[perm[1]]
					state[i0 + 4 * stride] = sign[2] * vals_r[perm[2]]
					state[i0 + 4 * stride + 1] = sign[2] * vals_i[perm[2]]
					state[i0 + 6 * stride] = sign[3] * vals_r[perm[3]]
					state[i0 + 6 * stride + 1] = sign[3] * vals_i[perm[3]]
				}
				base += block
			}
		}
		return
	}
	mut i := []int{len: 4, init: 0}
	mut xr := []f64{len: 4, init: 0.0}
	mut xi := []f64{len: 4, init: 0.0}
	unsafe {
		mut base := 0
		for base < n {
			for l in 0 .. stride {
				i[0] = 2 * (base + l)
				i[1] = i[0] + 2 * stride
				i[2] = i[0] + 4 * stride
				i[3] = i[0] + 6 * stride
				for r in 0 .. 4 {
					xr[r] = state[i[r]]
					xi[r] = state[i[r] + 1]
				}
				for row in 0 .. 4 {
					mut sr := 0.0
					mut si := 0.0
					for r in 0 .. 4 {
						sr += gr[row * 4 + r] * xr[r] - gi[row * 4 + r] * xi[r]
						si += gr[row * 4 + r] * xi[r] + gi[row * 4 + r] * xr[r]
					}
					state[i[row]] = sr
					state[i[row] + 1] = si
				}
			}
			base += block
		}
	}
}

// apply_gate_into applies a gate matrix to `targets` of `state` in place.
//
// The target subspace of every "base" index is disjoint from the others, so the
// gather/apply/scatter can run directly on the buffer without a copy.
fn apply_gate_into(mut state []f64, gate []f64, num_qubits int, targets []int) {
	if num_qubits <= 0 || targets.len == 0 {
		return
	}
	if targets.len == 1 {
		apply_gate_1q(mut state, gate, num_qubits, targets[0])
		return
	}
	if targets.len == 2 && (targets[0] - targets[1] == 1 || targets[1] - targets[0] == 1) {
		low := if targets[0] > targets[1] { targets[0] } else { targets[1] }
		high := if targets[0] > targets[1] { targets[1] } else { targets[0] }
		// The gate matrix rows are ordered with the first target most significant;
		// the adjacent helper expects the less significant target first, so swap
		// the basis labels when needed.
		if targets[0] > targets[1] {
			swapped := swap_2q_basis(gate)
			apply_gate_2q_adjacent(mut state, swapped, num_qubits, low)
		} else {
			apply_gate_2q_adjacent(mut state, gate, num_qubits, low)
		}
		_ = high
		return
	}
	n := 1 << num_qubits
	k := targets.len
	g := 1 << k
	mut off := []int{len: g, init: 0}
	mut tmask := 0
	for r in 0 .. g {
		mut o := 0
		for j in 0 .. k {
			if ((r >> (k - 1 - j)) & 1) == 1 {
				o |= 1 << (num_qubits - 1 - targets[j])
			}
		}
		off[r] = o
		tmask |= o
	}
	mut xr := []f64{len: g, init: 0.0}
	mut xi := []f64{len: g, init: 0.0}
	// Fast path: a permutation matrix (CNOT, CZ, SWAP, X, ...) just moves
	// amplitudes, so gather then scatter without any multiplications.
	mut perm := []int{len: g, init: -1}
	mut sign := []f64{len: g, init: 1.0}
	mut is_perm := true
	for row in 0 .. g {
		mut count := 0
		mut col := 0
		mut sgn := 1.0
		for c in 0 .. g {
			re := gate[2 * (row * g + c)]
			im := gate[2 * (row * g + c) + 1]
			if re * re + im * im > 1e-18 {
				count++
				col = c
				if math.abs(im) > 1e-15 {
					is_perm = false
				}
				sgn = if re < 0 { -1.0 } else { 1.0 }
			}
		}
		if count != 1 || !is_perm {
			is_perm = false
			break
		}
		perm[row] = col
		sign[row] = sgn
	}
	if is_perm {
		unsafe {
			for base in 0 .. n {
				if (base & tmask) != 0 {
					continue
				}
				for row in 0 .. g {
					idx := base | off[row]
					xr[row] = state[2 * idx]
					xi[row] = state[2 * idx + 1]
				}
				for row in 0 .. g {
					idx := base | off[row]
					state[2 * idx] = sign[row] * xr[perm[row]]
					state[2 * idx + 1] = sign[row] * xi[perm[row]]
				}
			}
		}
		return
	}
	unsafe {
		for base in 0 .. n {
			if (base & tmask) != 0 {
				continue
			}
			for r in 0 .. g {
				idx := base | off[r]
				xr[r] = state[2 * idx]
				xi[r] = state[2 * idx + 1]
			}
			for row in 0 .. g {
				mut sr := 0.0
				mut si := 0.0
				gate_row := 2 * row * g
				for r in 0 .. g {
					gr := gate[gate_row + 2 * r]
					gi := gate[gate_row + 2 * r + 1]
					sr += gr * xr[r] - gi * xi[r]
					si += gr * xi[r] + gi * xr[r]
				}
				idx := base | off[row]
				state[2 * idx] = sr
				state[2 * idx + 1] = si
			}
		}
	}
}

// swap_2q_basis returns the 4x4 matrix with the two qubit axes exchanged.
fn swap_2q_basis(m []f64) []f64 {
	mut out := []f64{len: 32, init: 0.0}
	for row in 0 .. 4 {
		r0 := (row >> 1) & 1
		r1 := row & 1
		new_row := (r1 << 1) | r0
		for col in 0 .. 4 {
			c0 := (col >> 1) & 1
			c1 := col & 1
			new_col := (c1 << 1) | c0
			out[2 * (new_row * 4 + new_col)] = m[2 * (row * 4 + col)]
			out[2 * (new_row * 4 + new_col) + 1] = m[2 * (row * 4 + col) + 1]
		}
	}
	return out
}

@[vc_fn]
pub fn apply_gate(state []f64, gate []f64, num_qubits int, targets []int) []f64 {
	mut out := state.clone()
	apply_gate_into(mut out, gate, num_qubits, targets)
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

// state_to_density returns the density matrix `|psi> <psi|` of a state vector.
@[vc_fn]
pub fn state_to_density(state []f64) []f64 {
	dim := state.len / 2
	mut out := []f64{len: 2 * dim * dim, init: 0.0}
	for i in 0 .. dim {
		re_i := state[2 * i]
		im_i := state[2 * i + 1]
		for j in 0 .. dim {
			re_j := state[2 * j]
			out[2 * (i * dim + j)] = re_i * re_j + im_i * state[2 * j + 1]
			out[2 * (i * dim + j) + 1] = im_i * re_j - re_i * state[2 * j + 1]
		}
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
		out_base := 2 * i * dim
		u_base := 2 * i * dim
		for k in 0 .. dim {
			ur := u[u_base + 2 * k]
			ui := u[u_base + 2 * k + 1]
			if ur == 0.0 && ui == 0.0 {
				continue
			}
			a_base := 2 * k * dim
			for j in 0 .. dim {
				ar := a[a_base + 2 * j]
				ai := a[a_base + 2 * j + 1]
				tmp[out_base + 2 * j] += ur * ar - ui * ai
				tmp[out_base + 2 * j + 1] += ur * ai + ui * ar
			}
		}
	}
	mut out := []f64{len: size, init: 0.0}
	for i in 0 .. dim {
		out_base := 2 * i * dim
		t_base := 2 * i * dim
		for k in 0 .. dim {
			tr := tmp[t_base + 2 * k]
			ti := tmp[t_base + 2 * k + 1]
			if tr == 0.0 && ti == 0.0 {
				continue
			}
			for j in 0 .. dim {
				// u^dagger[k][j] is the conjugate of u[j][k].
				ur := u[2 * j * dim + 2 * k]
				ui := -u[2 * j * dim + 2 * k + 1]
				out[out_base + 2 * j] += tr * ur - ti * ui
				out[out_base + 2 * j + 1] += tr * ui + ti * ur
			}
		}
	}
	return out
}

// channel_subspace_2x2_into is a specialised in-place 1-qubit channel.
fn channel_subspace_2x2_into(mut rho []f64, ops []f64, num_ops int, targets []int, num_qubits int) {
	dim := 1 << num_qubits
	bit := 1 << (num_qubits - 1 - targets[0])
	bitd := bit * dim
	unsafe {
		for rb in 0 .. dim {
			if (rb & bit) != 0 {
				continue
			}
			for cb in 0 .. dim {
				if (cb & bit) != 0 {
					continue
				}
				base := rb * dim + cb
				// Gather A = [[a00, a01], [a10, a11]].
				i00 := 2 * base
				i01 := i00 + 2 * bit
				i10 := i00 + 2 * bitd
				i11 := i10 + 2 * bit
				a00r := rho[i00]
				a00i := rho[i00 + 1]
				a01r := rho[i01]
				a01i := rho[i01 + 1]
				a10r := rho[i10]
				a10i := rho[i10 + 1]
				a11r := rho[i11]
				a11i := rho[i11 + 1]
				mut b00r := 0.0
				mut b00i := 0.0
				mut b01r := 0.0
				mut b01i := 0.0
				mut b10r := 0.0
				mut b10i := 0.0
				mut b11r := 0.0
				mut b11i := 0.0
				for o in 0 .. num_ops {
					k := ops[o * 8..(o + 1) * 8]
					k00r := k[0]
					k00i := k[1]
					k01r := k[2]
					k01i := k[3]
					k10r := k[4]
					k10i := k[5]
					k11r := k[6]
					k11i := k[7]
					// C = K A.
					c00r := k00r * a00r - k00i * a00i + k01r * a10r - k01i * a10i
					c00i := k00r * a00i + k00i * a00r + k01r * a10i + k01i * a10r
					c01r := k00r * a01r - k00i * a01i + k01r * a11r - k01i * a11i
					c01i := k00r * a01i + k00i * a01r + k01r * a11i + k01i * a11r
					c10r := k10r * a00r - k10i * a00i + k11r * a10r - k11i * a10i
					c10i := k10r * a00i + k10i * a00r + k11r * a10i + k11i * a10r
					c11r := k10r * a01r - k10i * a01i + k11r * a11r - k11i * a11i
					c11i := k10r * a01i + k10i * a01r + k11r * a11i + k11i * a11r
					// B += C K^dagger, with K^dagger[d][j] = conj(K[j][d]).
					b00r += c00r * k00r + c00i * k00i + c01r * k01r + c01i * k01i
					b00i += c00i * k00r - c00r * k00i + c01i * k01r - c01r * k01i
					b01r += c00r * k10r + c00i * k10i + c01r * k11r + c01i * k11i
					b01i += c00i * k10r - c00r * k10i + c01i * k11r - c01r * k11i
					b10r += c10r * k00r + c10i * k00i + c11r * k01r + c11i * k01i
					b10i += c10i * k00r - c10r * k00i + c11i * k01r - c11r * k01i
					b11r += c10r * k10r + c10i * k10i + c11r * k11r + c11i * k11i
					b11i += c10i * k10r - c10r * k10i + c11i * k11r - c11r * k11i
				}
				rho[i00] = b00r
				rho[i00 + 1] = b00i
				rho[i01] = b01r
				rho[i01 + 1] = b01i
				rho[i10] = b10r
				rho[i10 + 1] = b10i
				rho[i11] = b11r
				rho[i11 + 1] = b11i
			}
		}
	}
}

// channel_permutation_into applies a single permutation Kraus operator
// `rho -> P rho P^dagger` by permuting rows and columns.
fn channel_permutation_into(mut rho []f64, perm []int, targets []int, num_qubits int) {
	dim := 1 << num_qubits
	k := targets.len
	op_dim := perm.len
	mut off := []int{len: op_dim, init: 0}
	for r in 0 .. op_dim {
		mut o := 0
		for j in 0 .. k {
			if ((r >> (k - 1 - j)) & 1) == 1 {
				o |= 1 << (num_qubits - 1 - targets[j])
			}
		}
		off[r] = o
	}
	mut index_map := []int{len: dim, init: 0}
	unsafe {
		for i in 0 .. dim {
			mut sub := 0
			for j in 0 .. k {
				sub = (sub << 1) | ((i >> (num_qubits - 1 - targets[j])) & 1)
			}
			index_map[i] = (i & ~off[op_dim - 1]) | off[perm[sub]]
		}
	}
	mut tmp := rho.clone()
	unsafe {
		for i in 0 .. dim {
			mi := index_map[i] * dim
			for j in 0 .. dim {
				mj := index_map[j]
				src := 2 * (mi + mj)
				dst := 2 * (i * dim + j)
				tmp[dst] = rho[src]
				tmp[dst + 1] = rho[src + 1]
			}
		}
	}
	for i in 0 .. rho.len {
		rho[i] = tmp[i]
	}
}

// channel_subspace_into applies `sum_i K_i rho K_i^dagger` on the `targets` of a
// density matrix in place.
//
// It works block by block on the target subspace, so it never builds the full
// `dim`-by-`dim` embedded Kraus operators and runs in `O(dim**2 * op_dim)`.
fn channel_subspace_into(mut rho []f64, ops []f64, num_ops int, op_dim int, targets []int, num_qubits int) {
	if num_qubits <= 0 || targets.len == 0 {
		return
	}
	if num_ops == 1 {
		// A single permutation Kraus (a Clifford permutation) just permutes rows
		// and columns, so it can skip the block arithmetic entirely.
		mut perm := []int{len: op_dim, init: -1}
		mut is_perm := true
		for row in 0 .. op_dim {
			mut count := 0
			mut col := 0
			for c in 0 .. op_dim {
				re := ops[2 * (row * op_dim + c)]
				im := ops[2 * (row * op_dim + c) + 1]
				if re * re + im * im > 1e-18 {
					count++
					col = c
					// The fast path only permutes basis states; signed or complex
					// entries (e.g. the -1 of a CZ) must go through the block
					// arithmetic so the sign is preserved.
					if math.abs(im) > 1e-15 || math.abs(re - 1.0) > 1e-15 {
						is_perm = false
					}
				}
			}
			if count != 1 || !is_perm {
				is_perm = false
				break
			}
			perm[row] = col
		}
		if is_perm {
			channel_permutation_into(mut rho, perm, targets, num_qubits)
			return
		}
	}
	if op_dim == 2 && targets.len == 1 {
		channel_subspace_2x2_into(mut rho, ops, num_ops, targets, num_qubits)
		return
	}
	dim := 1 << num_qubits
	k := targets.len
	mut off := []int{len: op_dim, init: 0}
	mut tmask := 0
	for r in 0 .. op_dim {
		mut o := 0
		for j in 0 .. k {
			if ((r >> (k - 1 - j)) & 1) == 1 {
				o |= 1 << (num_qubits - 1 - targets[j])
			}
		}
		off[r] = o
		tmask |= o
	}
	size := 2 * op_dim * op_dim
	mut ar := []f64{len: op_dim * op_dim, init: 0.0}
	mut ai := []f64{len: op_dim * op_dim, init: 0.0}
	mut br := []f64{len: op_dim * op_dim, init: 0.0}
	mut bi := []f64{len: op_dim * op_dim, init: 0.0}
	unsafe {
	for rb in 0 .. dim {
		if (rb & tmask) != 0 {
			continue
		}
		for cb in 0 .. dim {
			if (cb & tmask) != 0 {
				continue
			}
			base := rb * dim + cb
			for b in 0 .. op_dim {
				row_off := off[b] * dim
				for d in 0 .. op_dim {
					idx := 2 * (base + row_off + off[d])
					ar[b * op_dim + d] = rho[idx]
					ai[b * op_dim + d] = rho[idx + 1]
				}
			}
			for idx in 0 .. op_dim * op_dim {
				br[idx] = 0.0
				bi[idx] = 0.0
			}
			for o in 0 .. num_ops {
				kop := ops[o * size..(o + 1) * size]
				for ri in 0 .. op_dim {
					for rj in 0 .. op_dim {
						mut sr := 0.0
						mut si := 0.0
						for b in 0 .. op_dim {
							kr := kop[2 * (ri * op_dim + b)]
							ki := kop[2 * (ri * op_dim + b) + 1]
							for d in 0 .. op_dim {
								dr := kop[2 * (rj * op_dim + d)]
								di := -kop[2 * (rj * op_dim + d) + 1]
								avr := ar[b * op_dim + d]
								avi := ai[b * op_dim + d]
								mr := avr * dr - avi * di
								mi := avr * di + avi * dr
								sr += kr * mr - ki * mi
								si += kr * mi + ki * mr
							}
						}
						br[ri * op_dim + rj] += sr
						bi[ri * op_dim + rj] += si
					}
				}
			}
			for ri in 0 .. op_dim {
				row_off := off[ri] * dim
				for rj in 0 .. op_dim {
					idx := 2 * (base + row_off + off[rj])
					rho[idx] = br[ri * op_dim + rj]
					rho[idx + 1] = bi[ri * op_dim + rj]
				}
			}
		}
	}
	}
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
	k := targets.len
	mut off := []int{len: op_dim, init: 0}
	mut tmask := 0
	for r in 0 .. op_dim {
		mut o := 0
		for j in 0 .. k {
			if ((r >> (k - 1 - j)) & 1) == 1 {
				o |= 1 << (num_qubits - 1 - targets[j])
			}
		}
		off[r] = o
		tmask |= o
	}
	mut out := []f64{len: 2 * dim * dim, init: 0.0}
	for base in 0 .. dim {
		if (base & tmask) != 0 {
			continue
		}
		for col in 0 .. op_dim {
			icol := base | off[col]
			op_base := 2 * col
			for row in 0 .. op_dim {
				irow := base | off[row]
				slot := 2 * (irow * dim + icol)
				out[slot] = operator[2 * (row * op_dim) + op_base]
				out[slot + 1] = operator[2 * (row * op_dim) + op_base + 1]
			}
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

// sample_outcome returns the index where the cumulative probability first
// exceeds the uniform draw `u`, i.e. one sample from a categorical
// distribution given its probabilities.
@[vc_fn]
pub fn sample_outcome(probabilities []f64, u f64) int {
	mut cumulative := 0.0
	for i, p in probabilities {
		cumulative += p
		if u < cumulative {
			return i
		}
	}
	return probabilities.len - 1
}

// is_identity_matrix reports (as 1 or 0) whether a `dim`-by-`dim` complex
// matrix is the identity within `atol`.
@[vc_fn]
pub fn is_identity_matrix(matrix []f64, dim int, atol f64) int {
	for i in 0 .. dim {
		for j in 0 .. dim {
			re := matrix[2 * (i * dim + j)]
			im := matrix[2 * (i * dim + j) + 1]
			expected := if i == j { 1.0 } else { 0.0 }
			if math.abs(re - expected) > atol || math.abs(im) > atol {
				return 0
			}
		}
	}
	return 1
}

// matrix_allclose reports (as 1 or 0) whether two `dim`-by-`dim` complex
// matrices agree within `atol`.
@[vc_fn]
pub fn matrix_allclose(a []f64, b []f64, dim int, atol f64) int {
	for i in 0 .. dim {
		for j in 0 .. dim {
			idx := 2 * (i * dim + j)
			if math.abs(a[idx] - b[idx]) > atol || math.abs(a[idx + 1] - b[idx + 1]) > atol {
				return 0
			}
		}
	}
	return 1
}
