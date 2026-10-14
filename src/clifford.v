module vqubit_native

// The Clifford (stabilizer) simulator of Aaronson and Gottesman.
//
// A stabilizer state is a binary tableau: `2 * n + 1` rows (the last one is
// scratch space) of `2 * n + 1` columns -- `n` x-bits, `n` z-bits and a phase
// bit. Rows `0..n-1` are the destabilizers, rows `n..2n-1` the stabilizers.
// Circuits arrive as a flat stream of `[code, a, b]` triples, so the whole
// simulation runs in one call with no per-gate crossings.

import rand
import math

// Instruction codes shared with the Python side.
const instr_h = 0
const instr_s = 1
const instr_sdg = 2
const instr_x = 3
const instr_y = 4
const instr_z = 5
const instr_cnot = 6
const instr_cz = 7
const instr_swap = 8
const instr_measure = 9

struct Tableau {
mut:
	n    int
	rows int
	x    [][]int
	z    [][]int
	r    []int
}

fn new_tableau(n int) Tableau {
	rows := 2*n + 1
	mut x := [][]int{}
	mut z := [][]int{}
	for _ in 0 .. rows {
		x << []int{len: n, init: 0}
		z << []int{len: n, init: 0}
	}
	mut r := []int{len: rows, init: 0}
	for i in 0 .. n {
		x[i][i] = 1
		z[n+i][i] = 1
	}
	return Tableau{
		n:    n
		rows: rows
		x:    x
		z:    z
		r:    r
	}
}

fn g(x1 int, z1 int, x2 int, z2 int) int {
	if x1 == 0 && z1 == 0 {
		return 0
	}
	if x1 == 1 && z1 == 1 {
		return z2 - x2
	}
	if x1 == 1 && z1 == 0 {
		return z2 * (2*x2 - 1)
	}
	return x2 * (1 - 2*z2)
}

fn (mut t Tableau) rowsum(h int, i int) {
	mut total := 2*t.r[h] + 2*t.r[i]
	for j in 0 .. t.n {
		total += g(t.x[i][j], t.z[i][j], t.x[h][j], t.z[h][j])
	}
	t.r[h] = if total % 4 == 0 { 0 } else { 1 }
	for j in 0 .. t.n {
		t.x[h][j] ^= t.x[i][j]
		t.z[h][j] ^= t.z[i][j]
	}
}

fn (mut t Tableau) hadamard(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i][a] * t.z[i][a]
		x := t.x[i][a]
		t.x[i][a] = t.z[i][a]
		t.z[i][a] = x
	}
}

fn (mut t Tableau) phase(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i][a] * t.z[i][a]
		t.z[i][a] ^= t.x[i][a]
	}
}

fn (mut t Tableau) pauli_x(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.z[i][a]
	}
}

fn (mut t Tableau) pauli_z(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i][a]
	}
}

fn (mut t Tableau) pauli_y(a int) {
	t.pauli_x(a)
	t.pauli_z(a)
}

fn (mut t Tableau) cnot(a int, b int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i][a] * t.z[i][b] * (t.x[i][b] ^ t.z[i][a] ^ 1)
		t.x[i][b] ^= t.x[i][a]
		t.z[i][a] ^= t.z[i][b]
	}
}

fn (mut t Tableau) cz(a int, b int) {
	t.hadamard(b)
	t.cnot(a, b)
	t.hadamard(b)
}

fn (mut t Tableau) swap(a int, b int) {
	t.cnot(a, b)
	t.cnot(b, a)
	t.cnot(a, b)
}

fn (mut t Tableau) measure(a int) int {
	n := t.n
	mut p := -1
	for i in n .. 2*n {
		if t.x[i][a] == 1 {
			p = i
			break
		}
	}
	if p != -1 {
		for i in 0 .. 2*n {
			if i != p && t.x[i][a] == 1 {
				t.rowsum(i, p)
			}
		}
		for j in 0 .. n {
			t.x[p-n][j] = t.x[p][j]
			t.z[p-n][j] = t.z[p][j]
		}
		t.r[p-n] = t.r[p]
		for j in 0 .. n {
			t.x[p][j] = 0
			t.z[p][j] = 0
		}
		t.z[p][a] = 1
		t.r[p] = rand.intn(2) or { 0 }
		return t.r[p]
	}
	scratch := 2 * n
	for j in 0 .. n {
		t.x[scratch][j] = 0
		t.z[scratch][j] = 0
	}
	t.r[scratch] = 0
	for i in 0 .. n {
		if t.x[i][a] == 1 {
			t.rowsum(scratch, i+n)
		}
	}
	return t.r[scratch]
}

fn (mut t Tableau) apply(code int, a int, b int) {
	match code {
		instr_h { t.hadamard(a) }
		instr_s { t.phase(a) }
		instr_sdg {
			t.phase(a)
			t.phase(a)
			t.phase(a)
		}
		instr_x { t.pauli_x(a) }
		instr_y { t.pauli_y(a) }
		instr_z { t.pauli_z(a) }
		instr_cnot { t.cnot(a, b) }
		instr_cz { t.cz(a, b) }
		instr_swap { t.swap(a, b) }
		else {}
	}
}

fn (t &Tableau) to_flat() []int {
	side := 2*t.n + 1
	mut out := []int{len: side * side, init: 0}
	for i in 0 .. side {
		for j in 0 .. t.n {
			out[i*side + j] = t.x[i][j]
			out[i*side + t.n + j] = t.z[i][j]
		}
		out[i*side + 2*t.n] = t.r[i]
	}
	return out
}

// clifford_run runs a Clifford circuit `repetitions` times and returns the
// measured bits, row-major (`repetitions` rows of one bit per measurement
// instruction, in order), with each measurement's invert applied.
@[vc_fn]
pub fn clifford_run(instructions []int, num_qubits int, invert_mask []int, repetitions int, seed int) []int {
	n := num_qubits
	mut num_meas := 0
	for i := 0; i * 3 + 2 < instructions.len; i++ {
		if instructions[i*3] == instr_measure {
			num_meas++
		}
	}
	mut out := []int{len: repetitions * num_meas, init: 0}
	if repetitions <= 0 || num_meas == 0 {
		return out
	}
	rand.seed([u32(seed), u32(seed >> 16)])
	for rep in 0 .. repetitions {
		mut t := new_tableau(n)
		mut m := 0
		for i := 0; i * 3 + 2 < instructions.len; i++ {
			code := instructions[i*3]
			a := instructions[i*3 + 1]
			b := instructions[i*3 + 2]
			if code == instr_measure {
				bit := t.measure(a)
				inv := if m < invert_mask.len { invert_mask[m] } else { 0 }
				out[rep*num_meas + m] = bit ^ inv
				m++
			} else {
				t.apply(code, a, b)
			}
		}
	}
	return out
}

// clifford_tableau returns the final tableau of a Clifford circuit, flattened
// into a `(2n+1)`-by-`(2n+1)` row-major matrix.
@[vc_fn]
pub fn clifford_tableau(instructions []int, num_qubits int) []int {
	mut t := new_tableau(num_qubits)
	for i := 0; i * 3 + 2 < instructions.len; i++ {
		t.apply(instructions[i*3], instructions[i*3 + 1], instructions[i*3 + 2])
	}
	return t.to_flat()
}

fn pauli_matrix(xbit int, zbit int) []f64 {
	if xbit == 1 && zbit == 0 {
		return [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
	}
	if xbit == 0 && zbit == 1 {
		return [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0]
	}
	if xbit == 1 && zbit == 1 {
		return [0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0]
	}
	return [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
}

fn complex_identity(dim int) []f64 {
	mut out := []f64{len: 2 * dim * dim, init: 0.0}
	for i in 0 .. dim {
		out[2*(i*dim + i)] = 1.0
	}
	return out
}

// clifford_state_vector reconstructs a state vector from a flattened tableau,
// up to a global phase.
@[vc_fn]
pub fn clifford_state_vector(tableau []int, num_qubits int) []f64 {
	n := num_qubits
	side := 2*n + 1
	dim := 1 << n
	mut acc := complex_identity(dim)
	for k in 0 .. n {
		row := n + k
		mut mat := []f64{len: 2, init: 0.0}
		mat[0] = 1.0
		mut mdim := 1
		for j in 0 .. n {
			xb := tableau[row*side + j]
			zb := tableau[row*side + n + j]
			mat = kron(mat, mdim, pauli_matrix(xb, zb), 2)
			mdim *= 2
		}
		sign := if tableau[row*side + 2*n] == 1 { -1.0 } else { 1.0 }
		// projector = I + sign * mat
		mut proj := []f64{len: 2 * dim * dim, init: 0.0}
		for i in 0 .. dim {
			for j in 0 .. dim {
				idx := 2 * (i * dim + j)
				proj[idx] = sign * mat[idx] + if i == j { 1.0 } else { 0.0 }
				proj[idx + 1] = sign * mat[idx + 1]
			}
		}
		acc = matmul(acc, proj, dim)
	}
	for col in 0 .. dim {
		mut norm := 0.0
		for i in 0 .. dim {
			re := acc[2*(i*dim + col)]
			im := acc[2*(i*dim + col) + 1]
			norm += re*re + im*im
		}
		if norm > 1e-18 {
			scale := 1.0 / math.sqrt(norm)
			mut out := []f64{len: 2 * dim, init: 0.0}
			for i in 0 .. dim {
				out[2*i] = acc[2*(i*dim + col)] * scale
				out[2*i + 1] = acc[2*(i*dim + col) + 1] * scale
			}
			return out
		}
	}
	return []f64{len: 2 * dim, init: 0.0}
}
