module vqubit_native

// The Clifford (stabilizer) simulator of Aaronson and Gottesman.
//
// A stabilizer state is a binary tableau: `2 * n + 1` rows (the last one is
// scratch space) of `2 * n + 1` columns -- `n` x-bits, `n` z-bits and a phase
// bit. Rows `0..n-1` are the destabilizers, rows `n..2n-1` the stabilizers.
//
// The tableau is stored flat (row-major, `(2n+1)` by `n`) so the same routines
// serve both the batched `clifford_run` and a `NativeTableau` whose `@[vc_buffer]`
// fields keep the state for the object's lifetime.

import rand
import math
import vcraft

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

// CliffordOp names the instructions the stabilizer simulator applies. The
// members are declared in the order of the `instr_*` constants above, so the
// integer wire format used by `clifford_run` and the enum stay in step.
pub enum CliffordOp {
	h
	s
	sdg
	x
	y
	z
	cnot
	cz
	swap
	measure
}

struct Tableau {
mut:
	n int
	x []int
	z []int
	r []int
}

fn new_tableau(n int) Tableau {
	rows := 2*n + 1
	mut x := []int{len: rows * n, init: 0}
	mut z := []int{len: rows * n, init: 0}
	mut r := []int{len: rows, init: 0}
	for i in 0 .. n {
		x[i*n + i] = 1
		z[(n+i)*n + i] = 1
	}
	return Tableau{
		n: n
		x: x
		z: z
		r: r
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
		total += g(t.x[i*t.n + j], t.z[i*t.n + j], t.x[h*t.n + j], t.z[h*t.n + j])
	}
	t.r[h] = if total % 4 == 0 { 0 } else { 1 }
	for j in 0 .. t.n {
		t.x[h*t.n + j] ^= t.x[i*t.n + j]
		t.z[h*t.n + j] ^= t.z[i*t.n + j]
	}
}

fn (mut t Tableau) hadamard(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i*t.n + a] * t.z[i*t.n + a]
		value := t.x[i*t.n + a]
		t.x[i*t.n + a] = t.z[i*t.n + a]
		t.z[i*t.n + a] = value
	}
}

fn (mut t Tableau) phase(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i*t.n + a] * t.z[i*t.n + a]
		t.z[i*t.n + a] ^= t.x[i*t.n + a]
	}
}

fn (mut t Tableau) pauli_x(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.z[i*t.n + a]
	}
}

fn (mut t Tableau) pauli_z(a int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i*t.n + a]
	}
}

fn (mut t Tableau) pauli_y(a int) {
	t.pauli_x(a)
	t.pauli_z(a)
}

fn (mut t Tableau) cnot(a int, b int) {
	for i in 0 .. 2*t.n {
		t.r[i] ^= t.x[i*t.n + a] * t.z[i*t.n + b] * (t.x[i*t.n + b] ^ t.z[i*t.n + a] ^ 1)
		t.x[i*t.n + b] ^= t.x[i*t.n + a]
		t.z[i*t.n + a] ^= t.z[i*t.n + b]
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

// measure measures qubit `a` in the Z basis, using `u` (a uniform draw in
// [0, 1)) for a random outcome, and returns the bit.
fn (mut t Tableau) measure(a int, u f64) int {
	n := t.n
	mut p := -1
	for i in n .. 2*n {
		if t.x[i*n + a] == 1 {
			p = i
			break
		}
	}
	if p != -1 {
		for i in 0 .. 2*n {
			if i != p && t.x[i*n + a] == 1 {
				t.rowsum(i, p)
			}
		}
		for j in 0 .. n {
			t.x[(p-n)*n + j] = t.x[p*n + j]
			t.z[(p-n)*n + j] = t.z[p*n + j]
		}
		t.r[p-n] = t.r[p]
		for j in 0 .. n {
			t.x[p*n + j] = 0
			t.z[p*n + j] = 0
		}
		t.z[p*n + a] = 1
		t.r[p] = if u < 0.5 { 0 } else { 1 }
		return t.r[p]
	}
	scratch := 2 * n
	for j in 0 .. n {
		t.x[scratch*n + j] = 0
		t.z[scratch*n + j] = 0
	}
	t.r[scratch] = 0
	for i in 0 .. n {
		if t.x[i*n + a] == 1 {
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
			out[i*side + j] = t.x[i*t.n + j]
			out[i*side + t.n + j] = t.z[i*t.n + j]
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
				bit := t.measure(a, rand.f64())
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

// NativeTableau is a stabilizer state whose tableau lives in `@[vc_buffer]`
// fields for the object's lifetime, so a circuit can be applied one gate at a
// time without rebuilding the state.
@[vc_class]
pub struct NativeTableau {
mut:
	@[vc_buffer] x []int
	@[vc_buffer] z []int
	@[vc_buffer] r []int
	@[vc_field] num_qubits int
}

// new_native_tableau creates the tableau of `|0...0>`.
@[vc_fn]
pub fn new_native_tableau(num_qubits int) &NativeTableau {
	n := num_qubits
	rows := 2*n + 1
	mut tableau := &NativeTableau{
		x:          vcraft.buffer_new_int(rows * n)
		z:          vcraft.buffer_new_int(rows * n)
		r:          vcraft.buffer_new_int(rows)
		num_qubits: n
	}
	tableau.reset()
	return tableau
}

fn (t &NativeTableau) view() Tableau {
	return Tableau{
		n: t.num_qubits
		x: t.x
		z: t.z
		r: t.r
	}
}

// reset returns the tableau to `|0...0>`.
@[vc_method]
pub fn (mut t NativeTableau) reset() {
	n := t.num_qubits
	for k in 0 .. t.x.len {
		t.x[k] = 0
	}
	for k in 0 .. t.z.len {
		t.z[k] = 0
	}
	for k in 0 .. t.r.len {
		t.r[k] = 0
	}
	for i in 0 .. n {
		t.x[i*n + i] = 1
		t.z[(n+i)*n + i] = 1
	}
}

// apply runs one instruction (a `CliffordOp`) in place.
@[vc_method]
pub fn (mut t NativeTableau) apply(op CliffordOp, a int, b int) {
	mut tableau := t.view()
	tableau.apply(int(op), a, b)
}

// measure measures qubit `a` with a random outcome drawn from `u` in [0, 1).
@[vc_method]
pub fn (mut t NativeTableau) measure(a int, u f64) int {
	mut tableau := t.view()
	return tableau.measure(a, u)
}

// state_vector reconstructs the state, up to a global phase.
@[vc_method]
pub fn (t &NativeTableau) state_vector() []f64 {
	tableau := t.view()
	return clifford_state_vector(tableau.to_flat(), t.num_qubits)
}

// stabilizer_tableau returns the flattened tableau, for inspection.
@[vc_method]
pub fn (t &NativeTableau) stabilizer_tableau() []int {
	tableau := t.view()
	return tableau.to_flat()
}
