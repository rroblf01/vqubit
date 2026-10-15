module vqubit_native

// State-vector and density-matrix state that lives in `@[vc_buffer]` fields.
//
// A simulation keeps one of these for its lifetime and applies gates to it in
// place, so the amplitudes are never converted to Python objects between calls.

import vcraft

// NativeState is the state vector of a circuit, as a flat interleaved complex
// buffer (`[re0, im0, re1, im1, ...]`).
@[vc_class]
pub struct NativeState {
mut:
	@[vc_buffer] state []f64
	@[vc_field] num_qubits int
}

// new_native_state creates the state `|0...0>` of `num_qubits` qubits.
@[vc_fn]
pub fn new_native_state(num_qubits int) &NativeState {
	dim := 1 << num_qubits
	return &NativeState{
		state:      vcraft.buffer_new_f64(2 * dim)
		num_qubits: num_qubits
	}
}

// set_amplitudes overwrites the state with `values` (interleaved complex).
@[vc_method]
pub fn (mut s NativeState) set_amplitudes(values []f64) {
	for i in 0 .. s.state.len {
		s.state[i] = if i < values.len { values[i] } else { 0.0 }
	}
}

// apply_gate applies a gate matrix to `targets`, in place.
@[vc_method]
pub fn (mut s NativeState) apply_gate(gate []f64, targets []int) {
	apply_gate_into(mut s.state, gate, s.num_qubits, targets)
}

// reset sets the state to `|0...0>` without building a Python list.
@[vc_method]
pub fn (mut s NativeState) reset() {
	for i in 0 .. s.state.len {
		s.state[i] = 0.0
	}
	if s.state.len > 0 {
		s.state[0] = 1.0
	}
}

// amplitudes returns the state as an interleaved complex buffer.
@[vc_method]
pub fn (s &NativeState) amplitudes() []f64 {
	return s.state.clone()
}

// probabilities returns the probability of each outcome of measuring `targets`.
@[vc_method]
pub fn (s &NativeState) probabilities(targets []int) []f64 {
	return state_probabilities(s.state, targets, s.num_qubits)
}

// measure samples `targets` with a uniform draw `u` in [0, 1), collapses the
// state and returns the outcome.
@[vc_method]
pub fn (mut s NativeState) measure(targets []int, u f64) int {
	probabilities := state_probabilities(s.state, targets, s.num_qubits)
	outcome := sample_outcome(probabilities, u)
	collapsed := collapse_state(s.state, targets, outcome, s.num_qubits)
	for i in 0 .. s.state.len {
		s.state[i] = collapsed[i]
	}
	return outcome
}

// NativeDensity is the density matrix of a circuit, as a flat interleaved
// complex buffer of `dim` by `dim` entries.
@[vc_class]
pub struct NativeDensity {
mut:
	@[vc_buffer] rho []f64
	@[vc_field] num_qubits int
}

// new_native_density creates the density matrix `|0...0><0...0|`.
@[vc_fn]
pub fn new_native_density(num_qubits int) &NativeDensity {
	dim := 1 << num_qubits
	return &NativeDensity{
		rho:        vcraft.buffer_new_f64(2 * dim * dim)
		num_qubits: num_qubits
	}
}

// set_from_amplitudes sets the density matrix to `|psi><psi|` of `values`.
@[vc_method]
pub fn (mut d NativeDensity) set_from_amplitudes(values []f64) {
	rho := state_to_density(values)
	for i in 0 .. d.rho.len {
		d.rho[i] = rho[i]
	}
}

// apply_channel applies the Kraus operators in `operators` (concatenated) to
// `targets`: `rho -> sum_i K_i rho K_i^dagger`.
@[vc_method]
pub fn (mut d NativeDensity) apply_channel(operators []f64, num_operators int, targets []int) {
	op_dim := op_dimension(operators.len / num_operators)
	channel_subspace_into(mut d.rho, operators, num_operators, op_dim, targets, d.num_qubits)
}

// reset sets the density matrix to `|0...0><0...0|` without a Python list.
@[vc_method]
pub fn (mut d NativeDensity) reset() {
	for i in 0 .. d.rho.len {
		d.rho[i] = 0.0
	}
	if d.rho.len > 0 {
		d.rho[0] = 1.0
	}
}

// density_matrix returns the density matrix as an interleaved complex buffer.
@[vc_method]
pub fn (d &NativeDensity) density_matrix() []f64 {
	return d.rho.clone()
}

// probabilities returns the probability of each outcome of measuring `targets`.
@[vc_method]
pub fn (d &NativeDensity) probabilities(targets []int) []f64 {
	return density_probabilities(d.rho, targets, d.num_qubits)
}

// measure samples `targets` with a uniform draw `u` in [0, 1), collapses the
// state and returns the outcome.
@[vc_method]
pub fn (mut d NativeDensity) measure(targets []int, u f64) int {
	probabilities := density_probabilities(d.rho, targets, d.num_qubits)
	outcome := sample_outcome(probabilities, u)
	collapsed := density_collapse(d.rho, targets, outcome, d.num_qubits)
	for i in 0 .. d.rho.len {
		d.rho[i] = collapsed[i]
	}
	return outcome
}

// op_dimension recovers the side of a square complex matrix from the length of
// its interleaved buffer.
fn op_dimension(complex_length int) int {
	mut d := 1
	for 2*d*d < complex_length {
		d++
	}
	return d
}
