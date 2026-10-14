module vqubit_native

// Numeric kernels for the KAK (Cartan) decomposition of a two-qubit unitary.
//
// The decomposition works with real symmetric matrices (diagonalized by a
// Jacobi sweep) and small complex matrices (determinants and inverses), and
// splits a local 4x4 gate into a global phase and two 2x2 factors. All of it
// runs here so no linear algebra happens in Python.

import math

// jacobi_eig diagonalizes a real symmetric `n`-by-`n` matrix (row-major) with
// the cyclic Jacobi method. It returns the `n` eigenvalues followed by the
// `n`-by-`n` eigenvector matrix (row-major, columns are eigenvectors).
@[vc_fn]
pub fn jacobi_eig(matrix []f64, n int) []f64 {
	mut a := matrix.clone()
	mut v := []f64{len: n * n, init: 0.0}
	for i in 0 .. n {
		v[i*n + i] = 1.0
	}
	for _ in 0 .. 100 {
		mut off := 0.0
		for i in 0 .. n {
			for j in 0 .. n {
				if i != j {
					off += a[i*n + j] * a[i*n + j]
				}
			}
		}
		if off < 1e-14 {
			break
		}
		for p in 0 .. n {
			for q in p + 1 .. n {
				apq := a[p*n + q]
				if math.abs(apq) < 1e-14 {
					continue
				}
				theta := (a[q*n + q] - a[p*n + p]) / (2 * apq)
				abs_theta := math.abs(theta)
				mut t := 1.0 / (abs_theta + math.sqrt(theta*theta + 1))
				if theta < 0 {
					t = -t
				}
				c := 1.0 / math.sqrt(t*t + 1)
				s := t * c
				for k in 0 .. n {
					akp := a[k*n + p]
					akq := a[k*n + q]
					a[k*n + p] = c*akp - s*akq
					a[k*n + q] = s*akp + c*akq
				}
				for k in 0 .. n {
					apk := a[p*n + k]
					aqk := a[q*n + k]
					a[p*n + k] = c*apk - s*aqk
					a[q*n + k] = s*apk + c*aqk
				}
				for k in 0 .. n {
					vkp := v[k*n + p]
					vkq := v[k*n + q]
					v[k*n + p] = c*vkp - s*vkq
					v[k*n + q] = s*vkp + c*vkq
				}
			}
		}
	}
	mut out := []f64{len: n + n*n, init: 0.0}
	for i in 0 .. n {
		out[i] = a[i*n + i]
	}
	for k in 0 .. n*n {
		out[n + k] = v[k]
	}
	return out
}

// simultaneous_diag jointly diagonalizes two commuting real symmetric
// matrices by a real orthogonal matrix. It returns the `n` complex eigenvalues
// (interleaved) followed by the eigenvector matrix (row-major, columns).
@[vc_fn]
pub fn simultaneous_diag(real []f64, imag []f64, n int) []f64 {
	first := jacobi_eig(real, n)
	mut values_a := []f64{len: n, init: 0.0}
	mut vectors := []f64{len: n*n, init: 0.0}
	for i in 0 .. n {
		values_a[i] = first[i]
	}
	for k in 0 .. n*n {
		vectors[k] = first[n + k]
	}
	mut order := []int{len: n, init: 0}
	for i in 0 .. n {
		order[i] = i
	}
	for i in 1 .. n {
		key := order[i]
		mut j := i - 1
		for j >= 0 && values_a[order[j]] > values_a[key] {
			order[j + 1] = order[j]
			j--
		}
		order[j + 1] = key
	}
	mut groups := [][]int{}
	mut current := [order[0]]
	for gi in 1 .. n {
		idx := order[gi]
		if math.abs(values_a[idx] - values_a[current[0]]) < 1e-9 {
			current << idx
		} else {
			groups << current
			current = [idx]
		}
	}
	groups << current
	for group in groups {
		m := group.len
		if m < 2 {
			continue
		}
		mut sub := []f64{len: m * m, init: 0.0}
		for i in 0 .. m {
			for j in 0 .. m {
				mut acc := 0.0
				for r in 0 .. n {
					for s in 0 .. n {
						acc += vectors[r*n + group[i]] * imag[r*n + s] * vectors[s*n + group[j]]
					}
				}
				sub[i*m + j] = acc
			}
		}
		sub_eig := jacobi_eig(sub, m)
		mut rotated := []f64{len: n * m, init: 0.0}
		for r in 0 .. n {
			for j in 0 .. m {
				mut acc := 0.0
				for k in 0 .. m {
					acc += vectors[r*n + group[k]] * sub_eig[m + k*m + j]
				}
				rotated[r*m + j] = acc
			}
		}
		for j in 0 .. m {
			for r in 0 .. n {
				vectors[r*n + group[j]] = rotated[r*m + j]
			}
		}
	}
	mut out := []f64{len: 2*n + n*n, init: 0.0}
	for j in 0 .. n {
		mut re := 0.0
		mut im := 0.0
		for r in 0 .. n {
			for s in 0 .. n {
				re += vectors[r*n + j] * real[r*n + s] * vectors[s*n + j]
				im += vectors[r*n + j] * imag[r*n + s] * vectors[s*n + j]
			}
		}
		out[2*j] = re
		out[2*j + 1] = im
	}
	for k in 0 .. n*n {
		out[2*n + k] = vectors[k]
	}
	return out
}

// matrix_det returns the determinant of a complex `n`-by-`n` matrix
// (interleaved) as a pair `[re, im]`.
@[vc_fn]
pub fn matrix_det(matrix []f64, n int) []f64 {
	mut a := matrix.clone()
	mut value_re := 1.0
	mut value_im := 0.0
	for col in 0 .. n {
		mut pivot := col
		mut best := math.sqrt(a[2*(col*n + col)] * a[2*(col*n + col)] + a[2*(col*n + col) + 1] * a[2*(col*n + col) + 1])
		for r in col + 1 .. n {
			mag := math.sqrt(a[2*(r*n + col)] * a[2*(r*n + col)] + a[2*(r*n + col) + 1] * a[2*(r*n + col) + 1])
			if mag > best {
				best = mag
				pivot = r
			}
		}
		if best < 1e-14 {
			return [0.0, 0.0]
		}
		if pivot != col {
			for c in 0 .. n {
				i1 := 2 * (col*n + c)
				i2 := 2 * (pivot*n + c)
				tr := a[i1]
				ti := a[i1 + 1]
				a[i1] = a[i2]
				a[i1 + 1] = a[i2 + 1]
				a[i2] = tr
				a[i2 + 1] = ti
			}
			value_re = -value_re
			value_im = -value_im
		}
		pre := a[2*(col*n + col)]
		pim := a[2*(col*n + col) + 1]
		nr := value_re*pre - value_im*pim
		ni := value_re*pim + value_im*pre
		value_re = nr
		value_im = ni
		den := pre*pre + pim*pim
		inv_re := pre / den
		inv_im := -pim / den
		for r in col + 1 .. n {
			fr := a[2*(r*n + col)]
			fi := a[2*(r*n + col) + 1]
			if fr == 0.0 && fi == 0.0 {
				continue
			}
			gfr := fr*inv_re - fi*inv_im
			gfi := fr*inv_im + fi*inv_re
			for c in col .. n {
				ar := a[2*(col*n + c)]
				ai := a[2*(col*n + c) + 1]
				a[2*(r*n + c)] -= gfr*ar - gfi*ai
				a[2*(r*n + c) + 1] -= gfr*ai + gfi*ar
			}
		}
	}
	return [value_re, value_im]
}

// matrix_inverse returns the inverse of a complex `n`-by-`n` matrix
// (interleaved) by Gauss-Jordan elimination.
@[vc_fn]
pub fn matrix_inverse(matrix []f64, n int) []f64 {
	mut a := matrix.clone()
	mut m := complex_identity_flat(n)
	for col in 0 .. n {
		mut pivot := col
		mut best := math.sqrt(a[2*(col*n + col)] * a[2*(col*n + col)] + a[2*(col*n + col) + 1] * a[2*(col*n + col) + 1])
		for r in col + 1 .. n {
			mag := math.sqrt(a[2*(r*n + col)] * a[2*(r*n + col)] + a[2*(r*n + col) + 1] * a[2*(r*n + col) + 1])
			if mag > best {
				best = mag
				pivot = r
			}
		}
		if best < 1e-14 {
			return []f64{}
		}
		if pivot != col {
			for c in 0 .. n {
				i1 := 2 * (col*n + c)
				i2 := 2 * (pivot*n + c)
				tr := a[i1]
				ti := a[i1 + 1]
				a[i1] = a[i2]
				a[i1 + 1] = a[i2 + 1]
				a[i2] = tr
				a[i2 + 1] = ti
				j1 := 2 * (col*n + c)
				j2 := 2 * (pivot*n + c)
				tr2 := m[j1]
				ti2 := m[j1 + 1]
				m[j1] = m[j2]
				m[j1 + 1] = m[j2 + 1]
				m[j2] = tr2
				m[j2 + 1] = ti2
			}
		}
		pre := a[2*(col*n + col)]
		pim := a[2*(col*n + col) + 1]
		den := pre*pre + pim*pim
		inv_re := pre / den
		inv_im := -pim / den
		for c in 0 .. n {
			ar := a[2*(col*n + c)]
			ai := a[2*(col*n + c) + 1]
			a[2*(col*n + c)] = ar*inv_re - ai*inv_im
			a[2*(col*n + c) + 1] = ar*inv_im + ai*inv_re
			mr := m[2*(col*n + c)]
			mi := m[2*(col*n + c) + 1]
			m[2*(col*n + c)] = mr*inv_re - mi*inv_im
			m[2*(col*n + c) + 1] = mr*inv_im + mi*inv_re
		}
		for r in 0 .. n {
			if r == col {
				continue
			}
			fr := a[2*(r*n + col)]
			fi := a[2*(r*n + col) + 1]
			if fr == 0.0 && fi == 0.0 {
				continue
			}
			for c in 0 .. n {
				ar := a[2*(col*n + c)]
				ai := a[2*(col*n + c) + 1]
				a[2*(r*n + c)] -= fr*ar - fi*ai
				a[2*(r*n + c) + 1] -= fr*ai + fi*ar
				mr := m[2*(col*n + c)]
				mi := m[2*(col*n + c) + 1]
				m[2*(r*n + c)] -= fr*mr - fi*mi
				m[2*(r*n + c) + 1] -= fr*mi + fi*mr
			}
		}
	}
	return m
}

fn complex_identity_flat(n int) []f64 {
	mut out := []f64{len: 2 * n * n, init: 0.0}
	for i in 0 .. n {
		out[2*(i*n + i)] = 1.0
	}
	return out
}

// solve3 solves a real 3-by-3 system `matrix @ x = rhs`.
@[vc_fn]
pub fn solve3(matrix []f64, rhs []f64) []f64 {
	mut a := []f64{len: 12, init: 0.0}
	for i in 0 .. 3 {
		for j in 0 .. 3 {
			a[i*4 + j] = matrix[i*3 + j]
		}
		a[i*4 + 3] = rhs[i]
	}
	for col in 0 .. 3 {
		mut pivot := col
		for r in col + 1 .. 3 {
			if math.abs(a[r*4 + col]) > math.abs(a[pivot*4 + col]) {
				pivot = r
			}
		}
		if pivot != col {
			for c in 0 .. 4 {
				tmp := a[col*4 + c]
				a[col*4 + c] = a[pivot*4 + c]
				a[pivot*4 + c] = tmp
			}
		}
		inv := 1.0 / a[col*4 + col]
		for c in 0 .. 4 {
			a[col*4 + c] *= inv
		}
		for r in 0 .. 3 {
			if r == col {
				continue
			}
			factor := a[r*4 + col]
			for c in 0 .. 4 {
				a[r*4 + c] -= factor * a[col*4 + c]
			}
		}
	}
	return [a[3], a[7], a[11]]
}

// kron_factor splits a 4x4 matrix that is `g * kron(a, b)` into `g`, `a` and
// `b`. It returns `[g_re, g_im, a (8), b (8)]`.
@[vc_fn]
pub fn kron_factor(matrix []f64) []f64 {
	mut best := 0
	mut best_mag := -1.0
	for i in 0 .. 4 {
		for j in 0 .. 4 {
			mag := matrix[2*(i*4 + j)] * matrix[2*(i*4 + j)] + matrix[2*(i*4 + j) + 1] * matrix[2*(i*4 + j) + 1]
			if mag > best_mag {
				best_mag = mag
				best = i * 4 + j
			}
		}
	}
	ar := best / 4
	br := best % 4
	mut f1 := []f64{len: 8, init: 0.0}
	mut f2 := []f64{len: 8, init: 0.0}
	for i in 0 .. 2 {
		for j in 0 .. 2 {
			p1 := ((ar >> 1) ^ i) * 2 + ((br >> 1) ^ j)
			p2 := ((ar & 1) ^ i) * 2 + ((br & 1) ^ j)
			s1 := (ar ^ (i << 1)) * 4 + (br ^ (j << 1))
			s2 := (ar ^ i) * 4 + (br ^ j)
			f1[2*p1] = matrix[2*s1]
			f1[2*p1 + 1] = matrix[2*s1 + 1]
			f2[2*p2] = matrix[2*s2]
			f2[2*p2 + 1] = matrix[2*s2 + 1]
		}
	}
	f1 = normalize_det2(f1)
	f2 = normalize_det2(f2)
	// global phase g = matrix[a][b] / (f1[..] * f2[..])
	ai := ar >> 1
	bi := br >> 1
	c1 := 2 * (ai*2 + bi)
	c2 := 2 * ((ar & 1)*2 + (br & 1))
	den_re := f1[c1] * f2[c2] - f1[c1 + 1] * f2[c2 + 1]
	den_im := f1[c1] * f2[c2 + 1] + f1[c1 + 1] * f2[c2]
	m_re := matrix[2*best]
	m_im := matrix[2*best + 1]
	den := den_re*den_re + den_im*den_im
	mut g_re := (m_re*den_re + m_im*den_im) / den
	mut g_im := (m_im*den_re - m_re*den_im) / den
	if g_re < 0 {
		for k in 0 .. 8 {
			f1[k] = -f1[k]
		}
		g_re = -g_re
		g_im = -g_im
	}
	mut out := []f64{len: 18, init: 0.0}
	out[0] = g_re
	out[1] = g_im
	for k in 0 .. 8 {
		out[2 + k] = f1[k]
		out[10 + k] = f2[k]
	}
	return out
}

fn normalize_det2(factor []f64) []f64 {
	mut out := factor.clone()
	determinant := matrix_det(out, 2)
	// complex square root of the determinant
	re := determinant[0]
	im := determinant[1]
	modulus := math.sqrt(re*re + im*im)
	if modulus < 1e-14 {
		return out
	}
	root_re := math.sqrt((modulus + re) / 2)
	mut root_im := math.sqrt((modulus - re) / 2)
	if im < 0 {
		root_im = -root_im
	}
	den := root_re*root_re + root_im*root_im
	for k in 0 .. 4 {
		vr := out[2*k]
		vi := out[2*k + 1]
		out[2*k] = (vr*root_re + vi*root_im) / den
		out[2*k + 1] = (vi*root_re - vr*root_im) / den
	}
	return out
}
