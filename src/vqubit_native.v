module vqubit_native

import vcraft

// greet returns a greeting for `name`.
//
// The doc comment becomes the Python `__doc__`, so it is worth writing.
@[vc_fn]
pub fn greet(name string) string {
	return 'Hello, ${name}!'
}

// add returns the sum of two integers.
@[vc_fn]
pub fn add(a int, b int) int {
	return a + b
}

// parse_int reads a decimal integer, refusing anything else.
//
// The exceptions are the ones Python's own int() raises, which is the point of
// raise_domain: a bare error() would make all three a RuntimeError and an
// `except ValueError` around this call would stop matching.
@[vc_fn]
pub fn parse_int(text string) !int {
	mut value := 0
	mut seen := false
	for ch in text {
		if ch < `0` || ch > `9` {
			return vcraft.raise_domain(.value_error, 'invalid literal for int()')
		}
		seen = true
		value = value * 10 + int(ch - `0`)
	}
	if !seen {
		return vcraft.raise_domain(.value_error, 'invalid literal for int()')
	}
	return value
}

// Counter is a class. Fields are scalars and become read/write attributes;
// methods become methods; a property becomes a read-only attribute.
//
// Fields must be scalars. A V string held in a Python object is a pointer V's
// collector cannot see, because it does not scan memory CPython allocated.
@[vc_class]
pub struct Counter {
mut:
	@[vc_field] value int
	@[vc_field] step int
}

// new_counter is the constructor. It takes no arguments; an instance is made by
// calling the type, and `Counter(1)` is a TypeError.
//
//
// `@[vc_fn]` is what makes the generator see this as a constructor. Without the
// annotation it is an ordinary private helper, `link_classes` never finds it, and
// the instances come out holding whatever the allocator left behind. The
// generator then removes it from the exported set, because a `new_*` is the
// type's `tp_new` rather than something Python can call.
@[vc_fn]
pub fn new_counter() &Counter {
	return &Counter{ step: 1 }
}

// increment adds step to value and returns the new total.
@[vc_method]
pub fn (mut c Counter) increment() int {
	c.value += c.step
	return c.value
}

// is_zero reports whether the value is still zero.
@[vc_method]
@[vc_property]
pub fn (c &Counter) is_zero() bool {
	return c.value == 0
}
