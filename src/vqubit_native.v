module vqubit_native

// version returns the version of the native extension.
//
// The Python package re-exports it as `vqubit.__version__`. The value is kept
// in step with `vcraft.toml` by a test, since the two are written by hand.
@[vc_fn]
pub fn version() string {
	return '0.2.0'
}
