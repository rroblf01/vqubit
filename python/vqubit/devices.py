"""The ``vqubit.devices`` namespace, mirroring ``cirq.devices``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.devices`` so a
deep import (``from cirq.devices import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	NO_NOISE,
	NOISE_MODEL_LIKE,
	UNCONSTRAINED_DEVICE,
	ConstantQubitNoiseModel,
	Device,
	DeviceMetadata,
	GridDeviceMetadata,
	GridQid,
	GridQubit,
	LineQid,
	LineQubit,
	LineTopology,
	NamedTopology,
	NoiseModel,
	NoiseModelFromNoiseProperties,
	NoiseProperties,
	OpIdentifier,
	SuperconductingQubitsNoiseProperties,
	TiltedSquareLattice,
	draw_gridlike,
	draw_placements,
	get_placements,
	is_valid_placement,
)

__all__ = [
	"NOISE_MODEL_LIKE",
	"NO_NOISE",
	"UNCONSTRAINED_DEVICE",
	"ConstantQubitNoiseModel",
	"Device",
	"DeviceMetadata",
	"GridDeviceMetadata",
	"GridQid",
	"GridQubit",
	"LineQid",
	"LineQubit",
	"LineTopology",
	"NamedTopology",
	"NoiseModel",
	"NoiseModelFromNoiseProperties",
	"NoiseProperties",
	"OpIdentifier",
	"SuperconductingQubitsNoiseProperties",
	"TiltedSquareLattice",
	"draw_gridlike",
	"draw_placements",
	"get_placements",
	"is_valid_placement",
]
