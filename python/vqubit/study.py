"""The ``vqubit.study`` namespace, mirroring ``cirq.study``.

vqubit's public API is flat, so every name lives on ``vqubit`` itself;
this module only re-exports the ones Cirq groups under ``cirq.study`` so a
deep import (``from cirq.study import ...``) keeps working.
"""

from __future__ import annotations

from . import (
	UNIT_SWEEP,
	Concat,
	ExpressionMap,
	Linspace,
	ListSweep,
	ParamDictType,
	ParamMappingType,
	ParamResolver,
	ParamResolverOrSimilarType,
	Points,
	Product,
	Result,
	ResultDict,
	Sweep,
	Sweepable,
	UnitSweep,
	Zip,
	ZipLongest,
	dict_to_product_sweep,
	dict_to_zip_sweep,
	flatten,
	flatten_with_params,
	flatten_with_sweep,
	list_of_dicts_to_zip,
	symbol,
	to_resolvers,
	to_sweep,
	to_sweeps,
)

__all__ = [
	"UNIT_SWEEP",
	"Concat",
	"ExpressionMap",
	"Linspace",
	"ListSweep",
	"ParamDictType",
	"ParamMappingType",
	"ParamResolver",
	"ParamResolverOrSimilarType",
	"Points",
	"Product",
	"Result",
	"ResultDict",
	"Sweep",
	"Sweepable",
	"UnitSweep",
	"Zip",
	"ZipLongest",
	"dict_to_product_sweep",
	"dict_to_zip_sweep",
	"flatten",
	"flatten_with_params",
	"flatten_with_sweep",
	"list_of_dicts_to_zip",
	"symbol",
	"to_resolvers",
	"to_sweep",
	"to_sweeps",
]
