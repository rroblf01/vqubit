"""Visualization helpers, mirroring the ``cirq.vis`` tools.

These functions and classes need Matplotlib (and therefore NumPy), which are
*not* dependencies of vqubit, so they import them lazily on first use. Importing
vqubit stays dependency-free.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from dataclasses import astuple, dataclass
from typing import Any


def _require_matplotlib() -> Any:
	try:
		import matplotlib.pyplot as plt
	except ImportError as exc:  # pragma: no cover - exercised only without matplotlib
		raise ImportError(
			"vqubit's visualization helpers need matplotlib; install it with 'pip install matplotlib'."
		) from exc
	return plt


def _require_numpy() -> Any:
	try:
		import numpy as np
	except ImportError as exc:  # pragma: no cover - exercised only without numpy
		raise ImportError("vqubit's visualization helpers need numpy; install it with 'pip install numpy'.") from exc
	return np


def relative_luminance(color: Any) -> float:
	"""Return the relative luminance of ``color`` per the W3C specification."""
	np = _require_numpy()
	rgb = np.asarray(color)[:3]
	rgb = np.where(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
	return float(rgb.dot([0.2126, 0.7152, 0.0722]))


# --- state histograms --------------------------------------------------------


def get_state_histogram(result: Any) -> Any:
	"""Compute a state histogram (a numpy array) from a result with repetitions."""
	np = _require_numpy()
	measurements = result.measurements
	num_qubits = sum(len(row[0]) for row in measurements.values()) if measurements else 0
	values = np.zeros(2**num_qubits)
	rows = [list(row) for row in measurements.values()]
	repetitions = len(rows[0]) if rows else 0
	for rep in range(repetitions):
		bits: list[int] = []
		for row in rows:
			bits.extend(row[rep])
		values[int("".join(str(b) for b in bits), 2) if bits else 0] += 1
	return values


def plot_state_histogram(
	data: Any,
	ax: Any = None,
	*,
	tick_label: Sequence[str] | None = None,
	xlabel: str | None = "qubit state",
	ylabel: str | None = "result count",
	title: str | None = "Result State Histogram",
) -> Any:
	"""Plot a state histogram from a result, a ``Counter`` or a sequence of values."""
	import collections

	plt = _require_matplotlib()
	np = _require_numpy()
	show_fig = not ax
	if ax is None:
		fig, ax = plt.subplots(1, 1)
	if hasattr(data, "measurements"):
		values = get_state_histogram(data)
	elif isinstance(data, collections.Counter):
		tick_label, counts = zip(*sorted(data.items()), strict=False)
		values = np.asarray(counts)
	else:
		values = np.array(data)
	if tick_label is None:
		tick_label = [str(i) for i in range(len(values))]
	ax.bar(np.arange(len(values)), values, tick_label=tick_label)
	if xlabel:
		ax.set_xlabel(xlabel)
	if ylabel:
		ax.set_ylabel(ylabel)
	if title:
		ax.set_title(title)
	if show_fig:
		fig.show()
	return ax


def integrated_histogram(
	data: Any,
	ax: Any = None,
	*,
	cdf_on_x: bool = False,
	axis_label: str = "",
	semilog: bool = True,
	median_line: bool = True,
	median_label: str | None = "median",
	mean_line: bool = False,
	mean_label: str | None = "mean",
	show_zero: bool = False,
	title: str | None = None,
	**kwargs: Any,
) -> Any:
	"""Plot the integrated histogram (cumulative distribution) of ``data``."""
	plt = _require_matplotlib()
	np = _require_numpy()
	show_plot = not ax
	if ax is None:
		fig, ax = plt.subplots(1, 1)

	if isinstance(data, Mapping):
		data = list(data.values())
	float_data = [float(d) for d in data if not np.isnan(float(d))]
	n = len(float_data)

	if not show_zero:
		bin_values = np.linspace(0, 1, n + 1)
		parameter_values = sorted(np.concatenate((np.array([0]), np.array(float_data))))
	else:
		bin_values = np.linspace(0, 1, n)
		parameter_values = sorted(float_data)
	plot_options = {"where": "post", "color": "b", "linestyle": "-", "lw": 1.0, "ms": 0.0}
	plot_options.update(kwargs)

	if cdf_on_x:
		ax.step(bin_values, parameter_values, **plot_options)
	else:
		ax.step(parameter_values, bin_values, **plot_options)

	set_semilog = ax.semilogy if cdf_on_x else ax.semilogx
	set_lim = ax.set_xlim if cdf_on_x else ax.set_ylim
	set_ticks = ax.set_xticks if cdf_on_x else ax.set_yticks
	set_line = ax.axhline if cdf_on_x else ax.axvline
	cdf_label = ax.set_xlabel if cdf_on_x else ax.set_ylabel
	ax_label = ax.set_ylabel if cdf_on_x else ax.set_xlabel

	if not title:
		title = f"N={n}"
	ax.set_title(title)
	if semilog:
		set_semilog()
	set_lim(0, 1)
	set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
	ax.grid(True)
	cdf_label("Integrated histogram")
	if axis_label:
		ax_label(axis_label)
	if "label" in plot_options:
		ax.legend()
	if median_line:
		set_line(
			float(np.median(float_data)),
			linestyle="--",
			color=plot_options["color"],
			alpha=0.5,
			label=median_label,
		)
	if mean_line:
		set_line(
			float(np.mean(float_data)),
			linestyle="-.",
			color=plot_options["color"],
			alpha=0.5,
			label=mean_label,
		)
	if show_plot:
		fig.show()
	return ax


# --- heatmaps ----------------------------------------------------------------


@dataclass
class Point:
	x: float
	y: float

	def __iter__(self) -> Any:
		return iter(astuple(self))


@dataclass
class PolygonUnit:
	"""Information about a single polygon unit of a heatmap."""

	polygon: Any
	value: float
	center: Point
	annot: str | None


class Heatmap:
	"""Distribution of a value on a 2D qubit lattice as a color map."""

	def __init__(self, value_map: Mapping[Any, Any], **kwargs: Any) -> None:
		self._value_map: dict[Any, Any] = {k if isinstance(k, tuple) else (k,): v for k, v in value_map.items()}
		self._validate_kwargs(kwargs)
		if "_config" not in self.__dict__:
			self._config: dict[str, Any] = {}
		self._config.update(
			{
				"plot_colorbar": True,
				"colorbar_position": "right",
				"colorbar_size": "5%",
				"colorbar_pad": "2%",
				"collection_options": {"cmap": "viridis"},
				"annotation_format": ".2g",
			}
		)
		self._config.update(kwargs)

	def _extra_valid_kwargs(self) -> list[str]:
		return []

	def _validate_kwargs(self, kwargs: Mapping[str, Any]) -> None:
		valid_kwargs = [
			"plot_colorbar",
			"colorbar_position",
			"colorbar_size",
			"colorbar_pad",
			"colorbar_options",
			"collection_options",
			"vmin",
			"vmax",
			"title",
			"annotation_map",
			"annotation_text_kwargs",
			"annotation_format",
			"highlighted_qubits",
			*self._extra_valid_kwargs(),
		]
		if any(k not in valid_kwargs for k in kwargs):
			invalid_args = ", ".join(k for k in kwargs if k not in valid_kwargs)
			raise ValueError(f"Received invalid argument(s): {invalid_args}")

	def update_config(self, **kwargs: Any) -> Heatmap:
		"""Add or modify the configuration arguments."""
		self._validate_kwargs(kwargs)
		self._config.update(kwargs)
		return self

	def _qubits_to_polygon(self, qubits: Any) -> tuple[Any, Point]:
		qubit = qubits[0]
		x, y = float(qubit.row), float(qubit.col)
		return (
			[(y - 0.5, x - 0.5), (y - 0.5, x + 0.5), (y + 0.5, x + 0.5), (y + 0.5, x - 0.5)],
			Point(y, x),
		)

	def _get_annotation_value(self, key: Any, value: Any) -> str | None:
		if self._config.get("annotation_map"):
			return self._config["annotation_map"].get(key)
		if self._config.get("annotation_format"):
			try:
				return format(value, self._config["annotation_format"])
			except (TypeError, ValueError):
				return format(float(value), self._config["annotation_format"])
		return None

	def _get_polygon_units(self) -> list[PolygonUnit]:
		units: list[PolygonUnit] = []
		for qubits, value in sorted(self._value_map.items()):
			polygon, center = self._qubits_to_polygon(qubits)
			units.append(
				PolygonUnit(
					polygon=polygon,
					center=center,
					value=float(value),
					annot=self._get_annotation_value(qubits, value),
				)
			)
		return units

	def _plot_colorbar(self, mappable: Any, ax: Any) -> Any:
		from mpl_toolkits import axes_grid1

		colorbar_ax = axes_grid1.make_axes_locatable(ax).append_axes(
			position=self._config["colorbar_position"],
			size=self._config["colorbar_size"],
			pad=self._config["colorbar_pad"],
		)
		position = self._config["colorbar_position"]
		orientation = "vertical" if position in ("left", "right") else "horizontal"
		colorbar = ax.figure.colorbar(
			mappable, colorbar_ax, ax, orientation=orientation, **self._config.get("colorbar_options", {})
		)
		colorbar_ax.tick_params(axis="y", direction="out")
		return colorbar

	def _write_annotations(self, centers_and_annot: list[Any], collection: Any, ax: Any) -> None:
		face_colors = collection.get_facecolor()
		for (center, annotation), facecolor in zip(centers_and_annot, face_colors, strict=False):
			if not annotation:
				continue
			x, y = center
			text_color = "black" if relative_luminance(facecolor) > 0.4 else "white"
			text_kwargs: dict[str, Any] = {"color": text_color, "ha": "center", "va": "center"}
			text_kwargs.update(self._config.get("annotation_text_kwargs", {}))
			ax.text(x, y, annotation, **text_kwargs)

	def _plot_on_axis(self, ax: Any) -> Any:
		import matplotlib.collections as mpl_collections

		np = _require_numpy()
		polygon_list = self._get_polygon_units()
		collection = mpl_collections.PolyCollection(
			[c.polygon for c in polygon_list], **self._config.get("collection_options", {})
		)
		collection.set_clim(self._config.get("vmin"), self._config.get("vmax"))
		collection.set_array(np.array([c.value for c in polygon_list]))
		ax.add_collection(collection)
		collection.update_scalarmappable()
		if self._config.get("annotation_map") or self._config.get("annotation_format"):
			self._write_annotations([(c.center, c.annot) for c in polygon_list], collection, ax)
		ax.set(xlabel="column", ylabel="row")
		if self._config.get("plot_colorbar"):
			self._plot_colorbar(collection, ax)
		rows = {q.row for qubits in self._value_map for q in qubits}
		cols = {q.col for qubits in self._value_map for q in qubits}
		min_row, max_row = min(rows), max(rows)
		min_col, max_col = min(cols), max(cols)
		min_xtick, max_xtick = np.floor(min_col), np.ceil(max_col)
		ax.set_xticks(np.arange(min_xtick, max_xtick + 1))
		min_ytick, max_ytick = np.floor(min_row), np.ceil(max_row)
		ax.set_yticks(np.arange(min_ytick, max_ytick + 1))
		ax.set_xlim((min_xtick - 0.6, max_xtick + 0.6))
		ax.set_ylim((max_ytick + 0.6, min_ytick - 0.6))
		if self._config.get("title"):
			ax.set_title(self._config["title"], fontweight="bold")
		return collection

	def plot(self, ax: Any = None, **kwargs: Any) -> tuple[Any, Any]:
		"""Plot the heatmap on ``ax`` (creating a figure when it is ``None``)."""
		plt = _require_matplotlib()
		show_plot = not ax
		if ax is None:
			fig, ax = plt.subplots(figsize=(8, 8))
		original_config = copy.deepcopy(self._config)
		self.update_config(**kwargs)

		highlighted_qubits = frozenset(kwargs.get("highlighted_qubits", ()))
		if highlighted_qubits:
			sorted_keys = sorted(self._value_map.keys())
			opts = self._config["collection_options"]
			opts.update(
				{
					"edgecolors": tuple(
						"red" if not highlighted_qubits.isdisjoint(qs) else opts.get("edgecolors", "grey")
						for qs in sorted_keys
					),
					"linestyles": tuple(
						"solid" if not highlighted_qubits.isdisjoint(qs) else opts.get("linestyles", "dashed")
						for qs in sorted_keys
					),
					"linewidths": tuple(
						4 if not highlighted_qubits.isdisjoint(qs) else opts.get("linewidths", 2) for qs in sorted_keys
					),
				}
			)

		collection = self._plot_on_axis(ax)
		if show_plot:
			fig.show()
		self._config = original_config
		return (ax, collection)


class TwoQubitInteractionHeatmap(Heatmap):
	"""Visualize interactions between neighboring qubits on a 2D grid."""

	def __init__(self, value_map: Mapping[Any, Any], **kwargs: Any) -> None:
		self._config: dict[str, Any] = {"coupler_margin": 0.03, "coupler_width": 0.6}
		super().__init__(value_map, **kwargs)

	def _extra_valid_kwargs(self) -> list[str]:
		return ["coupler_margin", "coupler_width"]

	def _qubits_to_polygon(self, qubits: Any) -> tuple[Any, Point]:
		coupler_margin = self._config["coupler_margin"]
		coupler_width = self._config["coupler_width"]
		cwidth = coupler_width / 2.0
		setback = 0.5 - cwidth
		row1, col1 = map(float, (qubits[0].row, qubits[0].col))
		row2, col2 = map(float, (qubits[1].row, qubits[1].col))
		if abs(row1 - row2) + abs(col1 - col2) != 1:
			raise ValueError(f"{qubits[0]}-{qubits[1]} is not supported because they are not nearest neighbors")
		if coupler_width <= 0:
			polygon: list[tuple[float, float]] = []
		elif row1 == row2:  # horizontal
			col1, col2 = min(col1, col2), max(col1, col2)
			col_center = (col1 + col2) / 2.0
			polygon = [
				(col1 + coupler_margin, row1),
				(col_center - setback, row1 + cwidth - coupler_margin),
				(col_center + setback, row1 + cwidth - coupler_margin),
				(col2 - coupler_margin, row2),
				(col_center + setback, row1 - cwidth + coupler_margin),
				(col_center - setback, row1 - cwidth + coupler_margin),
			]
		else:  # vertical
			row1, row2 = min(row1, row2), max(row1, row2)
			row_center = (row1 + row2) / 2.0
			polygon = [
				(col1, row1 + coupler_margin),
				(col1 + cwidth - coupler_margin, row_center - setback),
				(col1 + cwidth - coupler_margin, row_center + setback),
				(col2, row2 - coupler_margin),
				(col1 - cwidth + coupler_margin, row_center + setback),
				(col1 - cwidth + coupler_margin, row_center - setback),
			]
		return (polygon, Point((col1 + col2) / 2.0, (row1 + row2) / 2.0))

	def plot(self, ax: Any = None, **kwargs: Any) -> tuple[Any, Any]:
		"""Plot the interaction heatmap on ``ax``."""
		plt = _require_matplotlib()
		show_plot = not ax
		if ax is None:
			fig, ax = plt.subplots(figsize=(8, 8))
		original_config = copy.deepcopy(self._config)
		self.update_config(**kwargs)
		qubits = {q for qubits in self._value_map for q in qubits}
		collection_options: dict[str, Any] = {"cmap": "binary"}
		highlighted_qubits = frozenset(kwargs.get("highlighted_qubits", ()))
		if not highlighted_qubits:
			collection_options.update({"linewidths": 2, "edgecolors": "lightgrey", "linestyles": "dashed"})
		Heatmap(dict.fromkeys(qubits, 0.0)).plot(
			ax=ax,
			collection_options=collection_options,
			plot_colorbar=False,
			annotation_format=None,
			highlighted_qubits=highlighted_qubits,
		)
		collection = self._plot_on_axis(ax)
		if show_plot:
			fig.show()
		self._config = original_config
		return (ax, collection)


# --- density matrix ----------------------------------------------------------


def _plot_element_of_density_matrix(
	ax: Any, x: int, y: int, r: float, phase: float, show_rect: bool = False, show_text: bool = False
) -> None:
	from matplotlib import lines, patches

	plt = _require_matplotlib()
	np = _require_numpy()
	half = (1 / 1.1) * 0.5
	margin = 0.01
	opacity = 0.8 if not show_text else 0.4
	circle_out = plt.Circle((x + 0.5, y + 0.5), radius=half, fill=False, color="#333333")
	circle_in = plt.Circle((x + 0.5, y + 0.5), radius=r * half, fill=True, color="IndianRed", alpha=opacity)
	line = lines.Line2D(
		(x + 0.5, x + 0.5 + np.cos(phase) * half),
		(y + 0.5, y + 0.5 + np.sin(phase) * half),
		color="#333333",
		alpha=opacity,
	)
	ax.add_artist(circle_in)
	ax.add_artist(circle_out)
	ax.add_artist(line)
	if show_rect:
		rect = patches.Rectangle((x + margin, y + margin), 1.0 - 2 * margin, r * (1 - 2 * margin), alpha=0.25)
		ax.add_artist(rect)
	if show_text:
		plt.text(
			x + 0.5,
			y + 0.5,
			f"{np.round(r, decimals=2)}\n{np.round(phase * 180 / np.pi, decimals=2)} deg",
			horizontalalignment="center",
			verticalalignment="center",
		)


def plot_density_matrix(matrix: Any, ax: Any = None, *, show_text: bool = False, title: str | None = None) -> Any:
	"""Plot the magnitudes and phases of a density matrix (quirk-style)."""
	plt = _require_matplotlib()
	np = _require_numpy()
	plt.style.use("ggplot")
	padding = 0.001
	matrix = matrix.astype(np.complex128)
	num_qubits = int(np.log2(matrix.shape[0]))
	if ax is None:
		_, ax = plt.subplots(figsize=(10, 10))
	ax.set_xlim(0 - padding, 2**num_qubits + padding)
	ax.set_ylim(0 - padding, 2**num_qubits + padding)
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			_plot_element_of_density_matrix(
				ax,
				i,
				j,
				float(np.abs(matrix[i][-j - 1])),
				float(np.angle(matrix[i][-j - 1])),
				show_rect=(i == matrix.shape[1] - j - 1),
				show_text=show_text,
			)
	ticks = np.arange(0.5, matrix.shape[0])
	labels = [f"{'0' * (num_qubits - len(f'{i:b}'))}{i:b}" for i in range(matrix.shape[0])]
	ax.set_xticks(ticks)
	ax.set_xticklabels(labels, rotation=90)
	ax.set_yticks(ticks)
	ax.set_yticklabels(list(reversed(labels)))
	ax.set_facecolor("#eeeeee")
	if title is not None:
		ax.set_title(title)
	return ax


# --- KAK interaction scatter -------------------------------------------------


def scatter_plot_normalized_kak_interaction_coefficients(
	interactions: Any,
	*,
	include_frame: bool = True,
	ax: Any = None,
	**kwargs: Any,
) -> Any:
	"""Scatter-plot the normalized KAK interaction coefficients of two-qubit gates."""
	from ._linalg_public import kak_vector
	from ._protocols import unitary

	plt = _require_matplotlib()
	np = _require_numpy()
	show_plot = not ax
	if ax is None:
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1, projection="3d")

	def coord_transform(pts: Any) -> tuple[Any, Any, Any]:
		if len(pts) == 0:
			return np.array([]), np.array([]), np.array([])
		xs, ys, zs = np.transpose(pts)
		return xs, zs, ys

	if include_frame:
		envelope = [
			(0, 0, 0),
			(1, 1, 1),
			(1, 1, -1),
			(0, 0, 0),
			(1, 1, 1),
			(1, 0, 0),
			(0, 0, 0),
			(1, 1, -1),
			(1, 0, 0),
			(0, 0, 0),
			(1, 0, 0),
			(1, 1, 0),
			(0, 0, 0),
		]
		ax.plot(*coord_transform(envelope), c="black", linewidth=1)

	if not isinstance(interactions, np.ndarray):
		extracted = [a if isinstance(a, np.ndarray) else unitary(a) for a in interactions]
	else:
		extracted = [interactions]
	points = kak_vector(extracted) * 4 / np.pi
	ax.scatter(*coord_transform(points), **kwargs)
	ax.set_xlim(0, +1)
	ax.set_ylim(-1, +1)
	ax.set_zlim(0, +1)
	if show_plot:
		fig.show()
	return ax


__all__ = [
	"Heatmap",
	"Point",
	"PolygonUnit",
	"TwoQubitInteractionHeatmap",
	"get_state_histogram",
	"integrated_histogram",
	"plot_density_matrix",
	"plot_state_histogram",
	"relative_luminance",
	"scatter_plot_normalized_kak_interaction_coefficients",
]
