"""Tests for the plotting helpers (``vqubit.vis``).

They are skipped when matplotlib is not installed, since it is an optional
extra; when it is, each helper has to return the axes it drew on, so a caller can
compose figures.
"""

from __future__ import annotations

import collections
from collections.abc import Iterator

import pytest

import vqubit

matplotlib = pytest.importorskip("matplotlib", exc_type=ImportError)
matplotlib.use("Agg")  # a headless backend: nothing is displayed
numpy = pytest.importorskip("numpy", exc_type=ImportError)

from matplotlib import pyplot as plt  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
	yield
	plt.close("all")


def measured_result() -> vqubit.Result:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1), vqubit.measure(q0, q1, key="m"))
	return vqubit.Simulator(seed=0).run(circuit, repetitions=50)


def test_relative_luminance() -> None:
	assert vqubit.relative_luminance([0.0, 0.0, 0.0]) == pytest.approx(0.0)
	assert vqubit.relative_luminance([1.0, 1.0, 1.0]) == pytest.approx(1.0)
	assert 0.0 < vqubit.relative_luminance([0.5, 0.5, 0.5]) < 1.0


def test_get_state_histogram() -> None:
	histogram = vqubit.get_state_histogram(measured_result())
	assert float(histogram.sum()) == 50
	assert histogram.shape == (4,)
	# An empty result has an empty histogram.
	assert vqubit.get_state_histogram(vqubit.Result({})).shape == (1,)


def test_plot_state_histogram_from_a_result() -> None:
	ax = vqubit.plot_state_histogram(measured_result())
	assert ax.get_xlabel() == "qubit state"
	assert ax.get_ylabel() == "result count"
	assert ax.get_title() == "Result State Histogram"


def test_plot_state_histogram_from_a_counter_and_a_sequence() -> None:
	ax = vqubit.plot_state_histogram(collections.Counter({0: 3, 1: 5}))
	assert [label.get_text() for label in ax.get_xticklabels()] == ["0", "1"]
	# A plain sequence of counts, with no labels or titles.
	ax = vqubit.plot_state_histogram([1, 2, 3], tick_label=["a", "b", "c"], xlabel=None, ylabel=None, title=None)
	assert ax.get_xlabel() == ""
	assert ax.get_ylabel() == ""
	assert ax.get_title() == ""


def test_plot_state_histogram_on_a_given_axis() -> None:
	_figure, ax = plt.subplots()
	assert vqubit.plot_state_histogram([1, 2], ax=ax, title=None) is ax
	assert ax.get_title() == ""


def test_integrated_histogram() -> None:
	ax = vqubit.integrated_histogram([0.1, 0.2, 0.3, 0.4])
	assert ax.get_title() == "N=4"
	assert ax.get_ylabel() == "Integrated histogram"
	# The cumulative distribution can go on the x axis, with a log scale, and can
	# draw the mean and the median.
	ax = vqubit.integrated_histogram(
		{0: 0.2, 1: 0.4, 2: 0.6}, cdf_on_x=True, axis_label="fidelity", show_zero=True, mean_line=True
	)
	assert ax.get_xlabel() == "Integrated histogram"
	assert ax.get_ylabel() == "fidelity"
	ax = vqubit.integrated_histogram([0.5, 0.5], semilog=False, median_line=False, title="custom")
	assert ax.get_title() == "custom"


def test_integrated_histogram_skips_nan() -> None:
	ax = vqubit.integrated_histogram([0.1, float("nan"), 0.3])
	assert ax.get_title() == "N=2"


def test_plot_density_matrix() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	circuit = vqubit.Circuit(vqubit.H(q0), vqubit.CNOT(q0, q1))
	matrix = numpy.asarray(vqubit.DensityMatrixSimulator().simulate(circuit).density_matrix())
	ax = vqubit.plot_density_matrix(matrix, title="density")
	assert ax.get_title() == "density"
	# With the numbers written in, and without a title.
	ax = vqubit.plot_density_matrix(matrix, show_text=True)
	assert ax.get_title() == ""
	assert len(ax.texts) > 0


def test_heatmap_of_a_single_qubit() -> None:
	q0 = vqubit.GridQubit(0, 0)
	heatmap = vqubit.Heatmap({q0: 0.5}, title="one", vmin=0, vmax=1)
	assert heatmap.update_config(title="two") is heatmap
	ax, collection = heatmap.plot()
	assert ax.get_title() == "two"
	assert collection is not None
	# The same plot on a given axis, with a custom annotation map.
	_figure, ax = plt.subplots()
	plotted_ax, _ = vqubit.Heatmap({q0: 0.25}, annotation_map={q0: "x"}).plot(ax=ax)
	assert plotted_ax is ax


def test_heatmap_rejects_unknown_arguments() -> None:
	q0 = vqubit.GridQubit(0, 0)
	with pytest.raises(ValueError):
		vqubit.Heatmap({q0: 1.0}, not_an_option=True)
	with pytest.raises(ValueError):
		vqubit.Heatmap({q0: 1.0}).update_config(nope=1)


def test_heatmap_annotations_and_colorbar_options() -> None:
	q0, q1 = vqubit.GridQubit(0, 0), vqubit.GridQubit(0, 1)
	heatmap = vqubit.Heatmap(
		{q0: 0.5, q1: 1.0},
		plot_colorbar=False,
		annotation_format=".1f",
		highlighted_qubits=[q0],
		collection_options={"cmap": "gray"},
	)
	ax, collection = heatmap.plot()
	assert collection is not None
	assert ax is not None


def test_two_qubit_interaction_heatmap() -> None:
	q0, q1 = vqubit.GridQubit(0, 0), vqubit.GridQubit(0, 1)
	heatmap = vqubit.TwoQubitInteractionHeatmap({(q0, q1): 0.5})
	ax, collection = heatmap.plot()
	assert collection is not None
	assert ax is not None


def test_scatter_plot_of_kak_coefficients() -> None:
	q0, q1 = vqubit.LineQubit.range(2)
	interactions = [vqubit.CZ(q0, q1), vqubit.ISWAP(q0, q1)]
	ax = vqubit.scatter_plot_normalized_kak_interaction_coefficients(interactions)
	assert ax.name == "3d"
	assert len(ax.collections) > 0
	# Without the frame, and on a given axis.
	figure = plt.figure()
	ax = figure.add_subplot(1, 1, 1, projection="3d")
	assert vqubit.scatter_plot_normalized_kak_interaction_coefficients([], include_frame=False, ax=ax) is ax
