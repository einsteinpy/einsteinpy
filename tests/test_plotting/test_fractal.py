from io import StringIO
from unittest import mock

import plotly.graph_objects as go

from einsteinpy.plotting import fractal


@mock.patch.object(go.Figure, "show")
@mock.patch("sys.stdout", new_callable=StringIO)
def test_fractal_shows_figure(mock_stdout, mock_show):
    fig = fractal(divs=2, show=True)
    mock_show.assert_any_call()
    assert "plotly" in mock_stdout.getvalue()


def test_fractal_draws_figure():
    fig = fractal(divs=2, show=False)
    assert fig
