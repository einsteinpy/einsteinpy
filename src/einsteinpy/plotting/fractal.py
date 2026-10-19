"""Easter Egg.

This module was added before ``0.3.0`` release as an Easter Egg.
Contains a simple implementation of fractal.

"""
import numpy as np
from plotly import graph_objects as go

from einsteinpy.ijit import jit


@jit
def _julia(A, c, zabs_max, i, j, dims, x_range, y_range, iter_max):
    width, height = dims
    xmin, xmax = x_range
    ymin, ymax = y_range
    xwidth, ywidth = xmax - xmin, ymax - ymin
    it = 0
    z = complex((i / width) * xwidth + xmin + 1j * ((j / height) * ywidth + ymin))
    while abs(z) < zabs_max and it < iter_max:
        z = z**2 + c
        it += 1
    A[i, j] = (it / iter_max) ** 0.2


def _fractal_img(
    dims=(384, 384),
    angle=2.75,
    zabs_max=20,
    iter_max=64,
    x_range=(-1.6, 1.6),
    y_range=(-1.6, 1.6),
):
    width, height = dims
    A = np.zeros((width, height), dtype=float)
    c = 0.80 * np.exp(1j * angle)
    for i in range(width):
        for j in range(height):
            _julia(A, c, zabs_max, i, j, dims, x_range, y_range, iter_max)
    return A


def fractal(start_angle=(np.pi / 3), end_angle=np.pi, divs=12, show=True):
    """
    It's an Easter Egg.

    This function plots fractals arising out of a specific class of functions
    from Julia set. The plot also provides a slider to visualize various angles
    which can act as parameters for the complex function.

    Parameters
    ----------
    start_angle : float
        Starting angle for the slider. Defaults to ``pi/3``.
    end_angle : float
        Ending angle for the slider. Defaults to ``pi``.
    divs : int
        Number of divisions between ``start_angle`` and ``end_angle``.
        Defaults to ``12``.
    show : bool
        Shows plot if ``True``. Defaults to ``True``.

    Returns
    -------
    ~plotly.graph_objects.Figure
        Instance of plot created.

    """
    active = divs - 1
    fig = go.Figure()
    # Add traces, one for each slider step
    for step in np.linspace(start_angle, end_angle, num=divs):
        fig.add_trace(
            go.Heatmap(z=_fractal_img(angle=step, dims=(400, 400)), showscale=False)
        )

    fig.data[active].visible = True

    text = """
    A fractal is a never ending pattern.
    We are generating fractals out of Julia set with different angles as parameters.

    """
    print(text)

    steps = list()
    for i in range(len(fig.data)):
        step = dict(method="restyle", args=["visible", [False] * len(fig.data)])
        step["args"][1][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)

    sliders = [
        dict(
            active=active,
            currentvalue={
                "prefix": "Angle between %.2f & %.2f : " % (start_angle, end_angle)
            },
            pad={"t": divs},
            steps=steps,
        )
    ]

    fig.update_layout(
        sliders=sliders,
        width=600,
        height=600,
        xaxis=dict(
            autorange=True,
            showgrid=False,
            ticks="",
            showticklabels=False,
            fixedrange=True,
        ),
        yaxis=dict(
            autorange=True,
            showgrid=False,
            ticks="",
            showticklabels=False,
            fixedrange=True,
        ),
    )

    text = """
    Waiting for plotly object to load.

    """

    if show:
        print(text)
        fig.show()

    return fig
