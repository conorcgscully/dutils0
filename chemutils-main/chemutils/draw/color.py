import math

import numpy as np
import webcolors

Color = str


def parse_color(color: Color, /) -> tuple[float, float, float]:
    """Parse a color specification into a float tuple."""
    if color[0] == "#":
        r, g, b = webcolors.hex_to_rgb(color)
    else:
        r, g, b = webcolors.name_to_rgb(color)
    return (r / 255.0, g / 255.0, b / 255.0)


gradient = np.array(
    [
        webcolors.name_to_rgb("lightskyblue"),
        webcolors.name_to_rgb("turquoise"),
        webcolors.name_to_rgb("mediumspringgreen"),
        webcolors.name_to_rgb("yellowgreen"),
        webcolors.name_to_rgb("gold"),
        webcolors.name_to_rgb("darkorange"),
        webcolors.name_to_rgb("salmon"),
        webcolors.name_to_rgb("hotpink"),
        webcolors.name_to_rgb("violet"),
    ]
)
"""
Hand chosen gradient based on HTML color names.

Chosen to be of approximately equal lightness, and evenly
spaced in CAM02-UCS space.
"""


def get_colors(n: int) -> list[str]:
    """
    Get the list of N colors chosen from a predefined gradient.

    This is set up so when n is 1, it returns the central color
    of `gold` (the default for highlights).
    """
    width = min((n - 1) * 0.25, 1)
    points = np.linspace(0.5 - width / 2, 0.5 + width / 2, n)
    return [get_color_gradient(point) for point in points]


def get_color_gradient(fraction: float) -> str:
    """Get the color that corresponds to a fraction of the gradient."""
    num_viridis_points = len(gradient)
    fract = fraction * (num_viridis_points - 1)
    lower_index = math.floor(fract)
    upper_index = math.ceil(fract)
    if lower_index == upper_index:
        color = gradient[lower_index]
    else:
        fract_lower = upper_index - fract
        fract_upper = fract - lower_index
        color = [
            gradient[lower_index][0] * fract_lower + gradient[upper_index][0] * fract_upper,
            gradient[lower_index][1] * fract_lower + gradient[upper_index][1] * fract_upper,
            gradient[lower_index][2] * fract_lower + gradient[upper_index][2] * fract_upper,
        ]
    return webcolors.rgb_to_hex((int(color[0]), int(color[1]), int(color[2])))  # type: ignore


# 0, 1, 2, 3, 4, 5  (6 points)
