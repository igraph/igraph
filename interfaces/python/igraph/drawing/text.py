"""
Drawers for labels on plots.
"""

__all__ = ["TextDrawer"]
__license__ = "GPL"

from igraph.drawing.baseclasses import AbstractCairoDrawer

#####################################################################

class TextDrawer(AbstractCairoDrawer):
    """Class that draws text on a Cairo context.
    
    This class supports multi-line text unlike the original Cairo text
    drawing methods."""

    LEFT, CENTER, RIGHT = "left", "center", "right"

    def __init__(self, context, text="", halign="center", valign="center"):
        """Constructs a new instance that will draw the given `text` on
        the given Cairo `context`."""
        super(TextDrawer, self).__init__(context, (0, 0))
        self.text = text
        self.halign = halign
        self.valign = valign

    def draw(self, x = None, y = None):
        """Draws the text by setting up an appropriate path on the Cairo
        context and filling it. `x` and `y` denote the coordinates where the
        drawing should start. If they are both C{None}, the current position
        of the context will be used.
        
        Vertical alignment settings are not taken into account in this method
        as the text is not drawn within a box."""
        ctx = self.context

        if x is None or y is None:
            x, y = ctx.get_current_point()

        line_height = ctx.font_extents()[2]
        coords = []

        if self.halign == self.CENTER:
            # Centered alignment
            width = self.text_extents()[2]
            for idx, line in enumerate(self._iterlines()):
                xb, _, line_width, _, _, _ = ctx.text_extents(line)
                coords.append((x + (width - line_width) / 2. - xb,
                               y + idx * line_height))
        elif self.halign == self.RIGHT:
            # Right alignment
            x += self.text_extents()[2]
            for idx, line in enumerate(self._iterlines()):
                xb, _, width, _, _, _ = ctx.text_extents(line)
                coords.append((x - width - xb, y + idx * line_height))
        else:
            # Left alignment
            coords = [
                (x - ctx.text_extents(line)[0], y + idx * line_height)
                for idx, line in enumerate(self._iterlines())
            ]

        for (x, y), line in zip(coords, self._iterlines()):
            ctx.move_to(x, y)
            ctx.show_text(line)

    def _iterlines(self):
        """Iterates over the label line by line"""
        return iter(self._text.split("\n"))

    @property
    def text(self):
        """Returns the text to be drawn."""
        return self._text

    @text.setter
    def text(self, text):
        """Sets the text that will be drawn.
        
        If `text` is ``None``, it will be mapped to an empty string; otherwise,
        it will be converted to a string."""
        if text is None:
            self._text = ""
        else:
            self._text = str(text)

    def text_extents(self):
        """Returns the X-bearing, Y-bearing, width, height, X-advance and
        Y-advance of the text.
        
        For multi-line text, the X-bearing and Y-bearing correspond to the
        first line, while the X-advance is extracted from the last line.
        and the Y-advance is the sum of all the Y-advances. The width and
        height correspond to the entire bounding box of the text."""
        lines = self.text.split("\n")
        if len(lines) <= 1:
            return self.context.text_extents(self.text)

        x_bearing, y_bearing, width, height, x_advance, y_advance = \
            self.context.text_extents(lines[0])

        for line in lines[1:]:
            _, _, w, h, x_advance, ya = self.context.text_extents(line)
            width = max(width, w)
            height += h
            y_advance += ya

        return x_bearing, y_bearing, width, height, x_advance, y_advance

def test():
    """Testing routine for L{TextDrawer}"""
    import cairo
    import math

    text = "The quick brown fox\njumped over the\nlazy dog"
    width, height = (300, 400)

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    context = cairo.Context(surface)
    drawer = TextDrawer(context, text)

    context.set_source_rgb(1, 1, 1)
    context.rectangle(0, 0, width, height)
    context.fill()

    context.set_source_rgb(0.7, 0.7, 0.7)
    for i in range(100, width, 100):
        context.move_to(i, 0)
        context.line_to(i, height)
        context.stroke()
    for i in range(100, height, 100):
        context.move_to(0, i)
        context.line_to(width, i)
        context.stroke()

    for i, halign in enumerate(("left", "center", "right")):
        # Mark the reference point
        context.set_source_rgba(0, 0, 1, 0.5)
        x, y = i * 100, 20
        context.arc(x, y, 2, 0, 2 * math.pi)
        context.fill()

        # Draw the text
        context.set_source_rgb(0, 0, 0)
        drawer.halign = halign
        drawer.draw(i * 100, 20)

        # Mark the new reference point
        context.set_source_rgba(1, 0, 0, 0.5)
        x, y = context.get_current_point()
        context.arc(x, y, 2, 0, 2 * math.pi)
        context.fill()

    surface.write_to_png("test.png")

if __name__ == "__main__":
    test()

