"""
Drawers for labels on plots.
"""

from igraph.drawing.baseclasses import AbstractCairoDrawer

__all__ = ["TextDrawer"]
__license__ = "GPL"

__docformat__ = "restructuredtext en" 

#####################################################################

class TextDrawer(AbstractCairoDrawer):
    """Class that draws text on a Cairo context.
    
    This class supports multi-line text unlike the original Cairo text
    drawing methods."""

    LEFT, CENTER, RIGHT = "left", "center", "right"
    TOP, BOTTOM = "top", "bottom"

    def __init__(self, context, text="", halign="center", valign="center"):
        """Constructs a new instance that will draw the given `text` on
        the given Cairo `context`."""
        super(TextDrawer, self).__init__(context, (0, 0))
        self.text = text
        self.halign = halign
        self.valign = valign

    def draw(self):
        """Draws the text in the current bounding box of the drawer.
        
        Since the class itself is an instance of `AbstractCairoDrawer`, it
        has an attribute named ``bbox`` which will be used as a bounding
        box."""
        line_height = self.context.font_extents()[2]

        bbox = self.bbox
        left, top, width = bbox.left, bbox.top, bbox.width

        if self.valign == self.BOTTOM:
            # Bottom vertical alignment
            _, yb, _, total_height = self.text_extents()[:4]
            top += bbox.height - total_height - yb
        elif self.valign == self.CENTER:
            # Centered vertical alignment
            _, yb, _, total_height = self.text_extents()[:4]
            # total_height = self.text_extents()[3]
            top += (bbox.height - total_height - yb / 2. + line_height) / 2.
        else:
            # Top vertical alignment
            top += line_height

        return self.draw_at(left, top, width)

    def draw_at(self, x = None, y = None, width = None):
        """Draws the text by setting up an appropriate path on the Cairo
        context and filling it. `x` and `y` denote the coordinates where the
        drawing should start. If they are both C{None}, the current position
        of the context will be used.
        
        Vertical alignment settings are not taken into account in this method
        as the text is not drawn within a box.
        
        :Parameters:
          x: float or ``None``
            The X coordinate of the reference point where the drawing should
            start.
          y: float or ``None``
            The Y coordinate of the reference point where the drawing should
            start.
          width: float or ``None``
            The width of the box in which the text will be fitted. It matters
            only when the text is right-aligned or centered. The text will
            overflow the box if any of the lines is longer than the box width.
        """
        ctx = self.context

        if x is None or y is None:
            x, y = ctx.get_current_point()

        line_height = ctx.font_extents()[2]
        coords = []

        if self.halign == self.CENTER:
            # Centered alignment
            if width is None:
                width = self.text_extents()[2]
            for idx, line in enumerate(self._iterlines()):
                xb, _, line_width, _, _, _ = ctx.text_extents(line)
                coords.append((x + (width - line_width) / 2. - xb,
                               y + idx * line_height))
        elif self.halign == self.RIGHT:
            # Right alignment
            if width is None:
                width = self.text_extents()[2]
            x += width
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

        line_height = self.context.font_extents()[2]
        for line in lines[1:]:
            _, _, w, _, x_advance, ya = self.context.text_extents(line)
            width = max(width, w)
            height += line_height
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

    context.set_source_rgb(0.5, 0.5, 0.5)
    for i in range(100, width, 100):
        context.move_to(i, 0)
        context.line_to(i, height)
        context.stroke()
    for i in range(100, height, 100):
        context.move_to(0, i)
        context.line_to(width, i)
        context.stroke()
    context.set_source_rgb(0.75, 0.75, 0.75)
    context.set_line_width(1)
    for i in range(50, width, 100):
        context.move_to(i, 0)
        context.line_to(i, height)
        context.stroke()
    for i in range(50, height, 100):
        context.move_to(0, i)
        context.line_to(width, i)
        context.stroke()

    def mark_point(red, green, blue):
        """Marks the current point on the canvas by the given color"""
        x, y = context.get_current_point()
        context.set_source_rgba(red, green, blue, 0.5)
        context.arc(x, y, 2, 0, 2 * math.pi)
        context.fill()

    for i, halign in enumerate(("left", "center", "right")):
        context.move_to(i * 100, 20)

        # Mark the reference point
        mark_point(0, 0, 1)

        # Draw the text
        context.set_source_rgb(0, 0, 0)
        drawer.halign = halign
        drawer.draw_at(i * 100, 20)

        # Mark the new reference point
        mark_point(1, 0, 0)

    for i, halign in enumerate(("left", "center", "right")):
        for j, valign in enumerate(("top", "center", "bottom")):
            # Draw the text
            context.set_source_rgb(0, 0, 0)
            drawer.halign = halign
            drawer.valign = valign
            drawer.bbox = (i*100, j*100+100, i*100+100, j*100+200)
            drawer.draw()
            # Mark the new reference point
            mark_point(1, 0, 0)

    surface.write_to_png("test.png")

if __name__ == "__main__":
    test()

