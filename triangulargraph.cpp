#include <graphics.h>
#include <iostream>
  
// Driver code
int main()
{
    // gm is Graphics mode which
    // is a computer display
    // mode that generates
    // image using pixels.
    // DETECT is a macro
    // defined in "graphics.h"
    // header file
    int gd = DETECT, gm;
  
    // initgraph initializes
    // the graphics system
    // by loading a graphics
    // driver from disk
    initgraph(&gd, &gm, "");
  
    // Triangle
  
    // line for x1, y1, x2, y2
    line(150, 150, 450, 150);
  
    // line for x1, y1, x2, y2
    line(150, 150, 300, 300);
  
    // line for x1, y1, x2, y2
    line(450, 150, 300, 300);
  
    // closegraph function closes
    // the graphics mode and
    // deallocates all memory
    // allocated by graphics system
    getch();
  
    // Close the initialized gdriver
    closegraph();
}