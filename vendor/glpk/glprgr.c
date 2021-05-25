/* glprgr.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#define _GLPSTD_ERRNO
#define _GLPSTD_STDIO
#include "glpenv.h"
#include "glprgr.h"
#define xfault xerror

/***********************************************************************
*  NAME
*
*  rgr_write_bmp16 - write 16-color raster image in BMP file format
*
*  SYNOPSIS
*
*  #include "glprgr.h"
*  int rgr_write_bmp16(const char *fname, int m, int n, const char
*     map[]);
*
*  DESCRIPTION
*
*  The routine rgr_write_bmp16 writes 16-color raster image in
*  uncompressed BMP file format (Windows bitmap) to a binary file whose
*  name is specified by the character string fname.
*
*  The parameters m and n specify, respectively, the number of rows and
*  the numbers of columns (i.e. height and width) of the raster image.
*
*  The character array map has m*n elements. Elements map[0, ..., n-1]
*  correspond to the first (top) scanline, elements map[n, ..., 2*n-1]
*  correspond to the second scanline, etc.
*
*  Each element of the array map specifies a color of the corresponding
*  pixel as 8-bit binary number XXXXIRGB, where four high-order bits (X)
*  are ignored, I is high intensity bit, R is red color bit, G is green
*  color bit, and B is blue color bit. Thus, all 16 possible colors are
*  coded as following hexadecimal numbers:
*
*     0x00 = black         0x08 = dark gray
*     0x01 = blue          0x09 = bright blue
*     0x02 = green         0x0A = bright green
*     0x03 = cyan          0x0B = bright cyan
*     0x04 = red           0x0C = bright red
*     0x05 = magenta       0x0D = bright magenta
*     0x06 = brown         0x0E = yellow
*     0x07 = light gray    0x0F = white
*
*  RETURNS
*
*  If no error occured, the routine returns zero; otherwise, it prints
*  an appropriate error message and returns non-zero. */

static void put_byte(FILE *fp, int c)
{     fputc(c, fp);
      return;
}

static void put_word(FILE *fp, int w)
{     /* big endian */
      put_byte(fp, w);
      put_byte(fp, w >> 8);
      return;
}

static void put_dword(FILE *fp, int d)
{     /* big endian */
      put_word(fp, d);
      put_word(fp, d >> 16);
      return;
}

int rgr_write_bmp16(const char *fname, int m, int n, const char map[])
{     FILE *fp;
      int offset, bmsize, i, j, b, ret = 0;
      if (!(1 <= m && m <= 32767))
         xfault("rgr_write_bmp16: m = %d; invalid height\n", m);
      if (!(1 <= n && n <= 32767))
         xfault("rgr_write_bmp16: n = %d; invalid width\n", n);
      fp = fopen(fname, "wb");
      if (fp == NULL)
      {  xprintf("rgr_write_bmp16: unable to create `%s' - %s\n",
            fname, strerror(errno));
         ret = 1;
         goto fini;
      }
      offset = 14 + 40 + 16 * 4;
      bmsize = (4 * n + 31) / 32;
      /* struct BMPFILEHEADER (14 bytes) */
      /* UINT bfType */          put_byte(fp, 'B'), put_byte(fp, 'M');
      /* DWORD bfSize */         put_dword(fp, offset + bmsize * 4);
      /* UINT bfReserved1 */     put_word(fp, 0);
      /* UNIT bfReserved2 */     put_word(fp, 0);
      /* DWORD bfOffBits */      put_dword(fp, offset);
      /* struct BMPINFOHEADER (40 bytes) */
      /* DWORD biSize */         put_dword(fp, 40);
      /* LONG biWidth */         put_dword(fp, n);
      /* LONG biHeight */        put_dword(fp, m);
      /* WORD biPlanes */        put_word(fp, 1);
      /* WORD biBitCount */      put_word(fp, 4);
      /* DWORD biCompression */  put_dword(fp, 0 /* BI_RGB */);
      /* DWORD biSizeImage */    put_dword(fp, 0);
      /* LONG biXPelsPerMeter */ put_dword(fp, 2953 /* 75 dpi */);
      /* LONG biYPelsPerMeter */ put_dword(fp, 2953 /* 75 dpi */);
      /* DWORD biClrUsed */      put_dword(fp, 0);
      /* DWORD biClrImportant */ put_dword(fp, 0);
      /* struct RGBQUAD (16 * 4 = 64 bytes) */
      /* CGA-compatible colors: */
      /* 0x00 = black */         put_dword(fp, 0x000000);
      /* 0x01 = blue */          put_dword(fp, 0x000080);
      /* 0x02 = green */         put_dword(fp, 0x008000);
      /* 0x03 = cyan */          put_dword(fp, 0x008080);
      /* 0x04 = red */           put_dword(fp, 0x800000);
      /* 0x05 = magenta */       put_dword(fp, 0x800080);
      /* 0x06 = brown */         put_dword(fp, 0x808000);
      /* 0x07 = light gray */    put_dword(fp, 0xC0C0C0);
      /* 0x08 = dark gray */     put_dword(fp, 0x808080);
      /* 0x09 = bright blue */   put_dword(fp, 0x0000FF);
      /* 0x0A = bright green */  put_dword(fp, 0x00FF00);
      /* 0x0B = bright cyan */   put_dword(fp, 0x00FFFF);
      /* 0x0C = bright red */    put_dword(fp, 0xFF0000);
      /* 0x0D = bright magenta */ put_dword(fp, 0xFF00FF);
      /* 0x0E = yellow */        put_dword(fp, 0xFFFF00);
      /* 0x0F = white */         put_dword(fp, 0xFFFFFF);
      /* pixel data bits */
      b = 0;
      for (i = m - 1; i >= 0; i--)
      {  for (j = 0; j < ((n + 7) / 8) * 8; j++)
         {  b <<= 4;
            b |= (j < n ? map[i * n + j] & 15 : 0);
            if (j & 1) put_byte(fp, b);
         }
      }
      fflush(fp);
      if (ferror(fp))
      {  xprintf("rgr_write_bmp16: write error on `%s' - %s\n",
            fname, strerror(errno));
         ret = 1;
      }
fini: if (fp != NULL) fclose(fp);
      return ret;
}

/* eof */
