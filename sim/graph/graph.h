/****************************************************************

 graph.h

 =============================================================

 Copyright 1996-2018 Tom Barbalet. All rights reserved.

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or
 sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.

 This software and Noble Ape are a continuing work of Tom Barbalet,
 begun on 13 June 1996. No apes or cats were harmed in the writing
 of this software.

 ****************************************************************/

/*NOBLEMAKE VAR=""*/

#ifndef NOBLEAPE_GRAPH_H
#define NOBLEAPE_GRAPH_H

#ifndef	_WIN32

#include "../noble/noble.h"
#include "../universe/universe.h"
#include "../entity/entity.h"

#else

#include "..\noble\noble.h"
#include "..\universe.h\universe.h"
#include "..\entity\entity.h"

#endif

typedef struct
{
    n_byte r;
    n_byte g;
    n_byte b;
}n_rgb;

typedef struct
{
    n_vect2 * points;
    n_int no_of_points;
    n_int max_points;
}n_points;

void graph_vascular(n_genetics * being,
                    n_byte * buffer,
                    n_vect2* img,
                    n_int tx, n_int ty, n_int bx, n_int by,
                    n_byte thickness,
                    n_byte clear,
                    n_int shoulder_angle, n_int elbow_angle, n_int wrist_angle,
                    n_int hip_angle, n_int knee_angle,
                    n_byte show_skeleton_keypoints);

void graph_command(n_int gc_val);


void graph_erase(n_byte * buffer, n_vect2 * img);

/* draws a line */
void graph_line(n_byte * buffer,
                n_vect2 * img,
                n_vect2 * previous,
                n_vect2 * current,
                n_rgb * color,
                n_byte thickness);

void graph_curve(n_byte * buffer,
                 n_vect2 * img,
                 n_int x0, n_int y0,
                 n_int x1, n_int y1,
                 n_int x2, n_int y2,
                 n_rgb * color,
                 n_byte radius_percent,
                 n_uint start_thickness,
                 n_uint end_thickness);

void graph_fill_polygon(n_vect2 * points, n_int no_of_points,
                        n_byte r, n_byte g, n_byte b, n_byte transparency,
                        n_byte * buffer, n_vect2 * img);

void outline_points(const n_vect2 * source_points,
                    n_int no_of_source_points, n_int extra_points,
                    n_int x, n_int y,
                    n_byte mirror,
                    n_vect2 * scale,
                    n_int angle,
                    n_vect2 *axis,
                    n_vect2 *extra_1,
                    n_vect2 *extra_2,
                    n_vect2 *extra_3,
                    n_vect2 *extra_4,
                    n_points * collection);


#endif /* NOBLEAPE_GRAPH_H */

/*NOBLEMAKE END=""*/


