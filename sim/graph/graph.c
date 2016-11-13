/****************************************************************

 graph.c

 =============================================================

 Copyright 1996-2013 Tom Barbalet. All rights reserved.

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

#ifndef	_WIN32

#include "../noble/noble.h"
#include "../universe/universe.h"
#include "../entity/entity.h"

#else

#include "..\noble\noble.h"
#include "..\universe.h\universe.h"
#include "..\entity\entity.h"


#endif

#include <math.h>

#include "phosphene.h"
#include "graph.h"


static n_int  graph_state = GC_VASCULAR;
static n_byte graph_clear = 1;


void graph_erase(n_byte * buffer, n_vect2 * img)
{
    n_int i;
    for (i = 0; i < img->x *img->y *3; i++)
    {
        buffer[i] = 255;
    }
}

/* draws a line */
void graph_line(n_byte * buffer,
                       n_vect2 * img,
                       n_vect2 * previous,
                       n_vect2 * current,
                       n_rgb * color,
                       n_byte thickness)
{
    n_int i,max;
    n_vect2 delta;
    n_vect2 absdelta;
    
    
    vect2_subtract(&delta, current, previous);

    vect2_copy(&absdelta, &delta);
    
    
    
    if (absdelta.x < 0) absdelta.x = -delta.x;
    if (absdelta.y < 0) absdelta.y = -delta.y;
    
    max = absdelta.x;
    if (absdelta.y > max) max = absdelta.y;
    
    for (i=0; i<max; i++)
    {
        n_int xx = previous->x + (i*(current->x - previous->x)/max);
        if ((xx > -1) && (xx < img->x))
        {
            n_int yy = previous->y + (i*(current->y-previous->y)/max);
            if ((yy > -1) && (yy < img->y))
            {
                n_int n = (yy*img->x + xx)*3;
                buffer[n] = color->r;
                buffer[n+1] = color->g;
                buffer[n+2] = color->b;
            }
        }
    }
}


/**
 * @brief Draws a curve using three points
 * @param buffer Image buffer (three bytes per pixel)
 * @param img_width Width of the image
 * @param img_height Height of the image
 * @param x0 x coordinate of the start point
 * @param y0 y coordinate of the start point
 * @param x1 x coordinate of the middle point
 * @param y1 y coordinate of the middle point
 * @param x0 x coordinate of the end point
 * @param y0 y coordinate of the end point
 * @param r red
 * @param g green
 * @param b blue
 * @param radius_percent Radius of the curve as a percentage
 * @param start_thickness Thickness of the curve at the start point
 * @param end_thickness Thickness of the curve at the end point
 */
void graph_curve(n_byte * buffer,
                        n_vect2 * img,
                        n_int x0, n_int y0,
                        n_int x1, n_int y1,
                        n_int x2, n_int y2,
                        n_rgb * color,
                        n_byte radius_percent,
                        n_uint start_thickness,
                        n_uint end_thickness)
{
    n_int pts[8];
    
    n_vect2 current;
    n_vect2 previous = {0, 0};
    
    n_uint i;
    const n_uint divisions = 20;
    double c[5],d[5],f;
    
    /** turn three points into four using the curve radius */
    pts[0] = x0;
    pts[1] = y1;
    
    pts[2] = x1 + ((x0 - x1)*radius_percent/100);
    pts[3] = y1 + ((y0 - y1)*radius_percent/100);
    
    pts[4] = x1 + ((x2 - x1)*radius_percent/100);
    pts[5] = y1 + ((y2 - y1)*radius_percent/100);
    
    pts[6] = x2;
    pts[7] = y2;
    
    c[0] = (-pts[0*2] + 3 * pts[1*2] - 3 * pts[2*2] + pts[3*2]) / 6.0;
    c[1] = (3 * pts[0*2] - 6 * pts[1*2] + 3 * pts[2*2]) / 6.0;
    c[2] = (-3 * pts[0*2] + 3 * pts[2*2]) / 6.0;
    c[3] = (pts[0*2] + 4 * pts[1*2] + pts[2*2]) / 6.0;
    
    d[0] = (-pts[(0*2)+1] + 3 * pts[(1*2)+1] - 3 * pts[(2*2)+1] + pts[(3*2)+1]) / 6.0;
    d[1] = (3 * pts[(0*2)+1] - 6 * pts[(1*2)+1] + 3 * pts[(2*2)+1]) / 6.0;
    d[2] = (-3 * pts[(0*2)+1] + 3 * pts[(2*2)+1]) / 6.0;
    d[3] = (pts[(0*2)+1] + 4 * pts[(1*2)+1] + pts[(2*2)+1]) / 6.0;
    
    for (i = 0; i < divisions; i++)
    {
        f = (double)i / (double)divisions;
        current.x = (n_int)((c[2] + f * (c[1] + f * c[0])) * f + c[3]);
        current.y = (n_int)((d[2] + f * (d[1] + f * d[0])) * f + d[3]);
        
        if (i > 0)
        {
            graph_line(buffer, img,
                       &previous, &current,
                       color,
                       start_thickness +
                       ((end_thickness - start_thickness) * i / divisions));
        }
        vect2_copy(&previous, &current);
    }
    
}

#define  MAX_POLYGON_CORNERS 1000

/**
 * @brief Draw a filled polygon
 * @param points Array containing 2D points
 * @param no_of_points The number of 2D points
 * @param r Red
 * @param g Green
 * @param b Blue
 * @param transparency Degree of transparency
 * @param buffer Image buffer (3 bytes per pixel)
 * @param img_width Image width
 * @param img_height Image height
 */
void graph_fill_polygon(n_vect2 * points, n_int no_of_points,
                               n_byte r, n_byte g, n_byte b, n_byte transparency,
                               n_byte * buffer, n_vect2 * img)
{
    n_int nodes, nodeX[MAX_POLYGON_CORNERS], i, j, swap, n, x, y;
    n_int min_x = 99999, min_y = 99999;
    n_int max_x = -99999, max_y = -99999;
    
    for (i = 0; i < no_of_points; i++)
    {
        x = points[i].x;
        y = points[i].y;
        if ((x==9999) || (y==9999)) continue;
        if (x < min_x) min_x = x;
        if (y < min_y) min_y = y;
        if (x > max_x) max_x = x;
        if (y > max_y) max_y = y;
    }
    
    if (min_x < 0) min_x = 0;
    if (min_y < 0) min_y = 0;
    if (max_x >= img->x) max_x = img->x - 1;
    if (max_y >= img->y) max_y = img->y - 1;
    
    for (y = min_y; y <= max_y; y++)
    {
        /**  Build a list of nodes */
        nodes = 0;
        j = no_of_points-1;
        for (i = 0; i < no_of_points; i++)
        {
            if (((points[i].y < y) && (points[j].y >= y)) ||
                ((points[j].y < y) && (points[i].y >= y)))
            {
                nodeX[nodes++] =
                points[i].x + (y - points[i].y) *
                (points[j].x - points[i].x) /
                (points[j].y - points[i].y);
            }
            j = i;
            if (nodes == MAX_POLYGON_CORNERS) break;
        }
        
        /**  Sort the nodes, via a simple “Bubble” sort */
        i = 0;
        while (i < nodes-1)
        {
            if (nodeX[i] > nodeX[i+1])
            {
                swap = nodeX[i];
                nodeX[i] = nodeX[i+1];
                nodeX[i+1] = swap;
                if (i) i--;
            }
            else
            {
                i++;
            }
        }
        
        /**  Fill the pixels between node pairs */
        for (i = 0; i < nodes; i += 2)
        {
            if (nodeX[i] >= max_x) break;
            if (nodeX[i+1] > min_x)
            {
                /** range check */
                if (nodeX[i] <= min_x) nodeX[i] = min_x+1;
                if (nodeX[i+1] >= max_x) nodeX[i+1] = max_x-1;
                
                for (x = nodeX[i]; x < nodeX[i+1]; x++)
                {
                    n = ((y*img->x)+x)*3;
                    if (transparency == 0)
                    {
                        buffer[n] = b;
                        buffer[n+1] = g;
                        buffer[n+2] = r;
                    }
                    else
                    {
                        buffer[n]   = ((b*(255-transparency)) + (buffer[n]*transparency))/256;
                        buffer[n+1] = ((g*(255-transparency)) + (buffer[n+1]*transparency))/256;
                        buffer[n+2] = ((r*(255-transparency)) + (buffer[n+2]*transparency))/256;
                    }
                }
            }
        }
    }
}

/**
 * @brief Returns an array of 2D points used for drawing diagrams
 * @param source_points Array of 2D points which is the template
 * @param no_of_source_points Number of 2D points in the template
 * @param extra_points The number of points to be returned via the extra parameters
 * @param x The starting x coordinate
 * @param y The starting y coordinate
 * @param mirror Flip in the vertical axis
 * @param scale_width Horizontal scaling factor x1000
 * @param scale_length Length (vertical) scaling factor x1000
 * @param angle Rotation angle of the result
 * @param axis_x Returned x coordinate such that (x,y)-(axis_x,axis_y) defines the axis of the object
 * @param axis_y Returned y coordinate such that (x,y)-(axis_x,axis_y) defines the axis of the object
 * @param extra_x1
 * @param extra_y1
 * @param extra_x2
 * @param extra_y2
 * @param extra_x3
 * @param extra_y3
 * @param extra_x4
 * @param extra_y4
 * @param points Returned 2D points
 * @param no_of_points Number of returned 2D points
 * @param max_points The maximum number of points which may be returned
 */

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
                           n_points * collection)
{
    n_vect2  ds, location, vector;
    n_int    i, axis_length,point_length;
    float axis_angle, point_angle;
    float ang = angle*TWO_PI/7200;
    
    vect2_populate(&location, x, y);
    vect2_subtract(&ds, (n_vect2 *)&source_points[1], (n_vect2 *)source_points);
    vect2_multiplier(&ds, &ds, scale, 1, 1000);
    
    /** length of the object */
    axis_length = (n_int)math_root(vect2_dot(&ds, &ds, 1, 1));
    if (axis_length < 1) axis_length=1;
    
    /** invert around the vertical axis if needed */
    if (mirror != 0)
    {
        ds.x = -ds.x;
    }
    
    /** find the orientation angle of the axis */
    axis_angle = (float)acos(ds.x/(float)axis_length);
    
    if (ds.y < 0)
    {
        axis_angle = TWO_PI-axis_angle;
    }
    
    vect2_populate(&vector, (n_int)(axis_length*sin(ang+(TWO_PI/4)-axis_angle)),
                   (n_int)(axis_length*cos(ang+(TWO_PI/4)-axis_angle)));
    
    /** calculate the position of the end point of the axis */
    
    vect2_add(axis, &location, &vector);
    
    /** draw lines between each point */
    for (i = 2; i < no_of_source_points + 2 + extra_points; i++)
    {
        n_vect2 point;
        vect2_subtract(&ds, (n_vect2 *)&source_points[i], (n_vect2 *)source_points);
        vect2_multiplier(&ds, &ds, scale, 1, 1000);
        point_length = (n_int)math_root(vect2_dot(&ds, &ds, 1, 1));
        if (point_length < 1)
        {
            point_length = 1;
        }
        
        /** invert the line around the vertical axis if necessary */
        if (mirror != 0)
        {
            ds.x = -ds.x;
        }
        
        /** angle of the line */
        point_angle = (float)acos(ds.x/(float)point_length);
        if (ds.y < 0)
        {
            point_angle = (TWO_PI)-point_angle;
        }
        
        /** position of the end of the line */
        vect2_populate(&vector, (n_int)(point_length*sin(ang+point_angle-axis_angle)),
                                (n_int)(point_length*cos(ang+point_angle-axis_angle)));
        
        vect2_add(&point, &location, &vector);
        
        /** store the calculated point positions in an array */
        if (collection->no_of_points < collection->max_points)
        {
            if (i < no_of_source_points + 2)
            {
                vect2_copy(&collection->points[collection->no_of_points], &point);
                collection->no_of_points++;
            }
        }
        /* TODO
        else
        {
            (void)SHOW_ERROR("Maximum number of skeleton points reached");
        }
        */
        /** This is a crude way of keeping track of the last few points
         so that they can be returned by the function */
        vect2_copy(extra_1, extra_2);
        vect2_copy(extra_2, extra_3);
        vect2_copy(extra_3, extra_4);
        vect2_copy(extra_4, &point);
    }
    
    if (collection->no_of_points > -1 && collection->no_of_points < collection->max_points)
    {
        collection->points[collection->no_of_points].x = 9999;
        collection->points[collection->no_of_points].y = 9999;
        collection->no_of_points++;
    }
    else
    {
        SHOW_ERROR("Outside point range for drawing");
    }
}



/**
 * @brief Shows distribution of honor.  Note that beings are sorted in order of honor
 * @param sim Pointer to the simulation object
 * @param update type Whether to draw the background or foreground
 * @param buffer Image buffer (3 bytes per pixel)
 * @param img_width Image width
 * @param img_height Image height
 */
void graph_honor_distribution(noble_simulation * sim, n_byte update_type, n_byte * buffer, n_int img_width, n_int img_height)
{
    n_uint i,j,temp;
    n_int * honor_value;
    scope s;
    unsigned int channel = 0;
    unsigned int intensity_percent = 100;
    unsigned int grid_horizontal = 10;
    unsigned int grid_vertical = 8;
    
    s = create_scope((unsigned int)1);
    s.time_ms = (unsigned int)(sim->num);
    s.noise = 0.5;
    
    if (update_type == PHOSPHENE_DRAW_BACKGROUND) {
        scope_draw(&s, update_type, intensity_percent,
                   grid_horizontal, grid_vertical,
                   (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
        return;
    }
    
    honor_value = (n_int*)io_new(sim->num*sizeof(n_int));
    
    /** get the honor values */
    for (i = 0; i < sim->num; i++)
    {
        noble_being * local_being = &(sim->beings[i]);
        honor_value[i] = being_honor(local_being);
    }
    
    /** sort the honor values */
    for (i = 0; i < sim->num; i++)
    {
        for (j = i+1; j < sim->num; j++)
        {
            if (honor_value[i]<honor_value[j]) {
                temp = honor_value[i];
                honor_value[i] = honor_value[j];
                honor_value[j] = temp;
            }
        }
    }
    
    
    for (i = 0; i < sim->num; i++)
    {
        scope_update(&s, channel, (double)honor_value[i], 0.0, 256.0, (unsigned int)i);
    }
    
    io_free(honor_value);
    
    scope_draw(&s, update_type, intensity_percent,
               grid_horizontal, grid_vertical,
               (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
}

static n_int graph_being_score(noble_simulation * sim, noble_being * local_being, n_byte score_type)
{
    n_int nucleotide,i,score = 0;
    n_byte * bases = (n_byte*)being_genetics(local_being);
    
    switch (score_type)
    {
        case 0:
            for (i = 0; i < BRAINCODE_SIZE; i++)
            {
                score += being_braincode_external(local_being)[i] + being_braincode_internal(local_being)[i];
            }
            break;
        case 1:
            for (i = 0; i < CHROMOSOMES; i++)
            {
                for (nucleotide = 0; nucleotide < 8; nucleotide++)
                {
                    score += (bases[i]>>(nucleotide*2)) & 3;
                }
            }
            break;
    }
    return score;
}

/*
 Updates an index array which is used to sort beings in order of honor value
 */
static void graph_being_index(noble_simulation * sim, n_int *index, n_byte score_type)
{
#ifdef BRAINCODE_ON
#ifdef PARASITES_ON
    n_int score;
    n_uint i;
    n_byte * used = (unsigned char*)io_new(sim->num);
    for (i = 0; i < sim->num; i++)
    {
        used[i] = 0;
        index[i] = i;
    }
    /* sort by honor value */
    for (i = 0; i < sim->num; i++)
    {
        if (used[i]==0)
        {
            n_uint j;
            n_int max = -1;
            n_int idx = -1;
            for (j = 0; j < sim->num; j++)
            {
                if (used[j]==0)
                {
                    noble_being * local_being = &(sim->beings[j]);
                    score = graph_being_score(sim, local_being,score_type);
                    if (score>max)
                    {
                        max = score;
                        idx = j;
                    }
                }
            }
            if (idx>-1)
            {
                index[i] = idx;
                used[idx] = 1;
            }
        }
    }
    io_free((void*)used);
#endif
#endif
}

/*
 Displays the braincode programs (one per row) for each being in the population.
 Colors represent different instruction types, and individuals are sorted by honor.
 */
void graph_ideosphere(noble_simulation * sim, n_byte * buffer, n_int img_width, n_int img_height)
{
#ifdef BRAINCODE_ON
#ifdef PARASITES_ON
    n_int i,x,y,n,half_width,max_height=img_height;
    n_int *index = (n_int*)io_new((sim->num)*sizeof(n_int));
    noble_being * local_being;
    n_byte * code;
    
    graph_being_index(sim, index,0);
    
    half_width = img_width/2;
    n = 0;
    if (sim->num>0)
    {
        for (y = 0; y < max_height; y++)
        {
            i = index[y*(sim->num-1)/max_height];
            local_being = &(sim->beings[i]);
            for (x = 0; x < img_width; x++, n+=3)
            {
                if (x<half_width)
                {
                    i = (x * ((BRAINCODE_SIZE/BRAINCODE_BYTES_PER_INSTRUCTION)-1) / half_width)*BRAINCODE_BYTES_PER_INSTRUCTION;
                    code = being_braincode_internal(local_being);
                }
                else
                {
                    i = ((x-half_width) * ((BRAINCODE_SIZE/BRAINCODE_BYTES_PER_INSTRUCTION)-1) / half_width)*BRAINCODE_BYTES_PER_INSTRUCTION;
                    code = being_braincode_external(local_being);
                }
                buffer[n] = code[i];
                buffer[n+1] = code[i+1];
                buffer[n+2] = code[i+2];
            }
        }
    }
    
    io_free((void*)index);
#endif
#endif
}

/*
 Shows the genome for each individual in the population, with one genome per row.
 Individuals are sorted in order of honor.
 */
void graph_genepool(noble_simulation * sim, n_byte * buffer, n_int img_width, n_int img_height)
{
#ifdef PARASITES_ON
    n_int i,x,y,n,ch,idx,nucleotide,max_height=img_height;
    n_genetics * bases;
    n_int * index;
    noble_being * local_being;
    
    const n_byte col[] =
    {
        200,0,0,	0,200,0,	0,0,200,	200,200,0
    };
    
    index = (n_int*)io_new(sim->num*sizeof(n_int));
    graph_being_index(sim,index,1);
    
    n = 0;
    if (sim->num>0)
    {
        for (y = 0; y < max_height; y++)
        {
            i = index[y*(sim->num-1)/max_height];
            local_being = &(sim->beings[i]);
            bases = being_genetics(local_being);
            for (x = 0; x < img_width; x++, n+=3)
            {
                nucleotide = x * (CHROMOSOMES*16) / img_width;
                ch = nucleotide>>4;
                idx = ((bases[ch]>>((nucleotide-(16*ch))*2))&3)*3;
                buffer[n] = col[idx++];
                buffer[n+1] = col[idx++];
                buffer[n+2] = col[idx];
            }
        }
    }
    
    io_free((void*)index);
#endif
}

/*
 A matrix showing the relationships between beings.
 Green squares represent friendly relationships, red represent unfriendly relationships
 and black represents "don't care"
 Individuals are plotted on each axis in the same (honor sorted) order, such that a being's
 relationship with itself is along the diagonal axis
 */
void graph_relationship_matrix(noble_simulation * sim, n_byte * buffer, n_int img_width, n_int img_height)
{
#ifdef PARASITES_ON
    n_int j,x,y;
    n_uint i, k;
    n_int *index = (n_int*)io_new(sim->num*sizeof(n_int));
    n_vect2 img;

    graph_being_index(sim,index,1);
    
    img.x = img_width;
    img.y = img_height;
    graph_erase(buffer, &img);
    
    for (i = 0; i < sim->num; i++)
    {
        noble_being * local_being = &(sim->beings[index[i]]);
        social_link * graph = being_social(local_being);
        n_uint respect_threshold = social_respect_mean(sim, local_being);
        n_int tx = i*(img.x - 1)/sim->num;
        n_int bx = (i+1)*(img.x - 1)/sim->num;
        for (j = 0; j < SOCIAL_SIZE_BEINGS; j++)
        {
            if (!SOCIAL_GRAPH_ENTRY_EMPTY(graph, j))
            {
                if (j == 0)
                {
                    k = i;
                }
                else
                {
                    for (k = 0; k < sim->num; k++)
                    {
                        noble_being * local_being2 = &(sim->beings[index[k]]);
                        if ((being_family_name(local_being2) == graph[j].family_name[BEING_MET]) &&
                            (being_first_name(local_being2) == UNPACK_FAMILY_FIRST_NAME(graph[j].first_name[BEING_MET])))
                        {
                            break;
                        }
                    }
                }
                
                if (k < sim->num)
                {
                    n_int ty = k*(img.y-1)/sim->num;
                    n_int by = (k+1)*(img.y-1)/sim->num;
                    
                    for (y=ty; y<by; y++)
                    {
                        for (x=tx; x<bx; x++)
                        {
                            n_int v, n = ((y*img.x)+x)*3;
                            if (graph[j].friend_foe>=respect_threshold)
                            {
                                v = (graph[j].friend_foe - respect_threshold)*8;
                                if (v>255) v=255;
                                buffer[n] = 0;
                                buffer[n+1] = (n_byte)v;
                                buffer[n+2] = 0;
                            }
                            else
                            {
                                v = (respect_threshold - graph[j].friend_foe)*8;
                                if (v>255) v=255;
                                buffer[n] = (n_byte)v;
                                buffer[n+1] = 0;
                                buffer[n+2] = 0;
                            }
                        }
                    }
                }
            }
        }
    }
    
    io_free((void*)index);
#endif
}

/*
 Draws the number of antigens and antibodies in the population
 There are 256 possible antigens/antibodies which are along the horizontal axis.
 Antigens are shown in red and antibodies in green.
 */
void graph_pathogens(noble_simulation * sim, n_byte update_type, n_byte * buffer, n_int img_width, n_int img_height)
{
    n_c_uint i;
    n_c_uint * antibodies;
    n_c_uint * antigens;
#ifdef IMMUNE_ON
    noble_being * local_being;
    n_int j,p;
    n_uint max_val=1;
    noble_immune_system * immune;
    scope s;
    unsigned int intensity_percent = 100;
    unsigned int grid_horizontal = 10;
    unsigned int grid_vertical = 8;
#endif
    
    s = create_scope((unsigned int)1);
    s.time_ms = (unsigned int)256;
    s.no_of_traces = 2;
    s.noise = 200;
    
    if (update_type == PHOSPHENE_DRAW_BACKGROUND) {
        scope_draw(&s, update_type, intensity_percent,
                   grid_horizontal, grid_vertical,
                   (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
        return;
    }
    
    antibodies = (n_c_uint*)io_new(256*sizeof(n_c_uint));
    antigens = (n_c_uint*)io_new(256*sizeof(n_c_uint));
    
    for (i=0; i<256; i++)
    {
        antibodies[i]=0;
        antigens[i]=0;
    }
    
#ifdef IMMUNE_ON
    
    if (sim->num>0)
    {
        /* update histograms */
        for (p=0; p<256; p++)
        {
            for (i=0; i<sim->num; i++)
            {
                local_being = &(sim->beings[i]);
                immune = &(local_being->immune_system);
                for (j=0; j<IMMUNE_POPULATION; j++)
                {
                    antibodies[immune->shape_antibody[j]]++;
                }
                for (j=0; j<IMMUNE_ANTIGENS; j++)
                {
                    antigens[immune->shape_antigen[j]]++;
                }
            }
        }
        
        /* find the maximum value */
        for (p=0; p<256; p++)
        {
            if (antibodies[p]>max_val)
            {
                max_val=antibodies[p];
            }
            if (antigens[p]>max_val)
            {
                max_val=antigens[p];
            }
        }
        
        for (p=0; p<256; p++)
        {
            scope_update(&s, 0, (double)antibodies[p], 0.0, (double)max_val, (unsigned int)p);
            scope_update(&s, 1, (double)antigens[p], 0.0, (double)max_val, (unsigned int)p);
        }
    }
    
    scope_draw(&s, update_type, intensity_percent,
               grid_horizontal, grid_vertical,
               (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
#endif
    
    io_free((void*)antibodies);
    io_free((void*)antigens);
}

/*
 Shows a histogram of ages
 */
void graph_age_demographic(noble_simulation * sim, n_byte * buffer, n_vect2* img)
{
    const n_int max_age = AGE_OF_MATURITY*4;
    const n_int divisor = 4;
    n_int groups = max_age/divisor;
    n_int * age_group;
    n_int i,idx,max=1;
    n_int prev_x = 0, prev_y=img->y-1;
    n_int x,y;
    n_uint current_date;
    n_uint j;
    
    age_group = (n_int*)io_new(groups*sizeof(n_int));
    for (i=0; i<groups; i++)
    {
        age_group[i]=0;
    }
    
    graph_erase(buffer, img);
    
    current_date = TIME_IN_DAYS(sim->land->date);
    
    for (j = 0; j < sim->num; j++)
    {
        noble_being * local_being = &(sim->beings[j]);
        n_uint local_dob = being_dob(local_being);
        n_int age_days = current_date - local_dob;
        if (age_days >= max_age) age_days = max_age-1;
        idx = age_days/divisor;
        age_group[idx]++;
        if (age_group[idx] > max) max = age_group[idx];
    }
    max = max * 120/100;
    for (i = 0; i < groups; i++)
    {
        
        x = i*img->x/groups;
        y = img->y-1-(age_group[i]*img->y/max);
        /*graph_line(buffer,img,prev_x,prev_y,x,y,0,0,0,1);*/
        prev_x = x;
        prev_y = y;
    }
    io_free((void*)age_group);
}

/*
 Show a histogram of being heights
 */
void graph_heights(noble_simulation * sim, n_byte * buffer, n_vect2 * img)
{
    const n_int divisor = BEING_MAX_HEIGHT/16;
    n_int groups = BEING_MAX_HEIGHT/divisor;
    n_int * height_group;
    n_int i,idx,max=1;
    n_uint j;
    n_int prev_x = 0, prev_y=img->y-1;
    n_int x,y;
    
    height_group = (n_int*)io_new(groups*sizeof(n_int));
    for (i=0; i<groups; i++)
    {
        height_group[i]=0;
    }
    graph_erase(buffer, img);
    
    for (j = 0; j < sim->num; j++)
    {
        noble_being * local_being = &(sim->beings[j]);
        idx = GET_H(local_being)/divisor;
        height_group[idx]++;
        if (height_group[idx] > max) max = height_group[idx];
    }
    max = max * 120/100;
    for (i = 0; i < groups; i++)
    {
        x = i*img->x/groups;
        y = img->y-1-(height_group[i]*img->y/max);
        /*graph_line(buffer,img,prev_x,prev_y,x,y,0,0,0,1);*/
        prev_x = x;
        prev_y = y;
    }
    
    io_free((void*)height_group);
}

/* return the braincode standard deviation */

static n_uint braincode_standard_deviation(noble_simulation * sim, noble_being * local_being)
{
    n_uint sd = 0;
#ifdef BRAINCODE_ON
    n_int i,av=0,diff;
    
    for (i=0; i<BRAINCODE_SIZE; i++)
    {
        av += being_braincode_internal(local_being)[i];
        av += being_braincode_external(local_being)[i];
    }
    av /= (BRAINCODE_SIZE*2);
    
    for (i=0; i<BRAINCODE_SIZE; i++)
    {
        diff = (n_int)(being_braincode_internal(local_being)[i]) - av;
        if (diff<0) diff=-diff;
        sd += (n_uint)(diff);
        diff = (n_int)(being_braincode_external(local_being)[i]) - av;
        if (diff<0) diff=-diff;
        sd += (n_uint)(diff);
    }
#endif
    return sd;
}


/* return coordinates of the braincode system for phase space plot */
static void graph_braincode_coords(noble_simulation * sim, noble_being * local_being, n_uint * x, n_uint * y)
{
    n_int i;
    *x=0;
    *y=0;
    for (i=0; i<BRAINCODE_SIZE/2; i++)
    {
        *x = *x + being_braincode_internal(local_being)[i] +
        being_braincode_external(local_being)[i];
    }
    while (i < BRAINCODE_SIZE) {
        *y = *y + being_braincode_internal(local_being)[i] +
        being_braincode_external(local_being)[i];
        i++;
    }
}

/* return coordinates of the genome for phase space plot */
static void graph_genespace_coords(noble_being * local_being, n_uint * x, n_uint * y)
{
    n_int ch,b;
    n_genetics * genetics = being_genetics(local_being);
    *x=0;
    *y=0;
    for (ch=0; ch<CHROMOSOMES; ch++)
    {
        for (b=0; b<8; b++)
        {
            *x = *x + ((genetics[ch]>>(b*2))&3);
            *y = *y + ((genetics[ch]>>(16+(b*2)))&3);
        }
    }
}


static void graph_phasespace_dots(noble_simulation * sim, n_byte update_type, n_byte * buffer, n_int img_width, n_int img_height, n_byte graph_type)
{
#ifdef PARASITES_ON
    n_uint i,x=0,y=0;
    n_int min_x,max_x,min_y,max_y;
    n_int dx,dy;
    n_int av_x=0, av_y=0, av_dx=0, av_dy=0;
    const n_int min_variance = 8;
    scope s;
    unsigned int intensity_percent = 100;
    unsigned int grid_horizontal = 10;
    unsigned int grid_vertical = 8;
    
    if (sim->num == 0) {
        /* clear the image */
        n_vect2 img;
        img.x = img_width;
        img.y = img_height;
        graph_erase(buffer, &img);
        return;
    }
    
    s = create_scope((unsigned int)1);
    s.time_ms = (unsigned int)(sim->num);
    s.noise = 0.1;
    s.mode = PHOSPHENE_MODE_POINTS;
    
    if (update_type == PHOSPHENE_DRAW_BACKGROUND) {
        scope_draw(&s, update_type, intensity_percent,
                   grid_horizontal, grid_vertical,
                   (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
        return;
    }
    
    
    for (i=0; i<sim->num; i++)
    {
        switch(graph_type)
        {
            case 0:
                graph_braincode_coords(sim, &(sim->beings[i]), &x, &y);
                break;
            case 1:
                graph_genespace_coords(&(sim->beings[i]), &x, &y);
                break;
        }
        av_x += (n_int)x;
        av_y += (n_int)y;
    }
    av_x /= (n_int)sim->num;
    av_y /= (n_int)sim->num;
    
    for (i=0; i<sim->num; i++)
    {
        switch(graph_type)
        {
            case 0:
                graph_braincode_coords(sim, &(sim->beings[i]), &x, &y);
                break;
            case 1:
                graph_genespace_coords(&(sim->beings[i]), &x, &y);
                break;
        }
        
        dx = (n_int)x - av_x;
        if (dx < 0) dx = -dx;
        
        dy = (n_int)y - av_y;
        if (dy < 0) dy = -dy;
        
        av_dx += dx;
        av_dy += dy;
    }
    av_dx /= (n_int)sim->num;
    av_dy /= (n_int)sim->num;
    
    if (av_dx < min_variance) av_dx = min_variance;
    if (av_dy < min_variance) av_dy = min_variance;
    
    min_x = av_x - av_dx;
    max_x = av_x + av_dx;
    min_y = av_y - av_dy;
    max_y = av_y + av_dy;
    
    for (i=0; i<sim->num; i++)
    {
        switch(graph_type)
        {
            case 0:
                graph_braincode_coords(sim, &(sim->beings[i]), &x, &y);
                break;
            case 1:
                graph_genespace_coords(&(sim->beings[i]), &x, &y);
                break;
        }
        
        scope_update(&s, 0, x, min_x, max_x, (unsigned int)i);
        scope_update(&s, 1, y, min_y, max_y, (unsigned int)i);
    }
    
    scope_draw(&s, update_type, intensity_percent,
               grid_horizontal, grid_vertical,
               (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
#endif
}

static void graph_phasespace_density(noble_simulation * sim, n_byte * buffer, n_int img_width, n_int img_height, n_byte graph_type)
{
#ifdef BRAINCODE_ON
    const n_uint grid = 32;
    n_uint x=0,y=0;
    n_uint n0,n1,tx,ty,bx,by,xx,yy;
    n_uint density[32*32],max=1;
    n_byte r,b;
    n_int i;
    n_uint j;
    for (i = 0; i < img_width*img_height*3; i+=3)
    {
        n_vect2 img;
        img.x = img_width;
        img.y = img_height;
        graph_erase(buffer, &img);
    }
    
    io_erase((n_byte*)density, sizeof(n_uint)*32*32);
    
    for (j=0; j < sim->num; j++)
    {
        switch(graph_type)
        {
            case 0:
                graph_braincode_coords(sim, &(sim->beings[j]), &x, &y);
                x = x * (grid-1) / (256*BRAINCODE_SIZE);
                y = (grid-1) - (y * (grid-1) / (255*BRAINCODE_SIZE));
                break;
            case 1:
                graph_genespace_coords(&(sim->beings[j]), &x, &y);
                x = x * (grid-1) / (4*8*CHROMOSOMES);
                y = (grid-1) - (y * (grid-1) / (4*8*CHROMOSOMES));
                break;
        }
        density[y*grid+x]++;
    }
    
    for (j = 0; j < grid*grid; j++)
    {
        if (density[j] > max) max = density[j];
    }
    
    n0 = 0;
    for (y = 0; y < grid; y++)
    {
        for (x = 0; x < grid; x++, n0++)
        {
            if (density[n0] > 0)
            {
                r = (n_byte)(density[n0]*255/max);
                b = 255-r;
                
                tx = x * img_width / grid;
                bx = (x+1) * img_width / grid;
                ty = y * img_height / grid;
                by = (y+1) * img_height / grid;
                
                for (yy = ty; yy < by; yy++)
                {
                    n1 = ((yy*img_width)+tx)*3;
                    for (xx = tx; xx < bx; xx++, n1+=3)
                    {
                        buffer[n1] = r;
                        buffer[n1+1] = 0;
                        buffer[n1+2] = b;
                    }
                }
            }
        }
    }
#endif
}

void graph_socialsim(noble_simulation * sim, n_byte update_type, n_byte * buffer, n_int img_width, n_int img_height)
{
    n_uint i;
    n_int min_x=65536, max_x=-1, min_y=65536, max_y=-1;
    
    noble_being * local_being;
    scope s;
    unsigned int intensity_percent = 100;
    unsigned int grid_horizontal = 10;
    unsigned int grid_vertical = 10;
    
    if (sim->num == 0) {
        /* clear the image */
        n_vect2 img;
        img.x = img_width;
        img.y = img_height;
        graph_erase(buffer, &img);
        return;
    }
    
    s = create_scope((unsigned int)1);
    s.time_ms = (unsigned int)(sim->num);
    s.noise = 0.1;
    s.mode = PHOSPHENE_MODE_POINTS;
    
    if (update_type == PHOSPHENE_DRAW_BACKGROUND) {
        scope_draw(&s, update_type, intensity_percent,
                   grid_horizontal, grid_vertical,
                   (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
        return;
    }
    
    for (i = 0; i < sim->num; i++)
    {
        n_int coord_x, coord_y;
        
        local_being = &(sim->beings[i]);
        
        coord_x = local_being->social_x;
        coord_y = local_being->social_y;
        
        if (coord_x < min_x)
        {
            min_x = coord_x;
        }
        else if (coord_x > max_x)
        {
            max_x = coord_x;
        }
        if (coord_y < min_y)
        {
            min_y = coord_y;
        }
        else if (coord_y > max_y)
        {
            max_y = coord_y;
        }
    }
    
    if ((max_x <= min_x) || (max_y <= min_y)) return;
    
    for (i = 0; i < sim->num; i++)
    {
        int x, y;
        
        local_being = &(sim->beings[i]);
        
        x = local_being->social_x;
        y = local_being->social_y;
        
        scope_update(&s, 0, x, min_x, max_x, (unsigned int)i);
        scope_update(&s, 1, y, min_y, max_y, (unsigned int)i);
    }
    
    scope_draw(&s, update_type, intensity_percent,
               grid_horizontal, grid_vertical,
               (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
}

/* plot the places where beings met */
void graph_meet_places(noble_simulation * sim, n_byte update_type, n_byte * buffer, n_int img_width, n_int img_height)
{
    n_uint i,index,ctr=0;
    /** dimensions of APESPACE */
    n_int min_x=0, max_x=65535, min_y=0, max_y=65535;
    scope s;
    unsigned int intensity_percent = 100;
    unsigned int grid_horizontal = 10;
    unsigned int grid_vertical = 10;
    
    if (sim->num == 0) {
        /* clear the image */
        n_vect2 img;
        img.x = img_width;
        img.y = img_height;
        graph_erase(buffer, &img);
        return;
    }
    
    s = create_scope((unsigned int)1);
    s.time_ms = (unsigned int)(sim->num);
    s.noise = 0.1;
    s.mode = PHOSPHENE_MODE_POINTS;
    
    if (update_type == PHOSPHENE_DRAW_BACKGROUND) {
        scope_draw(&s, update_type, intensity_percent,
                   grid_horizontal, grid_vertical,
                   (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
        return;
    }
    
    /** count the number of locations */
    for (i = 0; i < sim->num; i++)
    {
		social_link * graph = being_social(&(sim->beings[i]));
        
        /** for each non-self social graph entry */
        for (index = 1; index < SOCIAL_SIZE_BEINGS; index++)
		{
			if (!SOCIAL_GRAPH_ENTRY_EMPTY(graph,index))
			{
                if (SOCIAL_GRAPH_ENTRY_LOCATION_EXISTS(graph,index))
                {
                    ctr++;
                }
			}
		}
    }
    
    if (ctr > 0) {
        s.time_ms = (unsigned int)ctr;
    }
    
    ctr = 0;
    for (i = 0; i < sim->num; i++)
    {
		social_link * graph = being_social(&(sim->beings[i]));
        
        /** for each non-self social graph entry */
        for (index = 1; index < SOCIAL_SIZE_BEINGS; index++)
		{
			if (!SOCIAL_GRAPH_ENTRY_EMPTY(graph,index))
			{
                if (SOCIAL_GRAPH_ENTRY_LOCATION_EXISTS(graph,index))
                {
                    scope_update(&s, 0, (int)graph[index].location[0],
                                 (int)min_x, (int)max_x, (unsigned int)ctr);
                    scope_update(&s, 1, (int)graph[index].location[1],
                                 (int)min_y, (int)max_y, (unsigned int)ctr);
                    ctr++;
                }
            }
		}
    }
    
    scope_draw(&s, update_type, intensity_percent,
               grid_horizontal, grid_vertical,
               (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
}

void graph_phasespace(noble_simulation * sim, n_byte update_type, n_byte * buffer, n_int img_width, n_int img_height, n_byte graph_type, n_byte data_type)
{
    if (graph_type==0)
    {
        graph_phasespace_dots(sim, update_type, buffer, img_width, img_height,data_type);
    }
    else
    {
        graph_phasespace_density(sim, buffer, img_width, img_height,data_type);
    }
}

/*
 Displays the braincode program for an individual
 */
void graph_braincode(noble_simulation * sim, noble_being * local_being, n_byte * buffer, n_int img_width, n_int img_height, n_byte clear)
{
#ifdef BRAINCODE_ON
#ifdef PARASITES_ON
    n_int i,x,y,n,half_width;
    n_byte * code;
    
    if (local_being!=0)
    {
        /* clear the image */
        if (clear!=0) for (i = 0; i < img_width*img_height*3; i++) buffer[i]=0;
        
        half_width = img_width/2;
        y = sim->land->time % img_height;
        n = y*img_width*3;
        for (x = 0; x < img_width; x++, n+=3)
        {
            if (x<half_width)
            {
                i = (x * ((BRAINCODE_SIZE/BRAINCODE_BYTES_PER_INSTRUCTION)-1) / half_width)*BRAINCODE_BYTES_PER_INSTRUCTION;
                code = being_braincode_internal(local_being);
            }
            else
            {
                i = ((x-half_width) * ((BRAINCODE_SIZE/BRAINCODE_BYTES_PER_INSTRUCTION)-1) / half_width)*BRAINCODE_BYTES_PER_INSTRUCTION;
                code = being_braincode_external(local_being);
            }
            buffer[n] = code[i];
            buffer[n+1] = code[i+1];
            buffer[n+2] = code[i+2];
        }
    }
    
#endif
#endif
}

/*
 Displays the preferences of the population
 */
void graph_preferences(noble_simulation * sim, n_byte update_type, n_byte * buffer, n_int img_width, n_int img_height)
{
    n_uint i;
    n_int p,x=0,y=0,half=PREFERENCES/2;
    n_int min_x=0, max_x=0, min_y=0, max_y=0;
    n_int av_x=0, av_y=0, av_dx=0, av_dy=0,dx,dy;
    const n_int dimension = 10000;
    const n_int min_variance = 20;
    noble_being * local_being;
    scope s;
    unsigned int intensity_percent = 100;
    unsigned int grid_horizontal = 10;
    unsigned int grid_vertical = 8;
    
    if (sim->num == 0) {
        /* clear the image */
        n_vect2 img;
        img.x = img_width;
        img.y = img_height;
        graph_erase(buffer, &img);
        
        return;
    }
    
    s = create_scope((unsigned int)1);
    s.time_ms = (unsigned int)(sim->num);
    s.noise = 0.1;
    s.mode = PHOSPHENE_MODE_POINTS;
    
    if (update_type == PHOSPHENE_DRAW_BACKGROUND) {
        scope_draw(&s, update_type, intensity_percent,
                   grid_horizontal, grid_vertical,
                   (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
        return;
    }
    
    for (i = 0; i < sim->num; i++)
    {
        local_being = &(sim->beings[i]);
        x = 0;
        for (p = 0; p < half; p++)
        {
            x += local_being->learned_preference[p];
        }
        av_x += x * dimension / (half*255);
        
        y = 0;
        while (p < PREFERENCES)
        {
            y += local_being->learned_preference[p];
            p++;
        }
        av_y += y * dimension / ((PREFERENCES-half)*255);
    }
    av_x /= (n_int)sim->num;
    av_y /= (n_int)sim->num;
    
    for (i = 0; i < sim->num; i++)
    {
        local_being = &(sim->beings[i]);
        x = 0;
        for (p = 0; p < half; p++)
        {
            x += local_being->learned_preference[p];
        }
        dx = (x * dimension / (half*255)) - av_x;
        if (dx < 0) dx = -dx;
        av_dx += dx;
        
        y = 0;
        while (p < PREFERENCES)
        {
            y += local_being->learned_preference[p];
            p++;
        }
        dy = (y * dimension / ((PREFERENCES-half)*255))-av_y;
        if (dy < 0) dy = -dy;
        av_dy += dy;
    }
    av_dx /= (n_int)sim->num;
    av_dy /= (n_int)sim->num;
    if (av_dx < min_variance) av_dx = min_variance;
    if (av_dy < min_variance) av_dy = min_variance;
    
    min_x = av_x - av_dx*2 - 1;
    max_x = av_x + av_dx*2 + 1;
    min_y = av_y - av_dy*2 - 1;
    max_y = av_y + av_dy*2 + 1;
    
    if ((max_x <= min_x) || (max_y <= min_y)) return;
    
    for (i = 0; i < sim->num; i++)
    {
        local_being = &(sim->beings[i]);
        x = 0;
        for (p = 0; p < half; p++)
        {
            x += local_being->learned_preference[p];
        }
        x = x * dimension / (half*255);
        scope_update(&s, 0, x, min_x, max_x, (unsigned int)i);
        
        y = 0;
        while (p < PREFERENCES)
        {
            y += local_being->learned_preference[p];
            p++;
        }
        y = y * dimension / ((PREFERENCES-half)*255);
        scope_update(&s, 1, y, min_y, max_y, (unsigned int)i);
    }
    
    scope_draw(&s, update_type, intensity_percent,
               grid_horizontal, grid_vertical,
               (unsigned char*)buffer, (unsigned int)img_width, (unsigned int)img_height);
    
}

void graph_command(n_int gc_val)
{
    if (gc_val == GC_CLEAR_BRAINCODE)
    {
        graph_clear = 1;
        graph_state = GC_BRAINCODE;
    }
    else
    {
        graph_state = gc_val;
    }
}

void  graph_draw(noble_simulation * local_sim, n_byte * graph, n_int dim_x, n_int dim_y)
{    
    switch (graph_state)
    {
        case GC_IDEOSPHERE:
            graph_ideosphere(local_sim, graph, dim_x, dim_y);
            break;
        case GC_BRAINCODE:
            if (local_sim->select != NO_BEINGS_FOUND)
            {
                graph_braincode(local_sim, &local_sim->beings[local_sim->select],graph, dim_x, dim_y, graph_clear);
                graph_clear = 0;
            }
            break;
        case GC_HONOR:
            graph_honor_distribution(local_sim, PHOSPHENE_DRAW_ALL, graph, dim_x, dim_y);
            break;
        case GC_PATHOGENS:
            graph_pathogens(local_sim, PHOSPHENE_DRAW_ALL, graph, dim_x, dim_y);
            break;
        case GC_RELATIONSHIPS:
            graph_relationship_matrix(local_sim, graph, dim_x, dim_y);
            break;
        case GC_GENEPOOL:
            graph_genepool(local_sim, graph, dim_x, dim_y);
            break;
        case GC_PREFERENCES:
            graph_preferences(local_sim, PHOSPHENE_DRAW_ALL, graph, dim_x, dim_y);
            break;
        case GC_SOCIALSIM:
            graph_socialsim(local_sim, PHOSPHENE_DRAW_ALL, graph, dim_x, dim_y);
            break;
        case GC_MEET_PLACES:
            graph_meet_places(local_sim, PHOSPHENE_DRAW_ALL, graph, dim_x, dim_y);
            break;
        case GC_PHASESPACE:
            graph_phasespace(local_sim, PHOSPHENE_DRAW_ALL, graph, dim_x, dim_y, 0, 0);
            break;
        case GC_VASCULAR:
        default:
            if (local_sim->select != NO_BEINGS_FOUND)
            {
                /** set this to a non-zero value to show key points on the skeleton
                 which may be useful for debugging */
                
                /* doesn't work with Ofast optimization */
                n_vect2 dim;
                n_byte show_skeleton_keypoints = 0;

                dim.x = dim_x;
                dim.y = dim_y;
                
                graph_vascular(&local_sim->beings[local_sim->select].genes, graph,
                               &dim,
                               dim_x*10/100,dim_y*10/100,
                               dim_x*40/100,dim_y*90/100,
                               1, 1,
                               30, 0, 20, 20, 0,
                               show_skeleton_keypoints);
                
            }
            break;
    }
}


