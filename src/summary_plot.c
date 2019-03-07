/* Copyright (C) 2006  Leonardo de Oliveira Martins
 *
 * leo at lbm ab a u-tokyo ac jp
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*! \file summary_plot.c
 * \brief Plots summary info in PDF format and makes a report
 *
 * Uses cairo graphics library to plot credible intervals, breakpoint locations, mosaic structures and
 * \f$\hat{d}_{SPR}\f$ heatmap.
 * Also prints a report to file.
 * 
 */

#include "summary.h"

/* remember that in cairo the origin (0,0) is on the superior left corner (thus y-axis is inverted) */ 

#define canvas_width 460.  /* in points; 1 inch = 72 points */
#define canvas_height 200. /* actual height is 2*canvas_height, we have two plots per page */
#define canvas_border 15. 

double x_range = canvas_width - 3. * canvas_border;
double y_range = canvas_height - 3. * canvas_border;

void show_text_centered_h (cairo_t *cr, char *utf8, double x, double y);
void show_text_centered_v (cairo_t *cr, char *utf8, double x, double y);
void plot_xy_margins (cairo_t *cr, double x, double y);
void plot_xy_distribution (cairo_t *cr, int *f2, int max_dist, empfreq freq, int sum_freq, int mode, summary sm);

void plot_y_margin (cairo_t *cr, double x0, double y0, double x1, double y1, double min_y, double max_y);
void plot_x_margin (cairo_t *cr, double x0, double y0, double x1, double y1, double min_x, double max_x);
void plot_histogram (cairo_t *cr, double x0, double y0, double x1, double y1, double *xv, double *yv, int n);

void
plot_post_recomb_freq (summary sm)
{
  cairo_surface_t *surface;
  cairo_t *cr;
	int i, j, sum_freq = 0, max_freq = 0, max_dist = 0, f2[sm->n_segments-1], mode;
  empfreq freq;
  char caption[128];

	surface = cairo_pdf_surface_create ("recomb_freq.pdf", canvas_width, 2 * canvas_height);
  cr = cairo_create (surface);
  cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size (cr, 8.0);

  mode = find_mode_int_weighted (sm->nCOP, sm->n_samples, sm->w_sample);

  // all samples
  for (i=0; i < sm->n_segments-1; i++) f2[i] = 0;
  freq = new_empfreq (sm->n_segments-1);

	for (j=0; j < sm->n_samples; j++)
		for (i=0; i < sm->n_segments-1; i++)
			if (sm->dSPR[i][j]) {
        /* frequency of recombination */
        freq->i[i].freq += sm->w_sample[j];
        sum_freq += sm->w_sample[j];
				if (freq->i[i].freq > max_freq) max_freq = freq->i[i].freq;
        /* recombination distance */
        f2[i] += sm->dSPR[i][j] * sm->w_sample[j];
        if (f2[i] > max_dist) max_dist = f2[i];
			}

  sprintf (caption, "all samples");
  show_text_centered_h (cr, caption, canvas_width/2., canvas_border/2.);
  plot_xy_margins (cr, (double)(sm->n_sites), (double)(max_dist)/(double)(sm->n_samples_original));
  plot_xy_distribution (cr, f2, max_dist, freq, sum_freq, mode, sm);

  // only samples with modal number of breakpoints
  for (i=0; i < sm->n_segments-1; i++) { // re-initialize values
    f2[i] = 0;
    freq->i[i].freq = 0;
    freq->i[i].idx=i;
  }
  sum_freq = max_dist = 0;

  max_freq = 0; // max_freq is actually the number of samples 
  for (j=0; j < sm->n_samples; j++) if (sm->nCOP[j] ==mode) {
    max_freq++; 
    for (i=0; i < sm->n_segments-1; i++)
      if (sm->dSPR[i][j]) {
        /* frequency of recombination */
        freq->i[i].freq += sm->w_sample[j];
        sum_freq += sm->w_sample[j];
        /* recombination distance */
        f2[i] += sm->dSPR[i][j] * sm->w_sample[j];
        if (f2[i] > max_dist) max_dist = f2[i];
      }
  }

  cairo_translate(cr, 0, canvas_height);

  sprintf (caption, "samples with modal ( %d ) number of breakpoints", mode);
  show_text_centered_h (cr, caption, canvas_width/2., canvas_border/2.);
  plot_xy_margins (cr, (double)(sm->n_sites), (double)(max_dist)/(double)(max_freq));
  plot_xy_distribution (cr, f2, max_dist, freq, sum_freq, mode, sm);

  del_empfreq (freq);
  cairo_show_page (cr);
  cairo_surface_destroy (surface);
  cairo_destroy (cr);
}


void
show_text_centered_h (cairo_t *cr, char *utf8, double x, double y)
{
  cairo_text_extents_t extents;

  cairo_text_extents (cr, utf8, &extents);
  x = x - (extents.width/2 + extents.x_bearing);
  y = y - (extents.height/2 + extents.y_bearing);
  cairo_move_to (cr, x, y);
  cairo_show_text (cr, utf8);
}

void
show_text_centered_v (cairo_t *cr, char *utf8, double x, double y)
{
  cairo_text_extents_t extents;

  cairo_save (cr);
  cairo_text_extents (cr, utf8, &extents);
  y = y + (extents.width/2 + extents.x_bearing);
  x = x - (extents.height/2 + extents.y_bearing);
  cairo_move_to (cr, x, y);
  cairo_rotate (cr, 3.1415926535 * 1.5);
  cairo_show_text (cr, utf8);
  cairo_restore (cr);
}

void
plot_xy_margins (cairo_t *cr, double x, double y)
{
  double vertical_x = 1.6 * canvas_border, 
         vertical_x_label = 0.8 * canvas_border,
         horiz_y = 1.4 * canvas_border + y_range,
         horiz_y_label = 2.2 * canvas_border + y_range,
         tick, frac;
  char label[32];


  cairo_save (cr);
	/* horizontal line (left-right) */
	cairo_move_to (cr, 2 * canvas_border, horiz_y);
  cairo_rel_line_to (cr, x_range, 0);
  /* tick marks for x axis */
  for (tick = 2 * canvas_border; tick <= 2 * canvas_border + x_range;  tick += x_range/4.) {
    cairo_move_to (cr, tick, horiz_y);
    cairo_rel_line_to (cr, 0, 0.2 * canvas_border);
  }
	/* vertical line ( top-bottom) */
	cairo_move_to (cr, vertical_x, canvas_border);
  cairo_rel_line_to (cr, 0, y_range);
  /* tick marks for y axis */
  for (tick = canvas_border; tick <= canvas_border + y_range;  tick += y_range/4.) {
    cairo_move_to (cr, vertical_x, tick);
    cairo_rel_line_to (cr, - 0.2 * canvas_border, 0);
  }
	cairo_set_line_width (cr, 0.1);
	cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
	cairo_stroke (cr);
 
  /* X labels */
  for (frac = 0.; frac <= 1.; frac += 0.25) {
    sprintf (label, "%.0f", frac * x );
    tick = 2 * canvas_border + frac * x_range;
    show_text_centered_h (cr, label, tick, horiz_y_label);
  }

  /* Y labels */
  for (frac = 0.; frac <= 1.; frac += 0.25) {
    sprintf (label, "%.3f", (1. - frac) * y );
    tick = canvas_border + frac * y_range;
    show_text_centered_h (cr, label, vertical_x_label, tick);
  }
  cairo_restore (cr);
}

void
plot_xy_distribution (cairo_t *cr, int *f2, int max_dist, empfreq freq, int sum_freq, int mode, summary sm)
{
  int i, partial_sum, q1[mode], q2[mode], q3[mode], i1=0, i2=0, i3=0;
  double rel_x, rel_y, shift_y, step_sum = (double)(sum_freq)/(double)(mode), 
         q1_sum = 0.025 * step_sum, 
         q2_sum = 0.5   * step_sum, 
         q3_sum = 0.975 * step_sum;
  char label[32];

  /* calculate centered credible intervals and median per breakpoint region */
  partial_sum = freq->i[0].freq;
  for(i=0; i < sm->n_segments-2; i++) {
    if (partial_sum > q1_sum) { q1[i1++] = freq->i[i].idx; q1_sum += step_sum; }
    if (partial_sum > q2_sum) { q2[i2++] = freq->i[i].idx; q2_sum += step_sum; }
    if (partial_sum > q3_sum) { q3[i3++] = freq->i[i].idx; q3_sum += step_sum; }
    partial_sum += freq->i[i+1].freq;
  }
  /* just in case (if we have some weird roundoff error or I made some mistake) */
  if (i1 < mode) for (i=i1; i<mode; i++) q1[i] = freq->i[sm->n_segments-2].idx;
  if (i2 < mode) for (i=i2; i<mode; i++) q2[i] = freq->i[sm->n_segments-2].idx;
  if (i3 < mode) for (i=i3; i<mode; i++) q3[i] = freq->i[sm->n_segments-2].idx;

  cairo_save (cr);
	/* recomb distance*/
	cairo_set_line_width (cr, x_range/(double)(sm->n_sites));
	cairo_set_source_rgb (cr, 0.9, 0.7, 0.7);
	for(i=0; i < sm->n_segments-1; i++) {
		rel_x = (double)(x_range * sm->breakpoint[i])/(double)(sm->n_sites);
    rel_y = (double)(y_range * f2[i])/(double)(max_dist);
		cairo_move_to (cr, 2. * canvas_border + rel_x, y_range + canvas_border);
		cairo_rel_line_to (cr, 0, -rel_y);
	}
	cairo_stroke (cr);

	/* coldspots */
	cairo_set_source_rgb (cr, 0.6, 0.6, 1.0);
	for(i=0; i < sm->n_coldspot; i++) {
		rel_x = (double)(x_range * sm->coldspot[i])/(double)(sm->n_sites);
    cairo_move_to (cr, 2. * canvas_border + rel_x, y_range + 1.05 * canvas_border);
		cairo_rel_line_to (cr, 0, 0.3 * canvas_border);
	}
	cairo_stroke (cr);

  /* credible intervals */
  sort_empfreq_decreasing (freq); /* sort from largest to smallest freq */
	cairo_set_source_rgb (cr, 1.0, 0.2, 0.2);
  partial_sum = freq->i[0].freq;
  for(i=0; (partial_sum < (int)(0.95 * (double)(sum_freq))) && (i < sm->n_segments-1); i++) {
    rel_x = (double)(x_range * sm->breakpoint[ freq->i[i].idx ])/(double)(sm->n_sites);
    cairo_move_to (cr, 2. * canvas_border + rel_x, y_range + 1.05 * canvas_border);
    cairo_rel_line_to (cr, 0, 0.3 * canvas_border);

    if (i<sm->n_segments-2) partial_sum += freq->i[i+1].freq;
  }
	cairo_stroke (cr);

  /* mode-based centered quantiles (median and credible intervals) */
  for(i=0; i < mode; i++) {
    shift_y = ((double)(i%4)/4.) * 0.5 * y_range;
    /* credible intervals as transparent rectangles */
    rel_x = (double)(x_range * sm->breakpoint[ q1[i] ])/(double)(sm->n_sites); // start site
    rel_y = (double)(x_range * sm->breakpoint[ q3[i] ])/(double)(sm->n_sites); // end site (note that this is _not_ Y)
    cairo_rectangle (cr, 2 * canvas_border + rel_x, canvas_border + shift_y, rel_y - rel_x, 0.125 * y_range);
    cairo_set_source_rgba (cr, 1., 0., 0., 0.6); // alpha channel (transparency)
    cairo_fill(cr);
    /* median values as vertical bars */
    rel_x = (double)(x_range * sm->breakpoint[ q2[i] ])/(double)(sm->n_sites);
    cairo_move_to (cr, 2 * canvas_border + rel_x, canvas_border + shift_y);
    cairo_rel_line_to (cr, 0., 0.125 * y_range);
    cairo_set_source_rgb (cr, 0.7, 0., 0.);
    cairo_stroke (cr);
  }
  for(i=0; i < mode; i++) {
    rel_x = (double)(x_range * sm->breakpoint[ q2[i] ])/(double)(sm->n_sites);
    sprintf (label, "%d", sm->breakpoint [ q2[i] ]);
    /* 4 levels between 0 and 0.25*y_range */
    shift_y = ((double)(i%4)/4.+0.25) * 0.5 * y_range;
    cairo_move_to (cr, 2 * canvas_border + rel_x, canvas_border + shift_y);
    cairo_set_source_rgb (cr, 0., 0., 0.);
    cairo_show_text (cr, label);
  }

  cairo_restore (cr);
}

void
plot_empirical_frequency (empfreq e, char *filename, char *title, double width, double height, double border)
{
  cairo_surface_t *surface;
  cairo_t *cr;
  int i, sum_freq = 0;
  double x[e->n], y[e->n];
  char file[strlen(filename+5)];

  sprintf (file, "%s.pdf", filename);
  surface = cairo_pdf_surface_create (file, width, height);
  cr = cairo_create (surface);
  cairo_select_font_face (cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size (cr, 8.0);

  for (i=0; i<e->n; i++) sum_freq += e->i[i].freq;
  for (i=0; i<e->n; i++) {
    x[i] = (double)(e->i[i].idx - e->min)/(double)(e->max - e->min);
    y[i] = (double)(e->i[i].freq)/(double)(e->i[0].freq); /* empfreq is ordered from highest to lowest frequency */
  }

  show_text_centered_h (cr, title, width/2., border/2.);
  plot_y_margin (cr, 0, border, 1.6*border, height- 2*border, 0., (double)(e->i[0].freq)/(double)(sum_freq));
  plot_x_margin (cr, 2*border, height - 1.6*border, width - border, height, (double)(e->min), (double)(e->max));
  plot_histogram (cr, 2*border, border, width - border, height - 2*border, x, y, e->n);

  cairo_show_page (cr);
  cairo_surface_destroy (surface);
  cairo_destroy (cr);
}

void
plot_y_margin (cairo_t *cr, double x0, double y0, double x1, double y1, double min_y, double max_y)
{
  double tick, frac;
  char label[32];

  cairo_save (cr);
	/* vertical line ( top-bottom) */
  cairo_move_to (cr, x1, y0);
  cairo_line_to (cr, x1, y1);
  /* tick marks for y axis */
  for (tick = y0; tick <= y1;  tick += (y1-y0)/4.) {
    cairo_move_to (cr, x1, tick);
    cairo_rel_line_to (cr, - 0.25 * (x1-x0), 0);
  }
	cairo_set_line_width (cr, 0.1);
	cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
	cairo_stroke (cr);

  /* Y labels */
  for (frac = 0.; frac <= 1.; frac += 0.25) {
    sprintf (label, "%.3f", (1. - frac) * (max_y - min_y) );
    tick = y0 + frac * (y1-y0);
    show_text_centered_h (cr, label, 0.5 * (x1+x0), tick);
  }
  cairo_restore (cr);
}

void
plot_x_margin (cairo_t *cr, double x0, double y0, double x1, double y1, double min_x, double max_x)
{
  double tick, frac;
  char label[32];

  cairo_save (cr);
	/* horizontal line (left-right) */
  cairo_move_to (cr, x0, y0);
  cairo_line_to (cr, x1, y0);
  /* tick marks for x axis */
  for (tick = x0; tick <= x1;  tick += (x1-x0)/4.) {
    cairo_move_to (cr, tick, y0);
    cairo_rel_line_to (cr, 0, 0.25 * (y1-y0));
  }
	cairo_set_line_width (cr, 0.1);
	cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
	cairo_stroke (cr);
 
  /* X labels */
  for (frac = 0.; frac <= 1.; frac += 0.25) {
    sprintf (label, "%.0f", min_x + frac * (max_x - min_x));
    tick = x0 + frac * (x1-x0);
    show_text_centered_h (cr, label, tick, 0.5 * (y1+y0));
  }
  cairo_restore (cr);
}

void
plot_histogram (cairo_t *cr, double x0, double y0, double x1, double y1, double *xv, double *yv, int n)
{
  int i;

  cairo_set_line_width (cr, (x1-x0)/(double)(n));
  cairo_set_source_rgb (cr, 0.9, 0.2, 0.2);
  for(i=0; i < n; i++) {
    cairo_move_to (cr, x0 + xv[i]*(x1-x0), y1);
    cairo_rel_line_to (cr, 0, yv[i]*(y0-y1)); // = -(y1-y0), since move upwards is negative 
  }
  cairo_stroke (cr);

  cairo_restore (cr);
}
