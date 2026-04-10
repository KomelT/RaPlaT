/****************************************************************************
 *
 * MODULE:       r.nr3gpp
 * AUTHOR(S):    RaPlaT community upgrade
 *
 * PURPOSE:      Calculates large-scale path loss for a single transmitter
 *               using simplified 3GPP TR 38.901 formulas (UMa/UMi).
 *
 * COPYRIGHT:    (C) 2026 Jozef Stefan Institute contributors
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <math.h>

/*
 * Simplified 3GPP TR 38.901 large-scale path loss models.
 * - "uma": Urban Macro
 * - "umi": Urban Micro Street Canyon
 *
 * Frequency is in MHz (converted to GHz internally), distance in meters.
 */
static FCELL calc_nr3gpp(double d2d_m, double freq_mhz, double tx_h_m, double rx_h_m, double radius_km, const char *scenario)
{
  FCELL x;
  double fc_ghz, d3d_m, los_db, nlos_db, pl_db;
  double dh_m = tx_h_m - rx_h_m;

  if (d2d_m < 10.0 || d2d_m > radius_km * 1000.0) {
    Rast_set_f_null_value(&x, 1);
    return x;
  }

  fc_ghz = freq_mhz / 1000.0;
  d3d_m = sqrt(d2d_m * d2d_m + dh_m * dh_m);

  if (strcmp(scenario, "uma") == 0) {
    los_db = 28.0 + 22.0 * log10(d3d_m) + 20.0 * log10(fc_ghz);
    nlos_db = 13.54 + 39.08 * log10(d3d_m) + 20.0 * log10(fc_ghz) - 0.6 * (rx_h_m - 1.5);
  } else if (strcmp(scenario, "umi") == 0) {
    los_db = 32.4 + 21.0 * log10(d3d_m) + 20.0 * log10(fc_ghz);
    nlos_db = 22.4 + 35.3 * log10(d3d_m) + 21.3 * log10(fc_ghz) - 0.3 * (rx_h_m - 1.5);
  } else {
    G_fatal_error(_("Unknown scenario: [%s]."), scenario);
  }

  pl_db = (los_db > nlos_db) ? los_db : nlos_db;
  x = (FCELL)pl_db;
  return x;
}

int main(int argc, char *argv[])
{
  double east;
  double north;
  double ant_height, frequency, radius;
  double rec_ant_height;

  struct Cell_head window;
  char *name;
  char *result;
  const char *mapset;
  void *inrast;
  unsigned char *outrast;
  int nrows, ncols;
  int row, col;
  int infd, outfd;
  int verbose;

  struct History history;
  struct GModule *module;

  struct Option *input, *opt_coord, *opt_ant_h, *opt_rx_h, *opt_radius, *opt_freq, *opt_scenario, *output;
  struct Flag *flag_quiet;

  G_gisinit(argv[0]);

  module = G_define_module();
  G_add_keyword(_("raster"));
  G_add_keyword(_("3gpp"));
  G_add_keyword(_("5g"));
  module->description = _("RaPlaT - simplified 3GPP NR path loss module (v10apr2026)");

  input = G_define_standard_option(G_OPT_R_INPUT);
  input->key = "input_dem";
  output = G_define_standard_option(G_OPT_R_OUTPUT);

  flag_quiet = G_define_flag();
  flag_quiet->key = 'q';
  flag_quiet->description = _("Quiet");

  opt_scenario = G_define_option();
  opt_scenario->key = "scenario";
  opt_scenario->type = TYPE_STRING;
  opt_scenario->required = NO;
  opt_scenario->options = "uma,umi";
  opt_scenario->answer = "uma";
  opt_scenario->description = _("Deployment scenario");

  opt_coord = G_define_option();
  opt_coord->key = "coordinate";
  opt_coord->type = TYPE_STRING;
  opt_coord->required = YES;
  opt_coord->key_desc = "x,y";
  opt_coord->description = _("Base station coordinates");

  opt_radius = G_define_option();
  opt_radius->key = "radius";
  opt_radius->type = TYPE_DOUBLE;
  opt_radius->required = NO;
  opt_radius->answer = "10";
  opt_radius->description = _("Computation radius [km]");

  opt_ant_h = G_define_option();
  opt_ant_h->key = "ant_height";
  opt_ant_h->type = TYPE_DOUBLE;
  opt_ant_h->required = NO;
  opt_ant_h->answer = "25";
  opt_ant_h->description = _("Transmitter antenna height [m]");

  opt_rx_h = G_define_option();
  opt_rx_h->key = "rx_ant_height";
  opt_rx_h->type = TYPE_DOUBLE;
  opt_rx_h->required = NO;
  opt_rx_h->answer = "1.5";
  opt_rx_h->description = _("Receiver antenna height [m]");

  opt_freq = G_define_option();
  opt_freq->key = "frequency";
  opt_freq->type = TYPE_DOUBLE;
  opt_freq->required = YES;
  opt_freq->description = _("Frequency [MHz]");

  if (G_parser(argc, argv))
    exit(EXIT_FAILURE);

  name = input->answer;
  result = output->answer;
  verbose = !flag_quiet->answer;
  G_scan_easting(opt_coord->answers[0], &east, G_projection());
  G_scan_northing(opt_coord->answers[1], &north, G_projection());
  sscanf(opt_ant_h->answer, "%lf", &ant_height);
  sscanf(opt_radius->answer, "%lf", &radius);
  sscanf(opt_freq->answer, "%lf", &frequency);
  sscanf(opt_rx_h->answer, "%lf", &rec_ant_height);

  if (frequency <= 0.0)
    G_fatal_error(_("Frequency must be positive"));
  if (ant_height <= 0.0 || rec_ant_height <= 0.0)
    G_fatal_error(_("Antenna heights must be positive"));
  if (radius <= 0.0)
    G_fatal_error(_("Radius must be positive"));

  mapset = G_find_raster(name, "");
  if (mapset == NULL)
    G_fatal_error(_("Raster map <%s> not found"), name);

  if (G_legal_filename(result) < 0)
    G_fatal_error(_("<%s> is an illegal file name"), result);

  if ((infd = Rast_open_old(name, mapset)) < 0)
    G_fatal_error(_("Unable to open raster map <%s>"), name);

  G_get_window(&window);

  if (east < window.west || east > window.east || north > window.north || north < window.south)
    G_fatal_error(_("Specified base station coordinates are outside current region bounds."));

  inrast = Rast_allocate_buf(FCELL_TYPE);
  nrows = Rast_window_rows();
  ncols = Rast_window_cols();
  outrast = Rast_allocate_buf(FCELL_TYPE);

  if ((outfd = Rast_open_new(result, FCELL_TYPE)) < 0)
    G_fatal_error(_("Unable to create raster map <%s>"), result);

  for (row = 0; row < nrows; row++) {
    FCELL f_in, f_out;
    double rec_east, rec_north, dist_2d;

    if (verbose)
      G_percent(row, nrows, 10);

    Rast_get_row(infd, inrast, row, FCELL_TYPE);

    for (col = 0; col < ncols; col++) {
      f_in = ((FCELL *)inrast)[col];
      if (isnan((double)f_in)) {
        Rast_set_f_null_value(&f_out, 1);
        ((FCELL *)outrast)[col] = f_out;
        continue;
      }

      rec_east = window.west + window.ew_res / 2.0 + (col * window.ew_res);
      rec_north = window.north - window.ns_res / 2.0 - (row * window.ns_res);
      dist_2d = sqrt(pow((east - rec_east), 2) + pow((north - rec_north), 2));

      f_out = calc_nr3gpp(dist_2d, frequency, ant_height, rec_ant_height, radius, opt_scenario->answer);
      ((FCELL *)outrast)[col] = f_out;
    }

    Rast_put_row(outfd, outrast, FCELL_TYPE);
  }

  G_free(inrast);
  G_free(outrast);
  Rast_close(infd);
  Rast_close(outfd);

  Rast_short_history(result, "raster", &history);
  Rast_command_history(&history);
  Rast_write_history(result, &history);

  exit(EXIT_SUCCESS);
}
