#include "asf.h"
#include "ceos.h"
#include "asf_meta.h"
#include "asf_convert.h"
#include "proj.h"
#include "asf_export.h"
#include "asf_import.h"
#include "asf_contact.h"
#include "asf_sar.h"
#include "asf_terrcorr.h"
#include "asf_geocode.h"
#include "asf_nan.h"
#include "ardop_defs.h"
#include "get_ceos_names.h"
#include "get_stf_names.h"
#include "dateUtil.h"
#include <ctype.h>
#include <string.h>
#include <sys/types.h> /* 'DIR' structure (for opendir) */
#include <dirent.h>    /* for opendir itself            */

int isCEOS(const char *input_file)
{
  char **inBandName = NULL, **inMetaName = NULL;
  int nBands, trailer;

  if (require_ceos_pair(input_file, &inBandName, &inMetaName,
            &nBands, &trailer))
    return TRUE;
  else
    return FALSE;
}

int isSTF(const char *input_file)
{
  char **inBandName = NULL, **inMetaName = NULL;

  if (require_stf_pair(input_file, inBandName, inMetaName))
    return TRUE;
  else
    return FALSE;
}

static void setup_tmp_dir(char *tmp_dir)
{
  // Create temporary test directory
  strcpy(tmp_dir, "tmp_meta-");
  strcat(tmp_dir, time_stamp_dir());
  create_clean_dir(tmp_dir);
}

meta_parameters *meta_read_cfg(const char *inName, convert_config *cfg)
{
  // Setup temporary directory for metadata file
  char *tmpDir = (char *) MALLOC(sizeof(char)*512);
  setup_tmp_dir(tmpDir);

  // Get the regular metadata structure
  char **inBandName = NULL, **inMetaName = NULL;
  int nBands=1, trailer;
  if (isCEOS(inName))
    require_ceos_pair(inName, &inBandName, &inMetaName,
              &nBands, &trailer);
  else if (isSTF(inName))
    require_stf_pair(inName, inBandName, inMetaName);
  meta_parameters *meta = meta_read(inMetaName[0]);

  // Assign temporary metadata file
  char *outMetaName = (char *) MALLOC(sizeof(char)*255);
  sprintf(outMetaName, "%s/tmp.meta", tmpDir);

  // Assign a couple values out of the configuration file
  int import_single_band = 0;
  int complex_flag = cfg->import->complex_slc;
  int nBandsOut = cfg->general->terrain_correct ? nBands+1 : nBands;
  char *bandExt = (char *) MALLOC(sizeof(char)*10);
  int ii;

  // Read radiometry
  radiometry_t radiometry;
  if (strncmp_case(cfg->import->radiometry, "AMPLITUDE_IMAGE", 15) == 0)
    radiometry = r_AMP;
  else if (strncmp_case(cfg->import->radiometry, "POWER_IMAGE", 11) == 0)
    radiometry = r_POWER;
  else if (strncmp_case(cfg->import->radiometry, "SIGMA_IMAGE", 11) == 0)
    radiometry = r_SIGMA;
  else if (strncmp_case(cfg->import->radiometry, "GAMMA_IMAGE", 11) == 0)
    radiometry = r_GAMMA;
  else if (strncmp_case(cfg->import->radiometry, "BETA_IMAGE", 10) == 0)
    radiometry = r_BETA;

  // Set image data type
  data_type_t data_type = meta->general->data_type;
  char *lutName = (char *) MALLOC(sizeof(char)*255);
  strcpy(lutName, cfg->import->lut);
  if (data_type >= COMPLEX_BYTE)
    meta->general->image_data_type = COMPLEX_IMAGE;
  else if (lutName)
    meta->general->image_data_type = LUT_IMAGE;
  else
    meta->general->image_data_type = AMPLITUDE_IMAGE;

  // Set data type and band count
  for (ii=0; ii<nBands; ii++) {
    int band = ii + 1;

    // Determine the band extension (band ID)
    if (strcmp_case(meta->general->sensor, "RSAT") == 0)
      strcpy(bandExt, "HH");
    else if (strcmp_case(meta->general->sensor, "ERS") == 0 ||
         strcmp_case(meta->general->sensor, "JERS") == 0)
      strcpy(bandExt, "VV");
    else if (strcmp_case(meta->general->sensor_name, "SAR") == 0 ||
         strcmp_case(meta->general->sensor_name, "PALSAR") == 0)
      bandExt = get_polarization(inBandName[ii]);
    else if (strcmp_case(meta->general->sensor_name, "AVNIR") == 0 ||
         strcmp_case(meta->general->sensor_name, "PRISM") == 0) {
      int band_number;
      band_number = get_alos_band_number(inBandName[ii]);
      if (band_number<9)
    sprintf(bandExt, "0%d", band_number);
      else
    sprintf(bandExt, "%d", band_number);
    }
    strcpy(meta->general->bands, "");

    if (meta->sar) {
      // Assign band names
      assign_band_names(meta, outMetaName, bandExt, band, nBands, nBandsOut,
            radiometry, complex_flag);
      if (complex_flag) {
        meta->general->data_type = COMPLEX_REAL32;
        meta->general->band_count = import_single_band ? 1 : band;
      }
      else if (data_type >= COMPLEX_BYTE) {
        meta->general->data_type = REAL32;
        meta->general->band_count = import_single_band ? 2 : band*2;
      }
      else {
        meta->general->data_type = REAL32;
        meta->general->band_count = import_single_band ? 1 : band;
      }
      meta_write(meta, outMetaName);
    }
    else if (meta->optical) {
      int band_number;
      band_number = get_alos_band_number(inBandName[ii]);
      if (band_number<9)
    sprintf(bandExt, "0%d", band_number);
      else
    sprintf(bandExt, "%d", band_number);
      //if (nBands > 1)
      //  asfPrintStatus("   Input band: %s\n", bandExt);
      if (band > 1) {
    if (strcmp_case(meta->general->sensor_name, "PRISM") != 0 &&
        (strcmp_case(meta->general->mode, "1A") != 0 ||
         strcmp_case(meta->general->mode, "1B1") != 0)) {
      meta_parameters *metaTmp=NULL;
      metaTmp = meta_read(outMetaName);
      strcat(meta->general->bands, metaTmp->general->bands);
      meta_free(metaTmp);
    }
      }
      if (strcmp(meta->general->bands, "") != 0)
    strcat(meta->general->bands, ",");
      strcat(meta->general->bands, bandExt);
    }
  }
  if (nBands == 2) {
    strcpy(meta->sar->polarization, "dual-pol");
    meta->general->image_data_type = POLARIMETRIC_IMAGE;
  }
  else if (nBands == 4) {
    strcpy(meta->sar->polarization, "quad-pol");
    meta->general->image_data_type = POLARIMETRIC_IMAGE;
  }

  // Clean up
  remove_dir(tmpDir);
  FREE(tmpDir);
  FREE(outMetaName);
  FREE(bandExt);

  return meta;
}

int findDemFile(char *fileName)
{
  int found = 0;

  // First check for the filename as-is
  if (fileExists(fileName)) {
    found = 1;
  }
  else {
    // Check for .img, .tif, .tiff, .TIF, and .TIFF extensions
    char *img_File = appendExt(fileName, ".img");
    char *tif_File = appendExt(fileName, ".tif");
    char *tiff_File = appendExt(fileName, ".tiff");
    char *TIF_File = appendExt(fileName, ".TIF");
    char *TIFF_File = appendExt(fileName, ".TIFF");

    if (fileExists(img_File)) {
      found = 1;
    }
    else if (fileExists(tif_File)) {
      found = 1;
      strcpy(fileName, tif_File);
    }
    else if (fileExists(tiff_File)) {
      found = 1;
      strcpy(fileName, tiff_File);
    }
    else if (fileExists(TIF_File)) {
      found = 1;
      strcpy(fileName, TIF_File);
    }
    else if (fileExists(TIFF_File)) {
      found = 1;
      strcpy(fileName, TIFF_File);
    }
    else {
      found = 0;
    }

    // Clean up
    FREE(img_File);
    FREE(tif_File);
    FREE(tiff_File);
    FREE(TIF_File);
    FREE(TIFF_File);
  }

  return found;
}

void check_return(int ret, char *msg)
{
  if (ret != 0)
    asfPrintError(msg);
}

void check_input(convert_config *cfg, char *processing_step, char *input)
{
  char **inBandName = NULL, **inMetaName = NULL;
  int nBands, trailer;

  if (strcmp_case(processing_step, "polarimetry") == 0) {
    if (isCEOS(input))
      require_ceos_pair(input, &inBandName, &inMetaName,
            &nBands, &trailer);
    else if (isSTF(input))
      require_stf_pair(input, inBandName, inMetaName);
    meta_parameters *meta = meta_create(inMetaName[0]);
    if (meta->sar) {
      meta_free(meta);
      // re-read meta with additional info
      meta = meta_read_cfg(inMetaName[0], cfg);
      // Pauli decomposition only works for complex quad-pol data
      if (cfg->polarimetry->pauli &&
      (meta->general->image_data_type != POLARIMETRIC_IMAGE ||
       strcmp_case(meta->sar->polarization, "QUAD-POL") != 0 ||
       meta->general->band_count < 8))
        asfPrintError("Pauli decomposition requires complex quad-pol data\n");
      // Sinclair decomposition ought to work on complex and detected
      // quad-pol data
      if (cfg->polarimetry->sinclair &&
      (meta->general->image_data_type != POLARIMETRIC_IMAGE ||
       strcmp_case(meta->sar->polarization, "QUAD-POL") != 0))
        asfPrintError("Sinclair decomposition requires quad-pol data\n");
    }
    else {
      asfPrintError("Polarimetry requires SAR data.\n");
    }
  }
}

// If a temporary directory has not been specified, create one using the time
// stamp as the name of the temporary directory
// extract the file basename from the in_name, extract the directory to put
// the tmp_dir in from the out_name
static void create_and_set_tmp_dir(char *in_name, char *out_name, char *tmp_dir)
{
  int length = strlen(in_name)+1;
  if ((strlen(out_name)+1) > length) {
    length = strlen(out_name)+1;
  }
  char *junk     = MALLOC(sizeof(char)*length);
  char *basename = MALLOC(sizeof(char)*length);
  char *out_dir  = MALLOC(sizeof(char)*length);

  split_dir_and_file(in_name, junk, basename);
  split_dir_and_file(out_name, out_dir, junk);
  FREE(junk);

  int tmp_len = strlen(tmp_dir);
  int out_len = strlen(out_dir);

  if (0==tmp_len) {
    if (0==out_len) {
      strcpy(tmp_dir, "");
    }
    else {
      strcpy(tmp_dir, out_dir);
      strcat(tmp_dir, DIR_SEPARATOR_STR);
    }
    strcat(tmp_dir, basename);
    strcat(tmp_dir, "-");
    strcat(tmp_dir, time_stamp_dir());
    create_clean_dir(tmp_dir);
  }
  else {
    DIR *dirp = opendir(tmp_dir);
    if (!dirp)
      create_clean_dir(tmp_dir);
    else
      closedir(dirp);
  }

  set_asf_tmp_dir(tmp_dir);

  FREE(basename);
  FREE(out_dir);
}

/* Tracking useful intermediate files, for exposure in the GUI */
static char *intermediates_file = NULL;
extern char *g_status_file;
static void save_intermediate(convert_config *cfg, char *tag, char *filename)
{
  // the "tag" that is passed in needs to match what is looked for
  // in the GUI: asf_convert_gui/file_list.c:move_to_completed_files_list()

  if (cfg->general->intermediates && g_status_file && strlen(g_status_file)>0)
  {
    if (!intermediates_file) {
      // the name of the intermediates file needs to match what the GUI
      // will look for in asf_convert_gui/execute.c:do_convert()
      char *status_file = STRDUP(g_status_file);
      char *ext = strstr(status_file, ".status");
      if (ext) *ext = '\0';
      intermediates_file = appendExt(status_file, ".files");
      FILE *fp = fopen(intermediates_file, "w");
      if (fp) {
        fprintf(fp, "List of useful intermediate files.\n");
        fclose(fp);
      }
      free(status_file);
    }
    FILE *fp = fopen(intermediates_file, "a");
    if (fp) {
      fprintf(fp, "%s: %s\n", tag, filename);
      fclose(fp);
    }
  }
}

/* Make a copy of the metdata file. */
static void copy_meta(char *src, char *dest)
{
  char *tmp1 = appendExt(src, ".meta");
  char *tmp2 = appendExt(dest, ".meta");
  fileCopy(tmp1, tmp2);
  free(tmp1);
  free(tmp2);
}

/* Returns TRUE if the filename has a TIFF extension */
static int has_tiff_ext(const char *f)
{
    char *ext = findExt(f);
    if (ext)
        return
            strcmp_case(ext, ".tif") == 0 ||
            strcmp_case(ext, ".tiff") == 0;
    else
        return FALSE;
}

static output_format_t get_format(convert_config *cfg)
{
  output_format_t format = JPEG;
  if (strncmp(uc(cfg->export->format), "TIFF", 4) == 0) {
    format = TIF;
  } else if (strncmp(uc(cfg->export->format), "GEOTIFF", 7) == 0) {
    format = GEOTIFF;
  } else if (strncmp(uc(cfg->export->format), "JPEG", 4) == 0) {
    format = JPEG;
  } else if (strncmp(uc(cfg->export->format), "PGM", 3) == 0) {
    format = PGM;
  } else if (strncmp(uc(cfg->export->format), "PNG", 3) == 0) {
    format = PNG;
  }
  return format;
}

static char *
convert_tiff(const char *tiff_file, char *what, convert_config *cfg,
             int save_converted)
{
    char *tiff_basename, imported[255], geocoded[255], status[255];
    tiff_basename = stripExt(tiff_file);
    sprintf(imported, "%s/imported_%s", cfg->general->tmp_dir, what);

    // want see "dem" as "DEM" in the status messages, and asf_import
    // requires an uppercase string for the image_data_type
    char *uc_what = MALLOC(strlen(what)+2);
    strcpy(uc_what, uc(what)); // uc() returns ptr to static mem => make copy

    // if user wants to save the converted file, do so, otherwise we
    // put it into the temporary directory along with everything else
    if (save_converted) {
        char *outfileDir = get_dirname(cfg->general->out_name);
        char *basename = get_basename(tiff_basename);
        // the naming scheme of the geocoded dem/mask must match what the
        // gui expects in settings_update_dem/mask.
        sprintf(geocoded, "%sgeocoded_%s_%s", outfileDir, what, basename);
        free(outfileDir);
        free(basename);
    }
    else {
        sprintf(geocoded, "%s/geocoded_%s", cfg->general->tmp_dir, what);
    }

    sprintf(status, "Importing %s...", uc_what);
    update_status(status);

    sprintf(status, "ingesting GeoTIFF %s (asf_import)\n", uc_what);
    check_return(
        asf_import(r_AMP, FALSE, FALSE, FALSE, FALSE, GENERIC_GEOTIFF, NULL,
                   NULL, what, NULL, NULL, -99, -99, 0, 0, -99, -99, 0,
                   NULL, NULL, NULL, FALSE, NULL, tiff_basename, imported),
        status);

    sprintf(status, "Geocoding %s...", uc_what);
    update_status(status);

    // Check to see if the TIFF is map-projected or not, and if not
    // then geocode it.  ALSO set the image_data_type field while
    // we're at it
    meta_parameters *mp = meta_read(imported);
    if (strncmp(uc_what, "DEM", 3) == 0) {
      mp->general->image_data_type = DEM;
    }
    else if (strncmp(uc_what, "MASK", 4) == 0) {
      mp->general->image_data_type = MASK;
    }
    meta_write(mp, imported);

    // Now we allow lat/lon pseudoprojected dems with terrain correction
    // So, we can usually get away with skipping the geocoding step for
    // the dem.  However, in the case where the user wants to save the
    // cut dem, we do geocode since in that case the cut dem will use the
    // the same projection as the final result (and the layover/shadow mask,
    // if they've also saved that)
    if ((!is_map_projected(mp) && !is_lat_lon_pseudo(mp)) ||
        (cfg->terrain_correct->save_terrcorr_dem && cfg->general->geocoding))
    {
      if (cfg->general->geocoding) {
          // use the user's projection, if we have it
          sprintf(status, "geocoding GeoTIFF %s (asf_geocode)\n", uc_what);
          check_return(
              asf_geocode_from_proj_file(
                  cfg->geocoding->projection, cfg->geocoding->force,
                  RESAMPLE_NEAREST_NEIGHBOR, 0.0, WGS84_DATUM, NAN,
                  NULL, imported, geocoded, cfg->geocoding->background),
              status);
      }
      else {
          // use UTM if no geocoding specified
          sprintf(status, "geocoding GeoTIFF %s (asf_geocode)\n", uc_what);
          check_return(
              asf_geocode_utm(RESAMPLE_NEAREST_NEIGHBOR, 0.0, WGS84_DATUM,
                              NAN, NULL, imported, geocoded, 0.0), status);
      }
    }
    else {
      // is map projected, so copy the imported file to the geocoded
      // filename for consistency
      copyImgAndMeta(imported, geocoded);
    }

    meta_free(mp);
    free(uc_what);
    free(tiff_basename);

    return STRDUP(geocoded);
}

static int geocode_airsar(convert_config *cfg, const char *projection_file,
                          int force_flag, resample_method_t resample_method,
                          double average_height, datum_type_t datum, double pixel_size,
                          char *band_id, char *inFile, char *outFile,
                          float background_val)
{
    asfPrintStatus("\nGeocoding AirSAR products...\n");
    if (cfg->airsar->c_vv)
    {
        update_status("Geocoding C-band DEM...");
        char *in_tmp = appendToBasename(inFile, "_c_dem");
        char *out_tmp = appendToBasename(outFile, "_c_dem");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding dem (asf_geocode)\n");
        free(in_tmp); free(out_tmp);

        update_status("Geocoding C-band coherence...");
        in_tmp = appendToBasename(inFile, "_c_coh");
        out_tmp = appendToBasename(outFile, "_c_coh");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding coherence (asf_geocode)\n");
        free(in_tmp); free(out_tmp);

        update_status("Geocoding C-band data...");
        in_tmp = appendToBasename(inFile, "_c_vv");
        out_tmp = appendToBasename(outFile, "_c_vv");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding C-band (asf_geocode)\n");
        free(in_tmp); free(out_tmp);
    }
    else {
        asfPrintStatus("Skipping geocoding of AirSAR C-band "
                       "interferometric data.\n");
    }

    if (cfg->airsar->l_vv)
    {
        update_status("Geocoding L-band DEM...");
        char *in_tmp = appendToBasename(inFile, "_l_dem");
        char *out_tmp = appendToBasename(outFile, "_l_dem");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding dem (asf_geocode)\n");
        free(in_tmp); free(out_tmp);

        update_status("Geocoding L-band coherence...");
        in_tmp = appendToBasename(inFile, "_l_coh");
        out_tmp = appendToBasename(outFile, "_l_coh");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding coherence (asf_geocode)\n");
        free(in_tmp); free(out_tmp);

        update_status("Geocoding L-band data...");
        in_tmp = appendToBasename(inFile, "_l_vv");
        out_tmp = appendToBasename(outFile, "_l_vv");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding C-band (asf_geocode)\n");
        free(in_tmp); free(out_tmp);
    }
    else {
        asfPrintStatus("Skipping geocoding of AirSAR L-band "
                       "interferometric data.\n");
    }

    if (cfg->airsar->p_pol) {
        update_status("Geocoding Polarimetric P-band...");
        char *in_tmp = appendToBasename(inFile, "_p");
        char *out_tmp = appendToBasename(outFile, "_p");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding P-band (asf_geocode)\n");
        free(in_tmp); free(out_tmp);
    }
    else {
        asfPrintStatus("Skipping geocoding of AirSAR Polarimetric "
                       "P-band image.\n");
    }

    if (cfg->airsar->l_pol) {
        update_status("Geocoding Polarimetric L-band...");
        char *in_tmp = appendToBasename(inFile, "_l");
        char *out_tmp = appendToBasename(outFile, "_l");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding L-band (asf_geocode)\n");
        free(in_tmp); free(out_tmp);
    }
    else {
        asfPrintStatus("Skipping geocoding of AirSAR Polarimetric "
                       "L-band image.\n");
    }

    if (cfg->airsar->c_pol) {
        update_status("Geocoding Polarimetric C-band...");
        char *in_tmp = appendToBasename(inFile, "_c");
        char *out_tmp = appendToBasename(outFile, "_c");

        asfPrintStatus("Geocoding: %s -> %s\n", in_tmp, out_tmp);
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, in_tmp, out_tmp, background_val),
              "geocoding C-band (asf_geocode)\n");
        free(in_tmp); free(out_tmp);
    }
    else {
        asfPrintStatus("Skipping geocoding of AirSAR Polarimetric "
                       "C-band image.\n");
    }

    return 0; // success
}

static int check_airsar(char *outFile, char *suffix)
{
  char *base = appendToBasename(outFile, suffix);

  char *full_img = appendExt(base, ".img");
  char *full_meta = appendExt(base, ".meta");

  int ret = fileExists(full_img) && fileExists(full_meta);

  free(full_img);
  free(full_meta);
  free(base);

  return ret;
}

static scale_t get_scale(convert_config *cfg);
static void do_export(convert_config *cfg, char *inFile, char *outFile);

int asf_convert_ext(int createflag, char *configFileName, int saveDEM)
{
  convert_config *cfg;
  char inFile[255], outFile[255], tmpFile[255];
  int is_airsar=0;

  // If requested, create a config file and exit (if the file does not exist),
  // otherwise read it
  if ( createflag==TRUE && !fileExists(configFileName) ) {
    init_convert_config(configFileName);
    return(EXIT_SUCCESS);
  }
  // Extend the configuration file if the file already exist
  else if ( createflag==TRUE && fileExists(configFileName) ) {
    cfg = read_convert_config(configFileName);
    check_return(write_convert_config(configFileName, cfg),
                 "Could not update configuration file");
    asfPrintStatus("   Initialized complete configuration file\n\n");
    FCLOSE(fLog);
    return(EXIT_SUCCESS);
  }
  else if (!fileExists(configFileName)) {
    asfPrintStatus("Couldn't find config file: %s\n", configFileName);
    FCLOSE(fLog);
    return EXIT_FAILURE;
  }
  else {
    cfg = read_convert_config(configFileName);
  }
  // Using PGM with polarimetry is not allowed -- all are color output
  // (except when using entropy/anisotropy/alpha -- that could be PGM,
  // if each is exported as a separate file)
  if (strncmp(uc(cfg->export->format), "PGM", 3) == 0 &&
      (cfg->general->polarimetry && !cfg->polarimetry->cloude_pottier_nc))
  {
    asfPrintError("Greyscale PGM output is not compatible with the color "
                  "polarimetric options.\n");
  }

  // Check for greyscale PGM output versus selection of a color option
  if (strncmp(uc(cfg->export->format), "PGM", 3) == 0 &&
      (strlen(cfg->export->rgb) > 0 ||
      cfg->export->truecolor        ||
      cfg->export->falsecolor)
     )
  {
    asfPrintWarning("Greyscale PGM output is not compatible with color options:\n"
        "(RGB, True Color or False Color) ...\n"
        "Defaulting to producing separate greyscale PGM files for each\n"
        "available band.\n");
    strcpy(cfg->export->rgb, "");
    cfg->export->truecolor = 0;
    cfg->export->falsecolor = 0;
  }

  // Check for greyscale single-band output versus selection of a color option
  if (strlen(cfg->export->band) > 0 &&
      (strlen(cfg->export->rgb) > 0 ||
          cfg->export->truecolor        ||
          cfg->export->falsecolor)
     )
  {
    asfPrintWarning("Single-band output is not compatible with color options:\n"
        "(RGB, True Color or False Color) ...\n"
        "Defaulting to producing a single greyscale file for the\n"
        "selected band (%s).\n", cfg->export->band);
    strcpy(cfg->export->rgb, "");
    cfg->export->truecolor = 0;
    cfg->export->falsecolor = 0;
  }

  // Check for and resolve input/output name clashes
  // - Only applies to importing ASF-internal format files combined
  //   with non-export processing that will not change the file extension,
  //   otherwise name clashes are already resolved by the change in file
  //   extension that occurs.
  // - When checking for same-name, do it in all caps since Windows is
  //   not case-sensitive in the file naming.
  //
  if ((!cfg->general->import ||
              (cfg->general->import && strncmp(uc(cfg->import->format), "ASF", 3) == 0))  &&
       !cfg->general->export                                                              &&
      (strlen(cfg->general->out_name) <= 0 ||
              (strlen(cfg->general->in_name) == strlen(cfg->general->out_name) &&
               strcmp(uc(cfg->general->in_name), uc(cfg->general->out_name)) == 0))       &&
       strlen(cfg->general->prefix) <= 0                                                  &&
       strlen(cfg->general->suffix) <= 0                                                   )
  {
    if (cfg->general->terrain_correct || cfg->general->geocoding) {
      sprintf(cfg->general->suffix, "%s%s",
              cfg->general->suffix,
              "_out");
    }
  }

  // Mosaicking
  if (cfg->general->mosaic) {
    // Requires configuration file with batch processing switched on.
    // GUI needs to generate temporary defaults values file and data list.
    // Input is batch configuration file that contains those two parameters.
    FILE *fList;
    char *mosaic_dir=NULL, line[255], tmpList[255];
    int batch = TRUE;

    // Check whether batch processing is required
    if (strncmp(uc(cfg->import->format), "ASF", 3) == 0)
      cfg->general->import = 0;
    if (!cfg->general->import && !cfg->general->terrain_correct)
      batch = FALSE;
    else
      mosaic_dir = MALLOC(sizeof(char)*255);

    if (batch) {
      // Create a temporary directory to collect intermediate processing
      // results
      char baseName[255], baseDir[255];
      char *p = findExt(configFileName);
      if (p)
        *p = '\0';
      split_dir_and_file(configFileName, baseDir, baseName);
      strcpy(mosaic_dir,"");
      create_and_set_tmp_dir(configFileName, configFileName, mosaic_dir);

      split_dir_and_file(configFileName, baseDir, baseName);

      // Generate temporary configuration files
      char tmp_dir[255];;
      convert_config *tmp_cfg=NULL;
      char tmpCfgName[255];
      int n_ok = 0, n_bad = 0;
      FILE *fBatch = FOPEN(cfg->general->batchFile, "r");
      sprintf(tmpList, "%s/data.lst", mosaic_dir);
      fList = FOPEN(tmpList, "w");

      strcpy(tmp_dir, cfg->general->tmp_dir);
      while (fgets(line, 255, fBatch) != NULL) {
        char batchItem[255], batchItemFile[255], batchItemDir[255];
        sscanf(line, "%s", batchItem);
        
        // Create temporary directory for processing
        char *p = findExt(batchItem);
        if (p) *p = '\0';
        split_dir_and_file(batchItem, batchItemDir, batchItemFile);
        create_and_set_tmp_dir(batchItem, batchItem, tmp_dir);
        
        // Generate temporary defaults values file
        char tmpDefaults[255];
        sprintf(tmpDefaults, "%s/tmp.defaults", tmp_dir);
        FILE *fDef = FOPEN(tmpDefaults, "w");
        fprintf(fDef, "import = %d\n", cfg->general->import);
        fprintf(fDef, "polarimetry = %d\n", cfg->general->polarimetry);
        fprintf(fDef, "terrain correction = %d\n",
                cfg->general->terrain_correct);
        fprintf(fDef, "geocoding = 0\n");
        fprintf(fDef, "export = 0\n");
        fprintf(fDef, "intermediates = %d\n", cfg->general->intermediates);
        fprintf(fDef, "quiet = 1\n");
        fprintf(fDef, "short configuration file = 1\n");
        if (cfg->general->import) {
          fprintf(fDef, "input format = %s\n", cfg->import->format);
          fprintf(fDef, "radiometry = %s\n", cfg->import->radiometry);
          fprintf(fDef, "output db = %d\n", cfg->import->output_db);
          fprintf(fDef, "multilook SLC = %d\n", cfg->import->multilook_slc);
        }
        if (cfg->general->polarimetry) {
          fprintf(fDef, "pauli = %d\n", cfg->polarimetry->pauli);
          fprintf(fDef, "sinclair = %d\n", cfg->polarimetry->sinclair);
          fprintf(fDef, "freeman durden = %d\n",
                  cfg->polarimetry->freeman_durden);
          fprintf(fDef, "cloude pottier = %d\n",
                  cfg->polarimetry->cloude_pottier);
          fprintf(fDef, "faraday correction = %d\n", cfg->polarimetry->farcorr);
        }
        if (cfg->general->terrain_correct) {
          fprintf(fDef, "pixel spacing = %lf\n", cfg->terrain_correct->pixel);
          fprintf(fDef, "digital elevation model = %s\n",
                  cfg->terrain_correct->dem);
          fprintf(fDef, "mask = %s\n", cfg->terrain_correct->mask);
          fprintf(fDef, "smooth dem holes =1\n");
          fprintf(fDef, "do radiometric = %d\n",
                  cfg->terrain_correct->do_radiometric);
          fprintf(fDef, "interpolate = 1\n");
        }
        FCLOSE(fDef);
        
        // Create temporary configuration file
        sprintf(tmpCfgName, "%s/%s.cfg", tmp_dir, batchItemFile);
        FILE *fConfig = FOPEN(tmpCfgName, "w");
        fprintf(fConfig, "asf_convert temporary configuration file\n\n");
        fprintf(fConfig, "[General]\n");
        fprintf(fConfig, "default values = %s\n", tmpDefaults);
        fprintf(fConfig, "input file = %s\n", batchItem);
        fprintf(fConfig, "output file = %s%s/%s\n", batchItemDir, mosaic_dir,
                batchItemFile);
        fprintf(fConfig, "tmp dir = %s\n", tmp_dir);
        FCLOSE(fConfig);
        
        // Extend the temporary configuration file
        tmp_cfg = read_convert_config(tmpCfgName);
        check_return(write_convert_config(tmpCfgName, tmp_cfg),
                     "Could not update configuration file");
        free_convert_config(tmp_cfg);
        
        // This is really quite a kludge-- we used to call the library
        // function here, now we shell out and run the tool directly, sort
        // of a step backwards, it seems.  Unfortunately, in order to keep
        // processing the batch even if an error occurs, we're stuck with
        // this method.  (Otherwise, we'd have to teach asfPrintError to
        // get us back here, to continue the loop.)
        asfPrintStatus("\nProcessing %s ...\n", batchItem);
        char cmd[1024];
        if (logflag) {
          sprintf(cmd, "%sasf_convert%s -log %s %s",
              get_argv0(), bin_postfix(), logFile, tmpCfgName);
        }
        else {
          sprintf(cmd, "%sasf_convert%s %s", get_argv0(), bin_postfix(),
                  tmpCfgName);
        }
        int ret = asfSystem(cmd);
        
        if (ret != 0) {
          asfPrintStatus("%s: failed\n", batchItem);
          ++n_bad;
        }
        else {
          asfPrintStatus("%s: ok\n", batchItem);
          fprintf(fList, "%s/%s.img\n", mosaic_dir, batchItem);
          ++n_ok;
        }
        
        strcpy(tmp_dir, cfg->general->tmp_dir);
      }
      FCLOSE(fBatch);
      FCLOSE(fList);
    }

    // Read file names to pass to mosaicking
    int ii, nFiles = 0;
    if (mosaic_dir)
      fList = FOPEN(tmpList, "r");
    else
      fList = FOPEN(cfg->general->batchFile, "r");
    while (fgets(line, 255, fList) != NULL)
      nFiles++;
    FCLOSE(fList);
    char *in_base_names[nFiles+1];
    for (ii=0; ii<nFiles; ii++)
      in_base_names[ii] = MALLOC(sizeof(char)*255);
    in_base_names[nFiles] = NULL;
    ii = 0;
    if (mosaic_dir)
      fList = FOPEN(tmpList, "r");
    else
      fList = FOPEN(cfg->general->batchFile, "r");
    while (fgets(line, 255, fList) != NULL) {
      line[strlen(line)-1] = '\0';
      sprintf(in_base_names[ii], line);
      ++ii;
    }
    FCLOSE(fList);

    // Mosaic imported (and terrain corrected) images
    double lat_min = -999, lon_min = -999;
    double lat_max = 999, lon_max = 999;
    project_parameters_t pp;
    projection_type_t proj_type;
    datum_type_t datum;
    resample_method_t resample_method;
    int multiband = 1;
    int band_num = 0;
    char *err=NULL;

    if (strcmp(cfg->geocoding->resampling, "NEAREST NEIGHBOR") == 0)
      resample_method = RESAMPLE_NEAREST_NEIGHBOR;
    else if (strcmp(cfg->geocoding->resampling, "BILINEAR") == 0)
      resample_method = RESAMPLE_BILINEAR;
    else if(strcmp(cfg->geocoding->resampling, "BICUBIC") == 0)
      resample_method = RESAMPLE_BICUBIC;

    if (mosaic_dir)
      sprintf(outFile, "%s/%s", mosaic_dir, cfg->general->out_name);
    else
      sprintf(outFile, "%s", cfg->general->out_name);
    if (!parse_proj_args_file(cfg->geocoding->projection, &pp, &proj_type, &datum, &err)) {
      asfPrintError("%s",err);
    }
    update_status("Mosaicking...");
    asf_mosaic(&pp, proj_type, cfg->geocoding->force, resample_method,
           cfg->geocoding->height, datum, cfg->geocoding->pixel,
           multiband, band_num, in_base_names, outFile,
           cfg->geocoding->background,
           lat_min, lat_max, lon_min, lon_max, cfg->mosaic->overlap);

    // Export mosaic
    if (cfg->general->export) {
      sprintf(inFile, "%s", outFile);
      sprintf(outFile, "%s", cfg->general->out_name);
      update_status("Exporting...");
      do_export(cfg, inFile, outFile);
    }
  }

  // Batch mode processing
  else if (strlen(cfg->general->batchFile) > 0) {
    convert_config *tmp_cfg=NULL;
    char tmp_dir[255];
    char tmpCfgName[255];
    char line[255];
    int n_ok = 0, n_bad = 0;
    FILE *fBatch = FOPEN(cfg->general->batchFile, "r");

    strcpy(tmp_dir, cfg->general->tmp_dir);
    while (fgets(line, 255, fBatch) != NULL) {
      char batchItem[255], batchItemFile[255], batchItemDir[255];
      sscanf(line, "%s", batchItem);

      // strip off known extensions
      char *p = findExt(batchItem);
      if (p) *p = '\0';

      split_dir_and_file(batchItem, batchItemDir, batchItemFile);

      char *tmpDir = MALLOC(sizeof(char)*(strlen(cfg->general->defaults)+1));
      char *tmpFile = MALLOC(sizeof(char)*(strlen(cfg->general->defaults)+1));
      split_dir_and_file(cfg->general->defaults, tmpDir, tmpFile);
      char cwd[10000];
      getcwd(cwd,10000);
      int defaultsLen;
      if (strlen(cfg->general->defaults) > strlen(cwd)+strlen(tmpFile)) {
        defaultsLen = strlen(cfg->general->defaults);
      } else {
        defaultsLen = strlen(cwd)+strlen(tmpFile);
      }
      char *defaults = MALLOC(sizeof(char)*defaultsLen+2);
      if (0==strlen(tmpDir)) {
        sprintf(defaults,"%s/%s",cwd,tmpFile);
      } else {
        strcpy(defaults,cfg->general->defaults);
      }
      FREE(tmpDir);
      FREE(tmpFile);

      // Create temporary configuration file
      create_and_set_tmp_dir(batchItem, cfg->general->default_out_dir, tmp_dir);
      sprintf(tmpCfgName, "%s/%s.cfg", tmp_dir, batchItemFile);


      FILE *fConfig = FOPEN(tmpCfgName, "w");
      fprintf(fConfig, "asf_mapready temporary configuration file\n\n");
      fprintf(fConfig, "[General]\n");
      fprintf(fConfig, "default values = %s\n", defaults);
      fprintf(fConfig, "input file = %s\n", batchItem);
      if (strlen(cfg->general->default_out_dir) == 0)
        fprintf(fConfig, "output file = %s%s%s\n",
                cfg->general->prefix, batchItemFile, cfg->general->suffix);
      else
        fprintf(fConfig, "output file = %s%c%s%s%s\n",
                cfg->general->default_out_dir, DIR_SEPARATOR,
                cfg->general->prefix, batchItemFile, cfg->general->suffix);
      fprintf(fConfig, "tmp dir = %s\n", tmp_dir);
      FCLOSE(fConfig);
      FREE(defaults);

      // Extend the temporary configuration file
      tmp_cfg = read_convert_config(tmpCfgName);
      check_return(write_convert_config(tmpCfgName, tmp_cfg),
                   "Could not update configuration file");
      free_convert_config(tmp_cfg);

      // Run asf_mapready for temporary configuration file

      // This is really quite a kludge-- we used to call the library
      // function here, now we shell out and run the tool directly, sort
      // of a step backwards, it seems.  Unfortunately, in order to keep
      // processing the batch even if an error occurs, we're stuck with
      // this method.  (Otherwise, we'd have to teach asfPrintError to
      // get us back here, to continue the loop.)
      asfPrintStatus("\nProcessing %s ...\n", batchItem);
      char cmd[1024];
      if (logflag) {
          sprintf(cmd, "%sasf_mapready%s -log %s %s",
                get_argv0(), bin_postfix(), logFile, tmpCfgName);
      }
      else {
          sprintf(cmd, "%sasf_mapready%s %s",
                get_argv0(), bin_postfix(), tmpCfgName);
      }
      int ret = asfSystem(cmd);

      if (ret != 0) {
          asfPrintStatus("%s: failed\n", batchItem);
          ++n_bad;
      } else {
          asfPrintStatus("%s: ok\n", batchItem);
          ++n_ok;
      }

      strcpy(tmp_dir, cfg->general->tmp_dir);
    }
    FCLOSE(fBatch);

    asfPrintStatus("\n\nBatch Complete.\n");
    asfPrintStatus("Successfully processed %d/%d file%s.\n", n_ok,
        n_ok + n_bad, n_ok + n_bad == 1 ? "" : "s");

    if (n_bad > 0)
        asfPrintStatus("  *** %d file%s failed. ***\n", n_bad,
            n_bad==1 ? "" : "s");
  }
  // Regular processing
  else {

    if (cfg->general->status_file && strlen(cfg->general->status_file) > 0)
      set_status_file(cfg->general->status_file);

    update_status("Processing...");

    // these are so we can tell how long processing took
    ymd_date start_date;
    hms_time start_time;
    get_current_date(&start_date, &start_time);

    char tmp[64];
    date_printTime(&start_time,0,':',tmp);
    asfPrintStatus("Starting at: %s\n", tmp);

    create_and_set_tmp_dir(cfg->general->in_name, cfg->general->out_name,
                           cfg->general->tmp_dir);
    save_intermediate(cfg, "Temp Dir", cfg->general->tmp_dir);

    // Check that input name isn't the same as the output name
    // (This can happen with CEOS Level 0-- both use .raw)
    if (strcmp_case(cfg->import->format, "CEOS") == 0 &&
        strcmp_case(cfg->general->in_name, cfg->general->out_name) == 0) {
        // the names (actually, the basenames) can match for L1 data
        // it is only the ".RAW" ones we want to flag
        char *RAW_file = appendExt(cfg->general->in_name, ".RAW");
        char *raw_file = appendExt(cfg->general->in_name, ".raw");
        if (fileExists(raw_file) || fileExists(RAW_file))
            asfPrintError("Input and Output names are the same: %s\n",
                          cfg->general->in_name);
        free(raw_file); free(RAW_file);
        RAW_file = appendStr(cfg->general->in_name, ".RAW");
        raw_file = appendStr(cfg->general->in_name, ".raw");
        if (fileExists(raw_file) || fileExists(RAW_file))
            asfPrintError("Input and Output names are the same: %s\n",
                          cfg->general->in_name);
        free(raw_file); free(RAW_file);
    }

    // Check if the user said the file is L1, but we really have a L0 file
    if (strcmp_case(cfg->import->format, "CEOS (1)") == 0) {
        char *raw_file = appendExt(cfg->general->in_name, ".raw");
        char *ldr_file = appendExt(cfg->general->in_name, ".ldr");
        if (fileExists(raw_file) && fileExists(ldr_file))
          asfPrintError("You selected CEOS Level 1 import, however the input "
                        "file appears to be Level 0:\n  %s\n", raw_file);
        free(raw_file); free(ldr_file);

        char *RAW_file = appendExt(cfg->general->in_name, ".RAW");
        char *LDR_file = appendExt(cfg->general->in_name, ".LDR");
        if (fileExists(RAW_file) && fileExists(LDR_file))
          asfPrintError("You selected CEOS Level 1 import, however the input "
                        "file appears to be Level 0:\n  %s\n", RAW_file);
        free(RAW_file); free(LDR_file);
    }

    // Check whether everything in the [Import] block is reasonable
    if (cfg->general->import) {

      // Import format
      if (strncmp(uc(cfg->import->format), "ASF", 3) != 0 &&
          strncmp(uc(cfg->import->format), "CEOS", 4) != 0 &&
          strncmp(uc(cfg->import->format), "STF", 3) != 0 &&
          strncmp(uc(cfg->import->format), "AIRSAR", 6) != 0 &&
          strncmp(uc(cfg->import->format), "BIL", 3) != 0 &&
          strncmp(uc(cfg->import->format), "GRIDFLOAT", 9) != 0 &&
          strncmp(uc(cfg->import->format), "GAMMA", 5) != 0 &&
          strncmp(uc(cfg->import->format), "GEOTIFF", 7) != 0) {
        asfPrintError("Chosen import format not supported\n");
      }

      is_airsar = strncmp(uc(cfg->import->format), "AIRSAR", 6) == 0;

      // Radiometry
      if (strncmp(uc(cfg->import->radiometry), "AMPLITUDE_IMAGE", 15) != 0 &&
          strncmp(uc(cfg->import->radiometry), "POWER_IMAGE", 11) != 0 &&
          strncmp(uc(cfg->import->radiometry), "SIGMA_IMAGE", 11) != 0 &&
          strncmp(uc(cfg->import->radiometry), "GAMMA_IMAGE", 11) != 0 &&
          strncmp(uc(cfg->import->radiometry), "BETA_IMAGE", 10) != 0) {
        asfPrintError("Chosen radiometry not supported\n");
      }

      // Look up table file existence check
      if (strlen(cfg->import->lut) > 0) {
        if (!fileExists(cfg->import->lut)) {
          asfPrintError("Look up table file does not exist\n");
        }
      }

      // When importing AirSAR data, don't allow terrain correction
      if (is_airsar && cfg->general->terrain_correct) {
        asfPrintError("Terrain correction of AirSAR data is not supported.\n");
      }

      // When importing a GeoTIFF, don't allow terrain correction.
      if (strncmp(uc(cfg->import->format), "GEOTIFF", 7) == 0 && 
          cfg->general->terrain_correct)
      {
        asfPrintError("Terrain correction of GeoTIFFs is not supported.\n");
      }

      // Precision state vector file check can only be done
      // from within asf_import


      // Dealing with single look complex data
      // Options -multilook and -complex are mutually exclusive
      if (cfg->import->complex_slc && cfg->import->multilook_slc)
        asfPrintError("Only single look complex data stored as "
                      "amplitude and phase can be multilooked.\n");

      // Get input file name ready
      strcpy(inFile, cfg->general->in_name);

      // Can skip import if the input is already asf internal.
      if (strncmp(uc(cfg->import->format), "ASF", 3) == 0) {
          cfg->general->import = 0;
      }
    }

    // Check whether everything in the [SAR processing] block is reasonable
    if (cfg->general->sar_processing) {

        // ADDED FOR 3.0 -- DO NOT SUPPORT L0 PROCESSING!
        // We actually think this is working in many cases, but it hasn't
        // been fully tested yet for all satellites/beams.
        asfPrintWarning("Processing from level 0 is not fully tested yet!\n");

      // Radiometry
      if (strncmp(uc(cfg->sar_processing->radiometry), "AMPLITUDE_IMAGE", 15) != 0 &&
          strncmp(uc(cfg->sar_processing->radiometry), "POWER_IMAGE", 11) != 0 &&
          strncmp(uc(cfg->sar_processing->radiometry), "SIGMA_IMAGE", 11) != 0 &&
          strncmp(uc(cfg->sar_processing->radiometry), "GAMMA_IMAGE", 11) != 0 &&
          strncmp(uc(cfg->sar_processing->radiometry), "BETA_IMAGE", 10) != 0) {
        asfPrintError("Chosen radiometry not supported\n");
      }

    }

    // Check whether everything in the [C2P] block is reasonable
    if (cfg->general->c2p) {
        asfPrintWarning("SLC ingest with c2p is still underdevelopment.\n");
        // make sure input data has the extension .cpx
        char *ext = findExt(cfg->general->in_name);
        if (strcmp_case(ext, ".img") != 0) {
            asfPrintWarning("Input data is not complex.  c2p flag ignored.\n");
            cfg->general->c2p = 0;
        }
    }

    // Check whether everything in the [Image stats] block is reasonable
    if (cfg->general->image_stats) {

      // Values
      if (strncmp(cfg->image_stats->values, "LOOK", 4) != 0 &&
          strncmp(cfg->image_stats->values, "INCIDENCE", 9) != 0 &&
          strncmp(cfg->image_stats->values, "RANGE", 5) != 0) {
        asfPrintError("Chosen values not supported\n");
      }
    }

    // Check whether everything in the [Polarimetry] block is reasonable
    if (cfg->general->polarimetry) {
      int pauli = cfg->polarimetry->pauli == 0 ? 0 : 1;
      int sinclair = cfg->polarimetry->sinclair == 0 ? 0 : 1;
      int cloude_pottier = cfg->polarimetry->cloude_pottier == 0 ? 0 : 1;
      int cloude_pottier_ext =
          cfg->polarimetry->cloude_pottier_ext == 0 ? 0 : 1;
      int cloude_pottier_nc = cfg->polarimetry->cloude_pottier_nc == 0 ? 0 : 1;
      int k_means_wishart = cfg->polarimetry->k_means_wishart == 0 ? 0 : 1;
      int freeman_durden = cfg->polarimetry->freeman_durden == 0 ? 0 : 1;
      int k_means_wishart_ext =
          cfg->polarimetry->k_means_wishart_ext == 0 ? 0 : 1;
      if (pauli + sinclair + cloude_pottier + cloude_pottier_ext +
          cloude_pottier_nc + k_means_wishart + k_means_wishart_ext +
          freeman_durden > 1)
        asfPrintError("More than one polarimetric processing scheme selected."
                      "\nOnly one of these options may be selected at a time."
                      "\n");
      if (pauli + sinclair + cloude_pottier + cloude_pottier_ext +
          cloude_pottier_nc + k_means_wishart + k_means_wishart_ext +
          freeman_durden + cfg->polarimetry->farcorr < 1)
        asfPrintError("No polarimetric processing scheme selected.\n");
    }

    // Check whether everything in the [Terrain correction] block is
    // reasonable
    if (cfg->general->terrain_correct) {
      // Reference DEM file existence check
      if (!findDemFile(cfg->terrain_correct->dem)) {
        asfPrintError("Reference DEM file '%s' does not exist\n",
                      cfg->terrain_correct->dem);
      }

      // Check for pixel size smaller than threshold ???

      // specified a mask and asked for an auto-mask
      int have_mask = cfg->terrain_correct->mask &&
          strlen(cfg->terrain_correct->mask) > 0;
      if (have_mask) {
          if (!findDemFile(cfg->terrain_correct->mask)) {
              asfPrintError("MASK file '%s' does not exist\n",
                            cfg->terrain_correct->mask);
          }

          if (cfg->terrain_correct->auto_mask_water) {
              asfPrintStatus("Mask File: %s\n", cfg->terrain_correct->mask);
              asfPrintStatus("Auto Water Mask: %d\n",
                             cfg->terrain_correct->auto_mask_water);
              asfPrintWarning(
                  "You cannot specify a mask for terrain correction "
                  "and also request an automatically generated mask.\n"
                  "Ignoring auto water mask -- will use the mask provided.\n");
          }
      }
    }

    // Check whether everything in the [Geocoding] block is reasonable
    if (cfg->general->geocoding) {

      // Projection file existence check
      if (!fileExists(cfg->geocoding->projection)) {
        asfPrintError("Projection parameter file '%s' does not exist\n",
                      cfg->geocoding->projection);
      }

      // Check for pixel size smaller than threshold ???

      // Datum
      if (meta_is_valid_string(cfg->geocoding->datum)          &&
          strlen(cfg->geocoding->datum) > 0                    &&
          strncmp(uc(cfg->geocoding->datum), "NONE", 4)   != 0 &&
          strncmp(uc(cfg->geocoding->datum), "WGS84", 5)  != 0 &&
          strncmp(uc(cfg->geocoding->datum), "NAD27", 5)  != 0 &&
          strncmp(uc(cfg->geocoding->datum), "NAD83", 5)  != 0 &&
          strncmp(uc(cfg->geocoding->datum), "HUGHES", 6) != 0) {
        asfPrintError("Chosen datum not supported\n");
      }

      // Resampling
      if (strncmp(uc(cfg->geocoding->resampling), "NEAREST_NEIGHBOR", 16) != 0 &&
          strncmp(uc(cfg->geocoding->resampling), "BILINEAR", 8) != 0 &&
          strncmp(uc(cfg->geocoding->resampling), "BICUBIC", 7) != 0) {
        asfPrintError("Chosen resampling method not supported\n");
      }

      // Check that the user didn't specify an average height, and
      // also is doing terrain correction
      if (cfg->general->terrain_correct &&
          !cfg->terrain_correct->refine_geolocation_only &&
          cfg->geocoding->height != 0)
      {
          asfPrintWarning("Since terrain correction is being applied, "
                          "asf_geocode will ignore the\nspecified "
                          "average height.\n");
      }
    }

    // Check whether everything in the [Terrain Correct] block is reasonable
    if (cfg->general->terrain_correct) {

        // specified a mask and asked for an auto-mask
        int have_mask = cfg->terrain_correct->mask &&
            strlen(cfg->terrain_correct->mask) > 0;

        if (have_mask && cfg->terrain_correct->auto_mask_water) {
            asfPrintStatus("Mask File: %s\n", cfg->terrain_correct->mask);
            asfPrintStatus("Auto Water Mask: %d\n",
                           cfg->terrain_correct->auto_mask_water);
            asfPrintWarning("You cannot specify a mask for terrain correction "
                          "and also request an automatically generated mask.\n"
                          "Ignoring auto water mask -- will use the mask "
                          "provided.\n");
        }

    }

    // Check whether everything in the [Export] block is reasonable
    if (cfg->general->export) {

      // Export format: ASF, TIFF, GEOTIFF, JPEG, PGM
      if (strncmp(uc(cfg->export->format), "ASF", 3) != 0 &&
          strncmp(uc(cfg->export->format), "TIFF", 4) != 0 &&
          strncmp(uc(cfg->export->format), "GEOTIFF", 7) != 0 &&
          strncmp(uc(cfg->export->format), "JPEG", 4) != 0 &&
          strncmp(uc(cfg->export->format), "PNG", 3) != 0 &&
          strncmp(uc(cfg->export->format), "PGM", 3) != 0) {
        asfPrintError("Chosen export format (%s) not supported\n",
                     cfg->export->format);
      }

      // Unset export flag when export format is ASF
      if (strncmp(uc(cfg->export->format), "ASF", 3) == 0) {
        cfg->general->export = 0;
      }

      // Scaling method for floating point GeoTIFFs
      if (strncmp(uc(cfg->export->format), "GEOTIFF", 7) == 0) {
        if (strlen(cfg->export->byte) < 1) {
          strcpy(cfg->export->byte, "NONE");
        }
        if (strncmp(uc(cfg->export->byte), "NONE", 4) != 0 &&
            strncmp(uc(cfg->export->byte), "SIGMA", 5) != 0 &&
            strncmp(uc(cfg->export->byte), "MINMAX", 6) != 0 &&
            strncmp(uc(cfg->export->byte), "TRUNCATE", 8) != 0 &&
            strncmp(uc(cfg->export->byte), "HISTOGRAM_EQUALIZE", 18) != 0) {
          asfPrintError("Chosen scaling method (%s) not supported\n",
                        cfg->export->byte);
        }
      }

      // If RGB Banding option is "ignore,ignore,ignore" then the
      // user has probably been using the gui, and didn't pick
      // anything for any of the RGB channels.
      int rgbBanding = strlen(cfg->export->rgb) > 0 ? 1 : 0;
      if (rgbBanding &&
          strcmp(uc(cfg->export->rgb), "IGNORE,IGNORE,IGNORE") == 0)
      {
          asfPrintError(
              "Exporting as RGB was selected, but no values for each RGB\n"
              "channels were selected.  Please choose which bands you wish\n"
              "to have placed into at least one of the RGB channels.\n");
      }

      // RGB Banding, True Color, and False Color are mutually-exclusive
      // options.  Check to make sure the config doesn't spec more than
      // one of these.
    }
    int rgbBanding = strlen(cfg->export->rgb) > 0 ? 1 : 0;
    int truecolor = cfg->export->truecolor == 0 ? 0 : 1;
    int falsecolor = cfg->export->falsecolor == 0 ? 0 : 1;
    if (rgbBanding + truecolor + falsecolor > 1)
    {
      asfPrintError("More than one color export mode was selected (rgb banding, truecolor\n"
          "or falsecolor)  Only one of these options may be selected at a time.\n");
    }

    if (!cfg->general->import && !cfg->general->sar_processing &&
        !cfg->general->polarimetry && !cfg->general->terrain_correct &&
    !cfg->general->geocoding && !cfg->general->export) {
      asfPrintError("Invalid configuration file found (%s) or no processing\n"
                    "is enabled, e.g. import, terrain correction, geocoding,\n"
                    "and export processing flags are all set to zero in the\n"
                    "configuration file.\n", configFileName);
    }

    // Check input
    if (cfg->general->polarimetry)
      check_input(cfg, "polarimetry", cfg->general->in_name);

    //---------------------------------------------------------------
    // Let's finally get to work
    if (strlen(cfg->general->out_name) == 0) {
      sprintf(cfg->general->out_name, "%s", cfg->general->in_name);
    }
    sprintf(outFile, "%s", cfg->general->out_name);

    // global variable-- if set, tells meta_write to also dump .hdr (ENVI) files
    dump_envi_header = cfg->general->dump_envi;

    // we add a "secret" band 0 containing amplitude data, if we are
    // generating something other than amplitude, and are going to
    // be terrain correcting
    int amp0_flag = FALSE;

    // The "saved radiometry" is the requested radiometry in the case where
    // import will not apply the calibration parameters, because the farcorr
    // code will be doing it -- in that case, import must use amplitude,
    // and "saved_radiometry" will be passed in to the FR correction func.
    radiometry_t saved_radiometry = r_AMP;

    // Call asf_import, if needed (=> input is not ASF Internal)
    if (cfg->general->import) {

      update_status("Importing...");

      radiometry_t radiometry;
      int db_flag = FALSE;
      int lut_flag = FALSE;
      input_format_t format_type;

      // Radiometry
      if (!cfg->general->sar_processing) {
          if (strncmp(uc(cfg->import->radiometry), "AMPLITUDE_IMAGE", 15) == 0) {
              radiometry = r_AMP;
          } else if (strncmp(uc(cfg->import->radiometry), "POWER_IMAGE", 11) == 0) {
              radiometry = r_POWER;
          } else if (strncmp(uc(cfg->import->radiometry), "SIGMA_IMAGE", 11) == 0) {
              radiometry = r_SIGMA;
          } else if (strncmp(uc(cfg->import->radiometry), "GAMMA_IMAGE", 11) == 0) {
              radiometry = r_GAMMA;
          } else if (strncmp(uc(cfg->import->radiometry), "BETA_IMAGE", 10) == 0) {
              radiometry = r_BETA;
          } else {
        radiometry = r_AMP;  // Default to something (silence compiler)
          }
      }

      // See above, this is the flag that adds a "secret" AMP band
      // for terrain correction
      amp0_flag = radiometry != r_AMP && cfg->general->terrain_correct;
      if (amp0_flag)
        asfPrintStatus("Adding Amplitude band, for terrain correction.\n");

      // LUT
      if (strlen(cfg->import->lut) > 0)
          lut_flag = TRUE;

      if (cfg->import->output_db)
          db_flag = TRUE;

      // When performing Faraday Rotation, we must import amplitude.
      // After the correction, the correct radiometry will be applied.
      // We include the db_flag in the saved_radiometry.
      if (cfg->general->polarimetry && cfg->polarimetry &&
          cfg->polarimetry->farcorr != FARCORR_OFF)
      {
        saved_radiometry = radiometry;
        if (db_flag) {
          if (saved_radiometry == r_SIGMA) saved_radiometry = r_SIGMA_DB;
          if (saved_radiometry == r_BETA)  saved_radiometry = r_BETA_DB;
          if (saved_radiometry == r_GAMMA) saved_radiometry = r_GAMMA_DB;
        }

        // now set radiometry to amplitude, and db off.
        radiometry = r_AMP;
        db_flag = FALSE;
      }

      // Generate a temporary output filename
      if (cfg->general->image_stats || cfg->general->detect_cr ||
          cfg->general->sar_processing || cfg->general->polarimetry ||
          cfg->general->terrain_correct || cfg->general->geocoding ||
          cfg->general->export) {
        sprintf(outFile, "%s/import", cfg->general->tmp_dir);
      }
      else {
        sprintf(outFile, "%s", cfg->general->out_name);
      }

      // Input Format Type
      if (strncmp_case(cfg->import->format, "CEOS", 4) == 0)
        format_type = CEOS;
      else if (strncmp_case(cfg->import->format, "STF", 3) == 0)
        format_type = STF;
      else if (strncmp_case(cfg->import->format, "GEOTIFF", 7) == 0)
        format_type = GENERIC_GEOTIFF;
      else if (strncmp_case(cfg->import->format, "BIL", 3) == 0)
        format_type = BIL;
      else if (strncmp_case(cfg->import->format, "GRIDFLOAT", 9) == 0)
        format_type = GRIDFLOAT;
      else if (strncmp_case(cfg->import->format, "AIRSAR", 6) == 0)
        format_type = AIRSAR;
      else if (strncmp_case(cfg->import->format, "GAMMA_MSP", 9) == 0)
        format_type = GAMMA_MSP;
      else if (strncmp_case(cfg->import->format, "GAMMA_ISP", 9) == 0)
        format_type = GAMMA_ISP;
      else if (strncmp_case(cfg->import->format, "VP", 2) == 0)
        format_type = VP;
      else {
        asfPrintError("Unknown Format: %s\n", cfg->import->format);
        format_type = CEOS; // actually this is not reached
      }

      // Call asf_import!
      check_return(asf_import(radiometry, db_flag,
                              cfg->import->complex_slc,
                              cfg->import->multilook_slc,
                              amp0_flag,
                              format_type,
                              NULL,
                              NULL,
                              MAGIC_UNSET_STRING,
                              cfg->import->lut, cfg->import->prc,
                              cfg->import->lat_begin, cfg->import->lat_end,
                              cfg->import->line, cfg->import->sample,
                              cfg->import->width, cfg->import->height,
                              cfg->general->intermediates, NULL, NULL, NULL,
                              cfg->import->ers2_gain_fix,
                              NULL,
                              cfg->general->in_name, outFile),
                              "ingesting data file (asf_import)\n");

      // For AirSAR data, let's see what we actually got.
      // Not all products are always present - update settings to
      // disable processing of products not present
      if (is_airsar) {
        if (cfg->airsar->l_pol && !check_airsar(outFile, "_l")) {
          asfPrintStatus("No L-band polarimetric AirSAR product.\n");
          cfg->airsar->l_pol = 0;
        }
        if (cfg->airsar->p_pol && !check_airsar(outFile, "_p")) {
          asfPrintStatus("No P-band polarimetric AirSAR product.\n");
          cfg->airsar->p_pol = 0;
        }
        if (cfg->airsar->c_pol && !check_airsar(outFile, "_c")) {
          asfPrintStatus("No C-band polarimetric AirSAR product.\n");
          cfg->airsar->c_pol = 0;
        }
        if (cfg->airsar->c_vv && !check_airsar(outFile, "_c_vv")) {
          asfPrintStatus("No C-band interferometric AirSAR product.\n");
          cfg->airsar->c_vv = 0;
        }
        if (cfg->airsar->l_vv && !check_airsar(outFile, "_l_vv")) {
          asfPrintStatus("No L-band interferometric AirSAR product.\n");
          cfg->airsar->l_vv = 0;
        }

        if (!cfg->airsar->l_pol && !cfg->airsar->p_pol &&
            !cfg->airsar->c_pol && !cfg->airsar->c_vv &&
            !cfg->airsar->l_vv)
          asfPrintError("No airsar products to process!\n");
      }

      if (!is_airsar) {
        // Make sure truecolor/falsecolor are only specified for optical data
        meta_parameters *meta = meta_read(outFile);
        if (meta->optical && cfg->general->terrain_correct) {
          asfPrintError("Terrain correction cannot be applied to optical images (...and\n"
              "orthorectification is not yet supported.)\n");
        }
        if (!meta->optical && (truecolor || falsecolor)) {
          asfPrintError("Cannot select True Color or False Color output "
                        "with non-optical data\n");
        }
        if (cfg->export->band &&
            strlen(cfg->export->band) > 0 &&
            get_band_number(meta->general->bands,
                            meta->general->band_count, cfg->export->band) < 0)
        {
          asfPrintError("Selected export band (%s) does not exist: \n"
                        "   Imported file: %s.img\n"
                        " Available bands: %s\n",
                        cfg->export->band, cfg->general->in_name,
                        meta->general->bands);
        }
        meta_free(meta);
      }
    }
    else {
      // skipping import ==> "outFile" (what import should have produced)
      // is really just the original input file (already ASF Internal)
      strcpy(outFile, cfg->general->in_name);
    }

    if (cfg->general->sar_processing) {
      if (is_airsar)
        asfPrintError("Cannot perform SAR Processing on AirSAR data.\n");

      update_status("Running ArDop...");

      // Check whether the input file is a raw image.
      // If not, skip the SAR processing step
      meta_parameters *meta = meta_read(outFile);
      if (meta->general->image_data_type == RAW_IMAGE) {

          sprintf(inFile, "%s", outFile);
          if (cfg->general->image_stats || cfg->general->detect_cr ||
              cfg->general->polarimetry || cfg->general->terrain_correct ||
              cfg->general->geocoding || cfg->general->export) {
            sprintf(outFile, "%s/sar_processing", cfg->general->tmp_dir);
          }
          else {
            sprintf(outFile, "%s", cfg->general->out_name);
          }

          struct INPUT_ARDOP_PARAMS *params_in;
          params_in = get_input_ardop_params_struct(inFile, outFile);
          int one = 1;

        // Radiometry
        if (strncmp(uc(cfg->sar_processing->radiometry), "POWER_IMAGE", 11) == 0) {
            params_in->pwrFlag = &one;
        } else if (strncmp(uc(cfg->sar_processing->radiometry), "SIGMA_IMAGE", 11) == 0) {
            params_in->sigmaFlag = &one;
        } else if (strncmp(uc(cfg->sar_processing->radiometry), "GAMMA_IMAGE", 11) == 0) {
            params_in->gammaFlag = &one;
        } else if (strncmp(uc(cfg->sar_processing->radiometry), "BETA_IMAGE", 10) == 0) {
            params_in->betaFlag = &one;
        }

        // When terrain correcting, add the deskew flag
        if (cfg->general->terrain_correct)
            params_in->deskew = &one;

        ardop(params_in);

        if (strcmp_case(cfg->sar_processing->radiometry, "AMPLITUDE_IMAGE") == 0)
          sprintf(outFile, "%s/sar_processing_amp", cfg->general->tmp_dir);
        else if (strcmp_case(cfg->sar_processing->radiometry, "POWER_IMAGE") == 0)
          sprintf(outFile, "%s/sar_processing_power", cfg->general->tmp_dir);
        else if (strcmp_case(cfg->sar_processing->radiometry, "SIGMA_IMAGE") == 0)
          sprintf(outFile, "%s/sar_processing_sigma", cfg->general->tmp_dir);
        else if (strcmp_case(cfg->sar_processing->radiometry, "GAMMA_IMAGE") == 0)
          sprintf(outFile, "%s/sar_processing_gamma", cfg->general->tmp_dir);
        else if (strcmp_case(cfg->sar_processing->radiometry, "BETA_IMAGE") == 0)
          sprintf(outFile, "%s/sar_processing_beta", cfg->general->tmp_dir);
        else
            asfPrintError("Unexpected radiometry: %s\n",
                          cfg->sar_processing->radiometry);

        // If we are not going to be terrain correcting, get the output
        // from ardop into ground range (terrain correction would do that
        // for us).
        if (!cfg->general->terrain_correct)
        {
            asfPrintStatus("Converting to ground range.\n");
            sprintf(inFile, "%s", outFile);
            sprintf(outFile, "%s_gr", inFile);
            sr2gr(inFile, outFile);
        }
      }
      else {
        asfPrintStatus("Image has already been processed - skipping SAR processing step");
      }
      meta_free(meta);
    }

    if (cfg->general->c2p) {
      char inDataName[512];
      char *baseName = (char *) MALLOC(sizeof(char)*512);
      baseName = get_basename(inFile);
      sprintf(inDataName, "%s.img", baseName);

        update_status("Converting Complex to Polar...");

        sprintf(inFile, "%s", outFile);
        if (cfg->general->polarimetry || cfg->general->terrain_correct ||
            cfg->general->geocoding || cfg->general->export)
        {
            sprintf(outFile, "%s/c2p", cfg->general->tmp_dir);
        }
        else
        {
            sprintf(outFile, "%s", cfg->general->out_name);
        }

        c2p(inDataName, outFile, cfg->c2p->multilook, TRUE);
    }

    if (cfg->general->image_stats) {
      char values[255];

      update_status("Running Image Stats...");

      // Values for statistics
      if (strncmp(uc(cfg->image_stats->values), "LOOK", 4) == 0) {
        sprintf(values, "-look");
      } else if (strncmp(uc(cfg->image_stats->values), "INCIDENCE", 9) == 0) {
        sprintf(values, "-incidence");
      } else if (strncmp(uc(cfg->image_stats->values), "RANGE", 5) == 0) {
        sprintf(values, "-range");
      }

      // Pass in command line
      sprintf(inFile, "%s", outFile);
      if (cfg->general->polarimetry || cfg->general->terrain_correct ||
          cfg->general->geocoding || cfg->general->export)
      {
        sprintf(outFile, "%s/image_stats", cfg->general->tmp_dir);
      }
      else {
        sprintf(outFile, "%s", cfg->general->out_name);
      }
      check_return(image_stats(inFile, outFile, values, cfg->image_stats->bins,
                               cfg->image_stats->interval),
                   "running statistics on data file (image_stats)\n");
    }

    if (cfg->general->detect_cr) {

      update_status("Detecting Corner Reflectors...");

      // Intermediate results
      if (cfg->general->intermediates) {
        cfg->detect_cr->chips = 1;
        cfg->detect_cr->text = 1;
      }

      // Pass in command line
      sprintf(inFile, "%s", outFile);
      if (cfg->general->polarimetry || cfg->general->polarimetry ||
          cfg->general->geocoding || cfg->general->export) {
        sprintf(outFile, "%s/detect_cr", cfg->general->tmp_dir);
      }
      else {
        sprintf(outFile, "%s", cfg->general->out_name);
      }
      check_return(detect_cr(inFile, cfg->detect_cr->cr_location, outFile,
                             cfg->detect_cr->chips, cfg->detect_cr->text),
                   "detecting corner reflectors (detect_cr)\n");
    }


    if (cfg->general->polarimetry) {

      int doing_pol = cfg->polarimetry->pauli
                      + cfg->polarimetry->sinclair
                      + cfg->polarimetry->cloude_pottier
                      + cfg->polarimetry->cloude_pottier_ext 
                      + cfg->polarimetry->cloude_pottier_nc
                      + cfg->polarimetry->freeman_durden
                      + cfg->polarimetry->k_means_wishart
                      + cfg->polarimetry->k_means_wishart_ext;
      int doing_far = cfg->polarimetry->farcorr;

      amp0_flag = cfg->general->terrain_correct;

      if (doing_far) {
        update_status("Applying Faraday rotation correction ...");

        // Pass in command line for faraday correction
        sprintf(inFile, "%s", outFile);
        if (cfg->general->terrain_correct || cfg->general->geocoding ||
            cfg->general->export)
        {
          sprintf(outFile, "%s/faraday_correction", cfg->general->tmp_dir);
        }
        else {
          sprintf(outFile, "%s", cfg->general->out_name);
        }
      
        int keep_flag = cfg->general->intermediates;
        int single_angle_flag = (FARCORR_MEAN == cfg->polarimetry->farcorr);
        radiometry_t rad = saved_radiometry;
        asfPrintStatus("\nApplying Faraday Rotation correction.\n");
        faraday_correct(inFile, outFile, keep_flag, single_angle_flag, rad);
        asfPrintStatus("Done.\n\n");

        sprintf(tmpFile, "%s/polarimetry_farcorr.img", cfg->general->tmp_dir);
        save_intermediate(cfg, "Faraday", tmpFile);
      }

      if (doing_pol) {
        update_status("Polarimetric processing ...");

        // Pass in command line for polarimetry
        sprintf(inFile, "%s", outFile);
        if (cfg->general->terrain_correct ||
            cfg->general->geocoding ||
            cfg->general->export)
        {
          sprintf(outFile, "%s/polarimetry", cfg->general->tmp_dir);
        }
        else {
          sprintf(outFile, "%s", cfg->general->out_name);
        }
  
        // Calculate polarimetric parameters
        if (cfg->polarimetry->pauli)
          cpx2pauli(inFile, outFile, cfg->general->terrain_correct);
        else if (cfg->polarimetry->sinclair) {
          // for sinclair, there are two possibilities: SLC & non-SLC
          meta_parameters *meta = meta_read(inFile);
          if (meta->general->band_count >= 8) {
            cpx2sinclair(inFile, outFile, cfg->general->terrain_correct);
          }
          else {
            // here, we don't need to do any processing -- we just need to
            // update the RGB Bands to Red=HH, Green=HV, Blue=VV
            strcpy(cfg->export->rgb, "SIGMA-HH,SIGMA-HV,SIGMA-VV");
            strcpy(outFile, inFile);
  
            // turn off the amp0_flag -- we don't need it in this case
            amp0_flag = FALSE;
          }
          meta_free(meta);
        }
        else if (cfg->polarimetry->cloude_pottier) {
          cpx2cloude_pottier8(inFile, outFile, cfg->general->terrain_correct);
        }
        else if (cfg->polarimetry->cloude_pottier_ext)
          cpx2cloude_pottier16(inFile, outFile, cfg->general->terrain_correct);
        else if (cfg->polarimetry->cloude_pottier_nc)
          cpx2entropy_anisotropy_alpha(inFile, outFile,
                                       cfg->general->terrain_correct);
        else if (cfg->polarimetry->freeman_durden)
          cpx2freeman_durden(inFile, outFile, cfg->general->terrain_correct);
        else if (cfg->polarimetry->k_means_wishart)
          asfPrintError("K-means Wishart clustering not supported yet.\n");
        else if (cfg->polarimetry->k_means_wishart_ext)
          asfPrintError("Extended K-means Wishart clustering not supported yet.\n");
        else if (cfg->polarimetry->lee_preserving)
          asfPrintError("Lee category preserving not supported yet.\n");
        else
          asfPrintError("Unsupported polarimetric processing technique.\n");

        if (cfg->polarimetry->cloude_pottier ||
            cfg->polarimetry->cloude_pottier_ext ||
            cfg->polarimetry->cloude_pottier_nc)
        {
          sprintf(tmpFile, "%s/polarimetry_ea_hist.img",
                  cfg->general->tmp_dir);
          save_intermediate(cfg, "Entropy-Alpha Histogram", tmpFile);          
          sprintf(tmpFile, "%s/polarimetry_class_map.img",
                  cfg->general->tmp_dir);
          save_intermediate(cfg, "Entropy-Alpha Class Map", tmpFile);          
        }
      }
    }

    if (cfg->general->terrain_correct) {

      // If the DEM is a GeoTIFF, we need to import it, and geocode it.
      if (has_tiff_ext(cfg->terrain_correct->dem)) {
          char * converted_dem =
              convert_tiff(cfg->terrain_correct->dem, "dem", cfg, saveDEM);
          strcpy(cfg->terrain_correct->dem, converted_dem);
          free(converted_dem);
      }

      // If the Mask is a GeoTIFF, we need to import it, and geocode it.
      if (cfg->terrain_correct->mask &&
          strlen(cfg->terrain_correct->mask) > 0 &&
          has_tiff_ext(cfg->terrain_correct->mask))
      {
          char * converted_mask =
              convert_tiff(cfg->terrain_correct->mask, "mask", cfg, saveDEM);
          strcpy(cfg->terrain_correct->mask, converted_mask);
          free(converted_mask);
      }

      // Generate filenames
      sprintf(inFile, "%s", outFile);
      sprintf(outFile, "%s/terrain_correct", cfg->general->tmp_dir);

      // Call asf_terrcorr!  Or refine_geolocation!
      if (cfg->terrain_correct->refine_geolocation_only) {
          update_status("Refining Geolocation...");
          check_return(
              refine_geolocation(inFile, cfg->terrain_correct->dem,
                                 cfg->terrain_correct->mask, outFile, FALSE,
                                 cfg->terrain_correct->auto_mask_water,
                                 cfg->terrain_correct->water_height_cutoff,
                                 !cfg->general->intermediates, NULL),
            "refining geolocation of the data file (refine_geolocation)\n");
      }
      else {
          update_status("Terrain Correcting...");
          check_return(
            asf_terrcorr_ext(inFile, cfg->terrain_correct->dem,
                             cfg->terrain_correct->mask, outFile,
                             cfg->terrain_correct->pixel,
                             !cfg->general->intermediates,
                             !cfg->terrain_correct->no_resampling,
                             FALSE, cfg->terrain_correct->interp,
                             TRUE, 20, TRUE, cfg->terrain_correct->fill_value,
                             cfg->terrain_correct->auto_mask_water,
                             cfg->terrain_correct->save_terrcorr_dem, FALSE,
                             cfg->terrain_correct->water_height_cutoff,
                             cfg->terrain_correct->do_radiometric,
                             cfg->terrain_correct->smooth_dem_holes,
                             NULL, cfg->terrain_correct->no_matching,
                 cfg->terrain_correct->range_offset,
                 cfg->terrain_correct->azimuth_offset),
            "terrain correcting data file (asf_terrcorr)\n");
      }

      // save the simulated sar image intermediate
      char *dem_basename = get_basename(cfg->terrain_correct->dem);
      sprintf(tmpFile, "%s/%s_sim_sar.img", cfg->general->tmp_dir,
              dem_basename);
      save_intermediate(cfg, "Simulated SAR", tmpFile);
      sprintf(tmpFile, "%s/import_slant.img", cfg->general->tmp_dir);
      save_intermediate(cfg, "Imported Slant Range", tmpFile);
      free(dem_basename);


      // If we added a "secret" AMP band to the beginning of the
      // file, we can remove it now
      if (amp0_flag) {
        asfPrintStatus("\nRemoving added amplitude band...\n");
        remove_band(outFile, 0, cfg->general->intermediates);
      }

      if (!cfg->general->export && !cfg->general->geocoding) {
          // if this was the last step, get the terrain corrected output
          // to the output directory -- as well as any other needed files
          renameImgAndMeta(outFile, cfg->general->out_name);

          // this is to get the thumbnail code all set -- it will use
          // "outFile" as the file to generate the thumbnail from
          strcpy(outFile, cfg->general->out_name);
      }
    }

    datum_type_t datum = WGS84_DATUM;

    if (cfg->general->geocoding) {
      update_status("Geocoding...");
      int force_flag = cfg->geocoding->force;
      resample_method_t resample_method = RESAMPLE_BILINEAR;
      double average_height = cfg->geocoding->height;
      double pixel_size = cfg->geocoding->pixel;
      float background_val = cfg->geocoding->background;

      // When terrain correcting, ignore average height -- the height
      // has already been corrected for.
      if (cfg->general->terrain_correct &&
          !cfg->terrain_correct->refine_geolocation_only &&
          cfg->geocoding->height != 0)
      {
          asfPrintWarning("Since terrain correction was applied, ignoring "
                          "average height specification\nfor geocoding.\n");
          average_height = 0.0;
      }

      // Datum
      if (strncmp(uc(cfg->geocoding->datum), "WGS84", 5) == 0) {
        datum = WGS84_DATUM;
      }
      else if (strncmp(uc(cfg->geocoding->datum), "NAD27", 5) == 0) {
        datum = NAD27_DATUM;
      }
      else if (strncmp(uc(cfg->geocoding->datum), "NAD83", 5) == 0) {
        datum = NAD83_DATUM;
      }
      else if (strncmp(uc(cfg->geocoding->datum), "HUGHES", 6) == 0 ) {
        datum = HUGHES_DATUM;
      }
      else {
        datum = WGS84_DATUM;
        asfPrintWarning("Unrecognized, missing, or unsupported datum found in configuration\n"
            "file (%s).  Defaulting to WGS-84 unless reading parameters from a\n"
            "projection definition (*.proj) file.  The proj file settings will\n"
            "override all else.\n", uc(cfg->geocoding->datum));
      }

      // Resampling method
      if (strncmp(uc(cfg->geocoding->resampling), "NEAREST_NEIGHBOR", 16) == 0) {
        resample_method = RESAMPLE_NEAREST_NEIGHBOR;
      }
      if (strncmp(uc(cfg->geocoding->resampling), "BILINEAR", 8) == 0) {
        resample_method = RESAMPLE_BILINEAR;
      }
      if (strncmp(uc(cfg->geocoding->resampling), "BICUBIC", 7) == 0) {
        resample_method = RESAMPLE_BICUBIC;
      }

      sprintf(inFile, "%s", outFile);
      if (cfg->general->export) {
        sprintf(outFile, "%s/geocoding", cfg->general->tmp_dir);
      }
      else {
        sprintf(outFile, "%s%s%s",
                cfg->general->prefix,
                cfg->general->out_name,
                cfg->general->suffix);
      }

      // Pass in command line
      if (is_airsar) {
        // airsar -- geocode only what was asked for
        check_return(geocode_airsar(cfg, cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum, pixel_size,
                NULL, inFile, outFile, background_val),
              "geocoding airsar (asf_geocode)\n");

      } else {
        // normal geocoding
        check_return(asf_geocode_from_proj_file(cfg->geocoding->projection,
                force_flag, resample_method, average_height, datum,
                pixel_size, NULL, inFile, outFile, background_val),
              "geocoding data file (asf_geocode)\n");
      }
    }

    if (cfg->general->export) {

      // Set up filenames
      sprintf(inFile, "%s", outFile);
      sprintf(outFile, "%s", cfg->general->out_name);

      if (!is_airsar)
      {
        // Do normal export
        update_status("Exporting... ");
        do_export(cfg, inFile, outFile);
      }
      else
      {
        // AirSAR export -- must export a bunch of different stuff

        // do multi-band stuff first
        asfPrintStatus("\nExporting AirSAR products...\n");
        if (cfg->airsar->p_pol) {
          update_status("Exporting polarimetric P-band...");
          char *in_tmp = appendToBasename(inFile, "_p");
          char *out_tmp = appendToBasename(outFile, "_p");

          asfPrintStatus("Exporting P-band: %s -> %s\n", in_tmp, out_tmp);
          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);
        }
        else {
          asfPrintStatus("Skipping export of AirSAR P-band data.\n");
        }

        if (cfg->airsar->l_pol) {
          update_status("Exporting polarimetric L-band...");
          char *in_tmp = appendToBasename(inFile, "_l");
          char *out_tmp = appendToBasename(outFile, "_l");

          asfPrintStatus("Exporting L-band: %s -> %s\n", in_tmp, out_tmp);
          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);
        }
        else {
          asfPrintStatus("Skipping export of AirSAR L-band data.\n");
        }

        if (cfg->airsar->c_pol) {
          update_status("Exporting polarimetric C-band...");
          char *in_tmp = appendToBasename(inFile, "_c");
          char *out_tmp = appendToBasename(outFile, "_c");

          asfPrintStatus("Exporting C-band: %s -> %s\n", in_tmp, out_tmp);
          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);
        }
        else {
          asfPrintStatus("Skipping export of AirSAR C-band data.\n");
        }

        // those were two multi-band images -- the rest are single.
        // so they must be export as greyscale, no matter what the user
        // actually asked for... temporarily reset the convert_config
        char *rgb = STRDUP(cfg->export->rgb);
        strcpy(cfg->export->rgb, "");

        if (cfg->airsar->c_vv) {
          update_status("Exporting C-band...");
          char *in_tmp = appendToBasename(inFile, "_c_vv");
          char *out_tmp = appendToBasename(outFile, "_c_vv");

          asfPrintStatus("Exporting C-band: %s -> %s\n", in_tmp, out_tmp);
          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);

          update_status("Exporting C-band DEM...");
          in_tmp = appendToBasename(inFile, "_c_dem");
          out_tmp = appendToBasename(outFile, "_c_dem");

          asfPrintStatus("Exporting C-band DEM: %s -> %s\n", in_tmp, out_tmp);

          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);

          update_status("Exporting C-band coherence...");
          in_tmp = appendToBasename(inFile, "_c_coh");
          out_tmp = appendToBasename(outFile, "_c_coh");

          asfPrintStatus("Exporting C-band coherence: %s -> %s\n",
                         in_tmp, out_tmp);
          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);
        }
        else {
          asfPrintStatus("Skipping export of AirSAR C-band interferometric "
                         "data.\n");
        }

        if (cfg->airsar->l_vv) {
          update_status("Exporting L-band...");
          char *in_tmp = appendToBasename(inFile, "_l_vv");
          char *out_tmp = appendToBasename(outFile, "_l_vv");

          asfPrintStatus("Exporting L-band: %s -> %s\n", in_tmp, out_tmp);
          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);

          update_status("Exporting L-band DEM...");
          in_tmp = appendToBasename(inFile, "_l_dem");
          out_tmp = appendToBasename(outFile, "_l_dem");

          asfPrintStatus("Exporting L-band DEM: %s -> %s\n", in_tmp, out_tmp);

          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);

          update_status("Exporting L-band coherence...");
          in_tmp = appendToBasename(inFile, "_l_coh");
          out_tmp = appendToBasename(outFile, "_l_coh");

          asfPrintStatus("Exporting L-band coherence: %s -> %s\n",
                         in_tmp, out_tmp);
          do_export(cfg, in_tmp, out_tmp);
          free(in_tmp); free(out_tmp);
        }
        else {
          asfPrintStatus("Skipping export of AirSAR L-band interferometric "
                         "data.\n");
        }

        // now put back the user's rgb settings
        strcpy(cfg->export->rgb, rgb);
        free(rgb);

        // airsar metadata...
        if (cfg->general->import) {
          // export c/l_vv band's metadata as the "official" metadata
          char *in_tmp=NULL;
          if (cfg->airsar->c_vv) {
            in_tmp = appendToBasename(inFile, "_c_vv");
          } else if (cfg->airsar->l_vv) {
            in_tmp = appendToBasename(inFile, "_l_vv");
          } else {
            // .meta should already be there, even if user did not request
            // it, because asf_import imports everything anyway
            // ==> however, it will be called "import_c_vv"
            char *dir = get_dirname(inFile);
            in_tmp = MALLOC(sizeof(char)*(strlen(dir)+32));

            if (strlen(dir) > 0)
              sprintf(in_tmp, "%simport_c_vv.meta", dir);
            else
              strcpy(in_tmp, "import_c_vv.meta");

            if (!fileExists(in_tmp)) {
              // try "import_l_vv", data may not have contained c-band
              // interferometric data...
              if (strlen(dir) > 0)
                sprintf(in_tmp, "%simport_l_vv.meta", dir);
              else
                strcpy(in_tmp, "import_l_vv.meta");

              if (!fileExists(in_tmp)) {
                // this is bad - we can't produce an output metadata file!
                asfPrintWarning("Failed to generate an output metadata file!\n");
                free(in_tmp);
                in_tmp = NULL;
              }
            }
          }
          if (in_tmp) {
            copy_meta(in_tmp, outFile);
            free(in_tmp);
          }
        }
      }
    }

    //---------------------------------------------------------------------
    // At this point the processing of the SAR image is done.
    // We'll now do some of the extra stuff the user may have asked for.

    // Generate a small thumbnail if requested.
    if (cfg->general->thumbnail) {
        asfPrintStatus("Generating Thumbnail image...\n");

        output_format_t format = PNG;

        // scale using the same scaling as we did for the output image,
        // unless it was a floating-point output (not possible for PNG),
        // in which case we will use SIGMA.
        scale_t scale = get_scale(cfg);
        if (scale == NONE) scale = SIGMA;

        meta_parameters *meta;
        double in_side_length, out_x_pixel_size, out_y_pixel_size;
        char *tmpFile;
        int i,n;

        tmpFile = (char *) MALLOC(sizeof(char)*512);

        if (!cfg->general->export)
            sprintf(inFile, "%s", outFile);

        if (is_airsar) {
          int found=FALSE;
          if (cfg->airsar->c_vv) {
            char *tmp = appendToBasename(inFile, "_c_vv");
            char *tmp_img = appendExt(tmp, ".img");
            if (fileExists(tmp_img)) {
              found=TRUE;
              strcpy(inFile, tmp);
            }
            free(tmp); free(tmp_img);
          }
          if (!found && cfg->airsar->l_vv) {
            char *tmp = appendToBasename(inFile, "_l_vv");
            char *tmp_img = appendExt(tmp, ".img");
            if (fileExists(tmp_img)) {
              found=TRUE;
              strcpy(inFile, tmp);
            }
            free(tmp); free(tmp_img);
          }
          if (!found && cfg->airsar->c_pol) {
            char *tmp = appendToBasename(inFile, "_c");
            char *tmp_img = appendExt(tmp, ".img");
            if (fileExists(tmp_img)) {
              found=TRUE;
              strcpy(inFile, tmp);
            }
            free(tmp); free(tmp_img);
          }
          if (!found && cfg->airsar->l_pol) {
            char *tmp = appendToBasename(inFile, "_l");
            char *tmp_img = appendExt(tmp, ".img");
            if (fileExists(tmp_img)) {
              found=TRUE;
              strcpy(inFile, tmp);
            }
            free(tmp); free(tmp_img);
          }
          if (!found && cfg->airsar->p_pol) {
            char *tmp = appendToBasename(inFile, "_p");
            char *tmp_img = appendExt(tmp, ".img");
            if (fileExists(tmp_img)) {
              found=TRUE;
              strcpy(inFile, tmp);
            }
            free(tmp); free(tmp_img);
          }
          if (!found) {
            strcpy(inFile, "");
            asfPrintWarning("No thumbnail available.\n");
          }
        }

        if (strlen(inFile) > 0) {
          // Calculate pixel size for generating right size thumbnail
          meta = meta_read(inFile);

          // Put the thumbnail in the intermediates directory, if it is
          // being kept, otherwise in the output directory.
          char *basename = get_basename(cfg->general->out_name);

          if (cfg->general->intermediates) {
            sprintf(outFile, "%s/%s_thumb.png",
                    cfg->general->tmp_dir, basename);
            sprintf(tmpFile, "%s/%s_thumb",
                    cfg->general->tmp_dir, basename);
          } else {
            char *tmp = appendToBasename(cfg->general->out_name, "_thumb");
            strcpy(tmpFile, tmp);
            strcpy(outFile, tmp);
            strcat(outFile, ".png");
            free(tmp);
          }

          if (strlen(cfg->export->lut) > 0) {
            // using a LUT -- we can't downsample in this case, because
            // that will screw up the behavior of TRUNCATE.  The burden
            // will be on the GUI to do the downsampling, but after we have
            // generated an RGB image.
            char *bands[2];
            bands[0] = meta->general->bands;
            bands[1] = NULL;

            check_return(
              asf_export_bands(format, TRUNCATE, TRUE, 0, 0, 0, 0,
                               cfg->export->lut, inFile, outFile, bands),
              "exporting thumbnail (asf_export), using rgb look up table.\n");
          }
          else {
            // non-LUT case
            in_side_length =
              (meta->general->line_count > meta->general->sample_count) ?
              meta->general->line_count : meta->general->sample_count;

            // calculate the reduced pixel sizes
            out_x_pixel_size = meta->general->x_pixel_size*in_side_length/512.;
            out_y_pixel_size = meta->general->y_pixel_size*in_side_length/512.;

            check_return(
              resample_to_pixsiz(inFile, tmpFile, out_x_pixel_size,
                                 out_y_pixel_size),
              "resampling data to thumbnail size (resample)\n");

            if (strlen(cfg->export->rgb) > 0) {
              char *red, *green, *blue;
              if (split3(cfg->export->rgb, &red, &green, &blue, ',')) {
                char **bands = find_bands(inFile, 1, red, green, blue, &n);
                if (n > 0) {
                  check_return(
                    asf_export_bands(format, scale, TRUE, 0, 0, 0, 0, NULL,
                                     tmpFile, outFile, bands),
                    "exporting thumbnail data file (asf_export), banded\n");
                  for (i=0; i<n; ++i)
                    FREE(bands[i]);
                  FREE(bands);
                }
              }
            }
            else { // no rgb bands selected
              int true_color = cfg->export->truecolor == 0 ? 0 : 1;
              int false_color = cfg->export->falsecolor == 0 ? 0 : 1;
              if (meta->optical && (true_color || false_color))
              {
                if (meta->optical && (true_color || false_color)) {
                  // Multi-band optical data, exporting as true or
                  // false color single file
                  char **bands = extract_band_names(meta->general->bands,
                                                    meta->general->band_count);
                  if (meta->general->band_count >= 4 && bands != NULL) {
                    // The imported file IS a multiband file with enough bands,
                    // but the extract bands need to be ordered correctly
                    if (true_color) {
                      strcpy(bands[0], "03");
                      strcpy(bands[1], "02");
                      strcpy(bands[2], "01");
                      strcpy(bands[3], "");
                    }
                    else {
                      strcpy(bands[0], "04");
                      strcpy(bands[1], "03");
                      strcpy(bands[2], "02");
                      strcpy(bands[3], "");
                    }
                    check_return(
                      asf_export_bands(format, NONE, TRUE,
                                       true_color, false_color, 0, 0, NULL,
                                       tmpFile, outFile, bands),
                      "exporting thumbnail (asf_export), color banded.\n");
                    for (i=0; i<meta->general->band_count; ++i)
                      FREE (bands[i]);
                    FREE(bands);
                  }
                }
              }
              else { // not a true or false color optical image
                if (meta->general->band_count == 1) {
                  check_return(asf_export(format, scale, tmpFile, outFile),
                               "exporting thumbnail data file (asf_export)\n");
                }
                else {
                  // just the first band -- set the others to NULL
                  char **bands = find_single_band(tmpFile, "all", &n);
                  for (i=1; i<n; ++i) {
                    FREE(bands[i]); bands[i] = NULL;
                  }
                  check_return(
                    asf_export_bands(format, scale, FALSE, 0, 0,
                                     0, 0, NULL, tmpFile, outFile, bands),
                    "exporting thumbnail data file (asf_export)\n");
                  // strip off the band name at the end!
                  char *banded_name =
                    MALLOC(sizeof(char)*(strlen(outFile)+32));
                  if (cfg->general->intermediates) {
                    sprintf(banded_name, "%s/%s_thumb_%s.png",
                            cfg->general->tmp_dir, basename, bands[0]);
                    sprintf(outFile, "%s/%s_thumb.png",
                            cfg->general->tmp_dir, basename);
                  }
                  else {
                    sprintf(banded_name, "%s_thumb_%s.png",
                            cfg->general->out_name, bands[0]);
                    char *tmp = appendToBasename(cfg->general->out_name,
                                                 "_thumb");
                    strcpy(outFile, tmp);
                    strcat(outFile, ".png");
                    free(tmp);
                  }
                  fileRename(banded_name, outFile);
                  FREE(bands[0]);
                  FREE(bands);
                  FREE(banded_name);
                }
              }
            }
          }

          meta_free(meta);
          free(basename);
        }
    }

    // Process the clipped DEM if requested
    if (cfg->terrain_correct->save_terrcorr_dem) {
        // We know the name of the cut DEM in the temporary directory
        char *tmp;
        if (cfg->terrain_correct->smooth_dem_holes)
          tmp = appendToBasename(cfg->terrain_correct->dem, "_tc_smooth_cut");
        else
          tmp = appendToBasename(cfg->terrain_correct->dem, "_cut");
        char *tmp2 = get_basename(tmp);
        sprintf(inFile, "%s/%s", cfg->general->tmp_dir, tmp2);
        free(tmp);
        free(tmp2);

        // output name will be the SAR image's name with a "_dem" added
        // to the basename
        tmp = appendToBasename(cfg->general->out_name, "_dem");
        strcpy(outFile, tmp);
        free(tmp);

        //Never re-geocode the DEM -- assume user has already put it into
        //their favored projection (since at this time we require that
        //DEMs be geocoded for terrain correction ingest).  So, proceed
        //directly to export.
        if (cfg->general->export) {
            update_status("Exporting clipped DEM...");

            check_return(
                asf_export(get_format(cfg), SIGMA, inFile, outFile),
                "exporting clipped dem (asf_export)\n");
        }
        else {
            // User requested that we save the clipped DEM, but chose not
            // to export.  So... just move the clipped DEM out of the tmp dir
            renameImgAndMeta(inFile, outFile);
        }

        save_intermediate(cfg, "Clipped DEM", outFile);
    }

    // Process the layover/shadow mask if requested
    if (cfg->terrain_correct->save_terrcorr_layover_mask) {
        if (cfg->general->geocoding) {
            update_status("Geocoding layover mask...");
            sprintf(inFile, "%s/terrain_correct_mask",cfg->general->tmp_dir);
            sprintf(outFile, "%s/layover_mask_geocoded",cfg->general->tmp_dir);
            check_return(
                asf_geocode_from_proj_file(
                    cfg->geocoding->projection, cfg->geocoding->force,
                    RESAMPLE_NEAREST_NEIGHBOR, cfg->geocoding->height,
                    datum, cfg->geocoding->pixel, NULL, inFile, outFile,
                    MASK_INVALID_DATA),
                "geocoding clipped DEM (asf_geocode)\n");
        }
        else {
            // no geocoding ... just prepare the 'outFile' param for export
            sprintf(outFile, "%s/terrain_correct_mask", cfg->general->tmp_dir);
        }

        sprintf(inFile, "%s", outFile);
        char *tmp = appendToBasename(cfg->general->out_name, "_layover_mask");
        strcpy(outFile, tmp);
        free(tmp);

        if (cfg->general->export) {
            update_status("Exporting layover mask...");

            meta_parameters *meta = meta_read(inFile);

            char *bands[2];
            bands[0] = meta->general->bands;
            bands[1] = NULL;

            check_return(
                asf_export_bands(get_format(cfg), TRUNCATE, 1, 0, 0, 0, 0,
                                 "layover_mask.lut", inFile, outFile, bands),
                "exporting layover mask (asf_export)\n");

            meta_free(meta);
        }
        else {
            // no export... just move the geocoded file out of the
            // temporary directory
            renameImgAndMeta(inFile, outFile);
        }

        save_intermediate(cfg, "Layover/Shadow Mask", outFile);
    }

    if (!cfg->general->intermediates) {
        remove_dir(cfg->general->tmp_dir);
    }

    // figure out how long the processing took, tell the user all about it
    ymd_date end_date;
    hms_time end_time;
    get_current_date(&end_date, &end_time);
    double elapsed = date_difference(&end_date, &end_time,
                                     &start_date, &start_time);

    date_printTime(&end_time,0,':',tmp);
    asfPrintStatus("Done at: %s\n", tmp);

    if (elapsed < 60) {
      asfPrintStatus("Elapsed time: %d seconds.\n", (int)elapsed);
    }
    else if (elapsed < 60*60) {
      int min = elapsed/60;
      int sec = elapsed-min*60;
      asfPrintStatus("Elapsed time: %d minute%s, %d second%s.\n",
                     min, min==1?"":"s", sec, sec==1?"":"s");
    }
    else if (elapsed < 60*60*24) {
      int hour = elapsed/(60*60);
      double rem = elapsed-hour*60*60;
      int min = rem/60;
      int sec = rem-min*60;
      asfPrintStatus("Elapsed time: %d hour%s, %d minute%s, %d second%s.\n",
                     hour, hour==1?"":"s", min, min==1?"":"s",
                     sec, sec==1?"":"s");
    }
    else {
      int day = elapsed/(60*60*24);
      double rem = elapsed-day*60*60*24;
      int hour = rem/(60*60);
      rem -= hour*60*60;
      int min = rem/60;
      int sec = rem-min*60;
      asfPrintStatus("Elapsed time: %d day%s, %d hour%s, %d minute%s, "
                     "%d second%s.\n",
                     day, day==1?"":"s", hour, hour==1?"":"s",
                     min, min==1?"":"s", sec, sec==1?"":"s");
    }

    // Don't change this message unless you also change the code in
    // asf_convert_gui/execute.c to look for a different successful
    // completion string.  GUI currently detects successful processing
    // by looking for this message in the log file.... (yeah, I know.)
    asfPrintStatus("\nSuccessful completion!\n\n");
  }

  update_status("Done");
  free_convert_config(cfg);
  clear_status_file();

  return(EXIT_SUCCESS);
}

int asf_convert(int createflag, char *configFileName)
{
    return asf_convert_ext(createflag, configFileName, FALSE);
}

static scale_t get_scale(convert_config *cfg)
{
  scale_t scale = SIGMA;

  // Byte scaling
  if (strncmp(uc(cfg->export->byte), "TRUNCATE", 8) == 0) {
    scale = TRUNCATE;
  } else if (strncmp(uc(cfg->export->byte), "MINMAX", 6) == 0) {
    scale = MINMAX;
  } else if (strncmp(uc(cfg->export->byte), "SIGMA", 5) == 0) {
    scale = SIGMA;
  } else if (strncmp(uc(cfg->export->byte), "HISTOGRAM_EQUALIZE", 18) == 0) {
    scale = HISTOGRAM_EQUALIZE;
  } else if (strncmp(uc(cfg->export->byte), "NONE", 4) == 0) {
    scale = NONE;
  }

  return scale;
}

static void do_export(convert_config *cfg, char *inFile, char *outFile)
{
  int true_color = cfg->export->truecolor == 0 ? 0 : 1;
  int false_color = cfg->export->falsecolor == 0 ? 0 : 1;
  output_format_t format = get_format(cfg);
  scale_t scale = get_scale(cfg);

  // Move the .meta file out of temporary status
  // Don't need to do this if we skipped import, we'd already have .meta
  int i,is_airsar = strncmp_case(cfg->import->format, "AIRSAR", 6)==0;
  if (cfg->general->import && !is_airsar)
    copy_meta(inFile, outFile);

  meta_parameters *meta = meta_read(inFile);
  if (strlen(cfg->export->lut) > 0) {
    if (meta->general->band_count != 1) {
      asfPrintWarning("RGB Look-up-table not allowed for multi-band images."
                      " Ignored.\n");
      strcpy(cfg->export->lut,"");
    }
    else {
      if (strlen(cfg->export->rgb) > 0) {
        asfPrintWarning("RGB Banding option not allowed with RGB look up table."
                        " Ignored.\n");
      }
      if (scale != TRUNCATE) {
        asfPrintWarning("Scale option %s not allowed with RGB look-up-table."
                        " Using TRUNCATE\n",
                        cfg->export->byte);
      }
      asfPrintStatus("Exporting using look-up-table: %s\n", cfg->export->lut);

      char *bands[2];
      bands[0] = meta->general->bands;
      bands[1] = NULL;

      check_return(
        asf_export_bands(format, TRUNCATE, TRUE, 0, 0, 0, 0, cfg->export->lut,
                         inFile, outFile, bands),
        "exporting data file (asf_export), using rgb look up table.\n");
    }
  }

  if (strlen(cfg->export->lut)==0) {
    // non look-up-table case
    if (strlen(cfg->export->rgb) > 0) {
      // user has requested banding
      char *red,  *green, *blue;
      if (split3(cfg->export->rgb, &red, &green, &blue, ',')) {
        int num_found;
        char **bands = find_bands(inFile, TRUE, red, green, blue,
                                  &num_found);
        if (num_found > 0) {
          asfPrintStatus("\nExporting RGB bands into single color file:\n"
                         "Red band  : %s\n"
                         "Green band: %s\n"
                         "Blue band : %s\n\n",
                         red, green, blue);
          check_return(asf_export_bands(format, scale, TRUE, 0, 0, 0, 0, NULL,
                                        inFile, outFile, bands),
                       "export data file (asf_export), banded.\n");
          for (i=0; i<num_found; i++)
            FREE(bands[i]);
          FREE(bands);
        } else {
          asfPrintError("The requested bands (%s) "
                        "were not available in the file: %s\n",
                        cfg->export->rgb, inFile);
          strcpy(cfg->export->rgb, "");
        }
        FREE(red); FREE(green); FREE(blue);
      } else {
        asfPrintError("Invalid listing of RGB bands: %s\n",
                      cfg->export->rgb);
        strcpy(cfg->export->rgb, "");
      }
    }

    if (strlen(cfg->export->rgb) == 0)
    {
      if (meta->optical && (true_color || false_color)) {
        // Multi-band optical data, exporting as true or false color single file
        asfPrintStatus("\nExporting %s file...\n\n\n",
                       true_color ? "True Color" : false_color ? "False Color" : "Unknown");
        if (true_color &&
            (strstr(meta->general->bands, "01") == NULL ||
             strstr(meta->general->bands, "02") == NULL ||
             strstr(meta->general->bands, "03") == NULL)
          )
        {
          asfPrintError("Imported file does not contain required color bands\n"
                        "necessary for true color output (03, 02, 01)\n");
        }
        if (false_color &&
            (strstr(meta->general->bands, "02") == NULL ||
             strstr(meta->general->bands, "03") == NULL ||
             strstr(meta->general->bands, "04") == NULL)
          )
        {
          asfPrintError("Imported file does not contain required color bands\n"
                        "necessary for false color output (04, 03, 02)\n");
        }
        if (scale != SIGMA) {
          asfPrintWarning("A byte conversion other than SIGMA was specified.  Since\n"
                          "True Color or False Color output was selected, the byte conversion\n"
                          "will be overridden with SIGMA, a 2-sigma contrast expansion.\n");
        }
        char **bands = extract_band_names(meta->general->bands, meta->general->band_count);
        if (meta->general->band_count >= 4 && bands != NULL) {
          // The imported file IS a multiband file with enough bands,
          // but the extract bands need to be ordered correctly
          if (true_color) {
            strcpy(bands[0], "03");
            strcpy(bands[1], "02");
            strcpy(bands[2], "01");
            strcpy(bands[3], "");
          }
          else {
            strcpy(bands[0], "04");
            strcpy(bands[1], "03");
            strcpy(bands[2], "02");
            strcpy(bands[3], "");
          }
          check_return(asf_export_bands(format, SIGMA, TRUE,
                                        true_color, false_color, 0, 0, NULL,
                                        inFile, outFile, bands),
                       "exporting data file (asf_export), color banded.\n");
          for (i=0; i<meta->general->band_count; ++i)
            FREE (bands[i]);
          FREE(bands);
        }
        else {
          asfPrintError("Cannot determine band names from imported metadata file:\n  %s\n",
                        inFile);
        }
      }
      else {
        if (meta->general->band_count != 1 && strlen(cfg->export->band) == 0) {
          // multi-band, exporting as separate greyscale files
          asfPrintStatus("\nExporting %d bands as separate greyscale files...\n",
                         meta->general->band_count);
          check_return(asf_export_bands(format, scale, FALSE,
                                        0, 0, 0, 0, NULL,
                                        inFile, outFile, NULL),
                       "exporting data file (asf_export), greyscale bands.\n");
        }
        else if (meta->general->band_count != 1 && strlen(cfg->export->band) > 0) {
          // multi-band, exporting single band in one greyscale file
          asfPrintStatus("\nExporting band \"%s\" in a single greyscale file...\n",
                         cfg->export->band);
          int num_bands_found = 0;
          char **band_names = find_single_band(inFile, cfg->export->band,
                                               &num_bands_found);
          if (num_bands_found != 1) {
            asfPrintError("Selected band for export not found.\n");
          }
          check_return(asf_export_bands(format, scale, FALSE,
                                        0, 0, 0, 0, NULL,
                                        inFile, outFile, band_names),
                       "exporting data file (asf_export), single selected greyscale band.\n");
          if (*band_names) FREE(*band_names);
          if (band_names) FREE(band_names);
        }
        else {
          // single band
          check_return(asf_export(format, scale, inFile, outFile),
                       "exporting data file (asf_export), single band.\n");
        }
      }
    }
  }
  meta_free(meta);
}

