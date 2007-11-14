#ifndef ASF_TERRCORR_H
#define ASF_TERRCORR_H

#include <asf_meta.h>

/**
   asf_terrcorr

   sarFile: Input SAR file
   demFile: Input DEM
   outFile: Output filename
   pixel_size: Desired pixel size, if negative will be left
               alone (or changed to match DEM if required
**/

int asf_terrcorr(char *sarFile, char *demFile,char *inMaskFile,
		 char *outFile, double pixel_size);

int refine_geolocation(char *sarFile, char *demFile, char *userMaskFile,
                       char *outFile, int update_flag, int auto_water_mask,
                       float mask_height_cutoff, int clean_files,
                       char **other_files_to_update_with_offsets);

int asf_check_geolocation(char *sarFile, char *demFile, char *inMaskFile,
			  char *simAmpFile, char *demSlant);

/**
   asf_terrcorr_ext

      Extended version of the above, with more options.  The above uses
      this one with sensible defaults.

   sarFile: Input SAR file
   demFile: Input DEM
   outFile: Output filename
   pixel_size: Desired pixel size, if negative will be left
               alone (or changed to match DEM if required
**/

int asf_terrcorr_ext(char *sarFile, char *demFile, char *inMaskFile,
                     char *outFile, double pixel_size, int clean_files,
                     int do_resample, int do_corner_matching, int do_interp,
                     int do_fftMatch_verification, int dem_grid_size, 
                     int do_terrain_correction, int fill_value,
                     int generate_water_mask, int save_clipped_dem,
                     int update_original_metadata_with_offsets,
                     float mask_height_cutoff, int doRadiometric,
                     int smooth_dem_holes,
                     char **other_files_to_update_with_offsets,
		     int no_matching, double range_offset,
		     double azimuth_offset);

void
clip_dem(meta_parameters *metaSAR, char *srFile, char *demFile,
         char *demClipped, char *what, char *otherFile, char *otherClipped,
         char *otherWhat, char *output_dir, int dem_grid_size,
         int clean_files, int *p_demHeight);

/**
   Functions private to terrain correction, not meant for general use.
**/

/* Prototypes from seedsquares.c */
int lay_seeds(int num_seeds, float *mask, long ns, long nl,
              int *x_tl_list, int *y_tl_list, 
              int *x_br_list, int *y_br_list,
              float *good_pct_list);

/* Prototypes from build_dem.c */
char *build_dem(meta_parameters *meta, const char *dem_cla_arg,
                const char *dir_for_tmp_dem);
int get_dem_chunk(char *dem_in, char *dem_out, meta_parameters *metaDEM,
                  meta_parameters *metaSAR);

#endif
