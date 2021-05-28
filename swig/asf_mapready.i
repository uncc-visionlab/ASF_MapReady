/* file : asf_mapready.i */

/* name of module to use*/
//%module(package="asf_mapready") asf_mapready
%module asf_mapready
%{
    /* Every thing in this file is being copied in  
     wrapper file. We include the C header file necessary 
     to compile the interface */
    // in src/asf/asf.h
    //#define MAGIC_UNSET_CHAR '?'
    //#define MAGIC_UNSET_STRING "???"
#include "../include/calibrate.h" 
#include "../include/asf_meta.h"
#include "../include/asf.h"
#include "../include/asf_import.h" 
#include "../include/calibrate.h" 
#include "../include/ceos.h"
#include "../include/ceos_io.h"
#include "../include/get_ceos_names.h"
#include "../include/geolocate.h"
#include "../src/libasf_import/decoder.h"
    /* variable declaration*/
    void ERS_readNextPulse(bin_state *s, iqType *iqBuf, char *inName, char *outName);
    void JRS_readNextCeosPulse(bin_state *s, iqType *iqBuf, char *inName, char *outName);
    void RSAT_readNextPulse(bin_state *s,iqType *iqBuf, char *inName, char *outName);
    bin_state * convertMetadata_ceos(char *inN, char *outN, int *nLines, readPulseFunc * readNextPulse);
    meta_parameters *import_ceos_raw(char *inDataName, char *inMetaName, char *outDataName,
         char *outMetaName, char *bandExt, int band, int nBands,
         radiometry_t radiometry, int line, int sample, int width, int height,
         int import_single_band);    
    //double myvar; 
    //extern META_DDR_STRUCT meta_ddr_structs[NUM_META_DDR_STRUCTS];
%}
%include "../include/calibrate.h"
%include "../include/asf_meta.h"
%include "../include/asf.h"
%include "../include/ceos.h"
%include "../include/ceos_io.h"
%include "../include/get_ceos_names.h"
%include "../include/geolocate.h"
//%ignore new_bin_state;
//%ignore delete_bin_state; 
%rename new_bin_state new_bin_state2;
%rename delete_bin_state delete_bin_state2;
%include "../src/libasf_import/decoder.h"

%include cpointer.i
%pointer_functions(int, intp);
%pointer_functions(readPulseFunc, readPulseFuncp);
%pointer_functions(meta_parameters, meta_parametersp);
%pointer_functions(double, doublep);
%include <carrays.i>
%array_functions(int, int_array1d);
%array_functions(double, double_array1d);
%array_functions(int *, int_array2d);
%array_functions(double *, double_array2d);
/* Some callback functions */
//%callback("%(uppercase)s");
//void ERS_readNextPulse(bin_state *s, iqType *iqBuf, char *inName, char *outName);
//void JRS_readNextPulse(bin_state *s, iqType *iqBuf, char *inName, char *outName);
//void RSAT_readNextPulse(bin_state *s,iqType *iqBuf, char *inName, char *outName);
//%nocallback;
//%constant(void (*readPulseFunc)(bin_state *s,iqType *iqBuf, char *inN, char *outN)) EERS_readNextCeosPulse = ERS_readNextCeosPulse;
//% constant void ERS_readNextPulse(bin_state *s, iqType *iqBuf, char *inName, char *outName);
//%ignore new_bin_state;
//%ignore delete_bin_state; 
//%rename(new_bin_state) renamed_new_bin_state; 
//%rename(new_bin_state) nnew_bin_state; 

/* explicitly list functions and variables to be interfaced */
//double myvar; 
int asf_import(radiometry_t radiometry, // r_AMP,R_SIGMA,r_BETA,r_GAMMA,r_POWER
        int db_flag, // TRUE if the output should be in decibels
        // only ok for radiometry=SIGMA,GAMMA,BETA
        int complex_flag, // TRUE if ingested SLC should in I/Q
        int multilook_flag, // TRUE is SLC should be multilooked while
        // being ingested.
        int azimuth_look_count,
        int range_look_count,
        int amp0_flag, // TRUE if we should generate a band 0
        // with amplitude data, with the rest
        // starting at band 1
        input_format_t format_type, // eg, STF, CEOS - etc
        char *band_id, // eg, "" (default for all bands), "VH", "03" - etc
        char *data_type, // data type for gamma ingest
        char *image_data_type, // "geocoded_image", "dem", or "mask"
        char *lutName, // NULL for no lookup table
        // otherwise, this is the lookup table filename
        char *prcPath, // NULL for not using precision orbit data
        // otherwise, this is the precision state vector
        // path
        double lowerLat, // -99 means not constrained
        double upperLat, // -99 means not constrained
        double lowerLon,
        double upperLon,
        int line, // start line subset - default set to 0
        int sample, // start sample subset - default set to 0
        int width, // -99 means no subsetting
        int height, // -99 means no subsetting
        int save_intermediates, // save intermediate files if any were created
        double *p_range_scale, // NULL for no scaling
        double *p_azimuth_scale, // NULL for no scaling
        double *p_correct_y_pixel_size, // NULL for no fixing
        int apply_ers2_gain_fix, // TRUE for correction of ers2 data
        char *inMetaNameOption, // NULL for normal metadata naming
        // otherwise, this is the meta file name
        char *inBaseName, // input file
        char *ancillary_file, // ancillary file (if needed for input file)
        char *colormapName, // colormap file
        char *slave_file, // slave metadata file
        char *interferogram_file, // interferogram file
        char *coherence_file, // coherence image file
        char *baseline_file, // baseline file
        int complex_gamm_file, // TRUE for complex GAMMA file
        char *uavsar_type, // data type for UAVSAR data
        int metaonly, // flat for generating XML metadata file only
        char *outBaseName // output file
        );
void ERS_readNextPulse(bin_state *s, iqType *iqBuf, char *inName, char *outName);
void JRS_readNextCeosPulse(bin_state *s, iqType *iqBuf, char *inName, char *outName);
void RSAT_readNextPulse(bin_state *s,iqType *iqBuf, char *inName, char *outName);
meta_parameters *import_ceos_raw(char *inDataName, char *inMetaName, char *outDataName,
         char *outMetaName, char *bandExt, int band, int nBands,
         radiometry_t radiometry, int line, int sample, int width, int height,
     int import_single_band);
bin_state *convertMetadata_ceos(char *inN, char *outN, int *nLines, readPulseFunc *readNextPulse);

/* or if we want to interface all functions then we can simply 
   include header file like this -  
   %include "gfg.h" 
 */
