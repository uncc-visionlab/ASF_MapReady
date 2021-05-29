import asf_mapready

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
def parseCEOS(inMetaName):
    dssr = asf_mapready.dataset_sum_rec();
    asf_mapready.get_dssr(inMetaName, dssr);
    return dssr;

def parseBinState(inMetaName):
    outMetaName = "testout.meta";
    # SWIGTYPE_p_int nLines = SWIGTYPE_p_int();        
    #SWIGTYPE_p_int nLines = asf_mapready.new_intp();
    # SWIGTYPE_p_p_f_p_bin_state_p_unsigned_char_p_char_p_char__void readNextPulse = SWIGTYPE_p_p_f_p_bin_state_p_unsigned_char_p_char_p_char__void();
    #SWIGTYPE_p_p_f_p_bin_state_p_unsigned_char_p_char_p_char__void readNextPulse = asf_mapready.new_readPulseFuncp();
    nLines = asf_mapready.new_intp();
    readNextPulse = asf_mapready.new_readPulseFuncp();
    s = asf_mapready.convertMetadata_ceos(inMetaName, outMetaName, nLines, readNextPulse);
    return s;

if __name__ == "__main__":
    inFile = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000.ldr"; # NULL for normal metadata naming
    dssr = parseCEOS(inFile);

    darray = dssr.crt_dopcen;
    print(darray);
    crt_dopcen = [asf_mapready.double_array1d_getitem(darray, 0),
        asf_mapready.double_array1d_getitem(darray, 1),
            asf_mapready.double_array1d_getitem(darray, 2)];
    print("dopCent = {", crt_dopcen[0], " , ", crt_dopcen[1], " , ", crt_dopcen[2], "}");
    inRawName = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000.raw"; # input file
    s = parseBinState(inRawName);
    print(s);
    
    radiometry = asf_mapready.r_AMP; # r_AMP,R_SIGMA,r_BETA,r_GAMMA,r_POWER
    db_flag = 0;   # TRUE if the output should be in decibels
    # only ok for radiometry=SIGMA,GAMMA,BETA
    complex_flag = 0; # TRUE if ingested SLC should in I/Q
    multilook_flag = 0; # TRUE is SLC should be multilooked while
    # being ingested.
    azimuth_look_count = 0;
    range_look_count = 0;
    amp0_flag = 0;      # TRUE if we should generate a band 0
    # with amplitude data, with the rest
    # starting at band 1
    format_type = asf_mapready.CEOS; # eg, STF, CEOS - etc
    band_id = ""; # eg, "" (default for all bands), "VH", "03" - etc
    data_type = "";    # data type for gamma ingest
    image_data_type = ""; # "geocoded_image", "dem", or "mask"
    lutName = ""; # NULL for no lookup table
    # otherwise, this is the lookup table filename
    prcPath = ""; # NULL for not using precision orbit data
    # otherwise, this is the precision state vector
    # path
    lowerLat = 0; # -99 means not constrained
    upperLat = 0; # -99 means not constrained
    lowerLon = 0;
    upperLon = 0;
    line = 0; # start line subset - default set to 0
    sample = 0; # start sample subset - default set to 0
    width = 0; # -99 means no subsetting
    height = 0; # -99 means no subsetting
    save_intermediates = 0; # save intermediate files if any were created
    p_range_scale = asf_mapready.new_doublep(); # NULL for no scaling
    p_azimuth_scale = asf_mapready.new_doublep(); # NULL for no scaling
    p_correct_y_pixel_size = asf_mapready.new_doublep(); # NULL for no fixing
    apply_ers2_gain_fix = 0; # TRUE for correction of ers2 data
    inMetaNameOption = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000.ldr"; # NULL for normal metadata naming
    # otherwise, this is the meta file name
    inBaseName = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000"; # input file
    ancillary_file = ""; # ancillary file (if needed for input file)
    colormapName = ""; # colormap file
    slave_file = ""; # slave metadata file
    interferogram_file = ""; # interferogram file
    coherence_file = ""; # coherence image file
    baseline_file = ""; # baseline file
    complex_gamma_file = 0; # TRUE for complex GAMMA file
    uavsar_type = ""; # data type for UAVSAR data
    metaonly = 0; # flat for generating XML metadata file only
    outBaseName = "test"; # output file

    asf_mapready.asf_import(radiometry, db_flag, complex_flag, multilook_flag,
                            azimuth_look_count, range_look_count,
                            amp0_flag, format_type, band_id, data_type, image_data_type,
                            lutName, prcPath, lowerLat, upperLat, lowerLon, upperLon,
                            line, sample, width, height, save_intermediates,
                            p_range_scale, p_azimuth_scale, p_correct_y_pixel_size,
                            apply_ers2_gain_fix, inMetaNameOption, inBaseName,
                            ancillary_file, colormapName, slave_file, interferogram_file,
                            coherence_file, baseline_file, complex_gamma_file,
                            uavsar_type, metaonly, outBaseName);