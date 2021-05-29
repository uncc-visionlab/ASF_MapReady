/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import swig.asf.*;

/**
 *
 * @author arwillis
 */
public class ASFMapReadyJava {

    public ASFMapReadyJava() throws IOException {
        loadNativeLibraryFromJar();
    }

    public static void loadNativeLibraryFromJar() throws IOException {
        File file;
        OutputStream os;
        System.out.println("Loading library ....");
        try ( InputStream is = ASFMapReadyJava.class.getResourceAsStream("/native_libs/libasf_mapready.so")) {
            file = File.createTempFile("lib", ".so");
            os = new FileOutputStream(file);
            byte[] buffer = new byte[1024];
            int length;
            while ((length = is.read(buffer)) != -1) {
                os.write(buffer, 0, length);
            }
        }
        os.close();

        System.load(file.getAbsolutePath());
        file.deleteOnExit();
        System.out.println("Library loaded.");
    }

    public static dataset_sum_rec parseCEOS(String inMetaName) {
        dataset_sum_rec dssr = new dataset_sum_rec();
        asf_mapready.get_dssr(inMetaName, dssr);
        return dssr;
    }

    public static bin_state parseBinState(String inMetaName) {
        String outMetaName = "testout.meta";
        bin_state s = new bin_state();
        //SWIGTYPE_p_int nLines = SWIGTYPE_p_int();        
        SWIGTYPE_p_int nLines = asf_mapready.new_intp();
        //SWIGTYPE_p_p_f_p_bin_state_p_unsigned_char_p_char_p_char__void readNextPulse = SWIGTYPE_p_p_f_p_bin_state_p_unsigned_char_p_char_p_char__void();
        SWIGTYPE_p_p_f_p_bin_state_p_unsigned_char_p_char_p_char__void readNextPulse = asf_mapready.new_readPulseFuncp();
        s = asf_mapready.convertMetadata_ceos(inMetaName, outMetaName, nLines, readNextPulse);
        return s;
    }

    public static meta_parameters asfImportCEOS(String inBaseName) {
        String inDataName = inBaseName + ".raw";
        String inMetaName = inBaseName + ".ldr";
        String outDataName = "out" + ".img";
        String outMetaName = "out" + ".meta";
        String bandExt = "";
        int band = 0;
        int nBands = 0;
        int line = 0;
        int sample = 0;
        int width = 0;
        int height = 0;
        int import_single_band = 0;
        radiometry_t radiometry = radiometry_t.r_AMP; // r_AMP,R_SIGMA,r_BETA,r_GAMMA,r_POWER
        meta_parameters meta =  asf_mapready.import_ceos_raw(inDataName, inMetaName, outDataName,
                outMetaName, bandExt, band, nBands, radiometry, line, sample, width, height,
                import_single_band);
        return meta;
    }

    public static void main(String[] args) throws IOException {
        loadNativeLibraryFromJar();

        String inFile = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000.ldr"; // NULL for normal metadata naming
        dataset_sum_rec dssr = parseCEOS(inFile);

        SWIGTYPE_p_double darray = dssr.getCrt_dopcen();
        double[] crt_dopcen = {asf_mapready.double_array1d_getitem(darray, 0),
            asf_mapready.double_array1d_getitem(darray, 1),
            asf_mapready.double_array1d_getitem(darray, 2)};
        System.out.println("dopCent = {" + crt_dopcen[0] + " , " + crt_dopcen[1] + " , " + crt_dopcen[2] + "}");

        String inRawName = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000.raw"; // input file
        bin_state s = parseBinState(inRawName);

        radiometry_t radiometry = radiometry_t.r_AMP; // r_AMP,R_SIGMA,r_BETA,r_GAMMA,r_POWER
        int db_flag = 0;   // TRUE if the output should be in decibels
        // only ok for radiometry=SIGMA,GAMMA,BETA
        int complex_flag = 0; // TRUE if ingested SLC should in I/Q
        int multilook_flag = 0; // TRUE is SLC should be multilooked while
        // being ingested.
        int azimuth_look_count = 0;
        int range_look_count = 0;
        int amp0_flag = 0;      // TRUE if we should generate a band 0
        // with amplitude data, with the rest
        // starting at band 1
        input_format_t format_type = input_format_t.CEOS; // eg, STF, CEOS - etc
        String band_id = ""; // eg, "" (default for all bands), "VH", "03" - etc
        String data_type = "";    // data type for gamma ingest
        String image_data_type = ""; // "geocoded_image", "dem", or "mask"
        String lutName = ""; // NULL for no lookup table
        // otherwise, this is the lookup table filename
        String prcPath = ""; // NULL for not using precision orbit data
        // otherwise, this is the precision state vector
        // path
        double lowerLat = 0; // -99 means not constrained
        double upperLat = 0; // -99 means not constrained
        double lowerLon = 0;
        double upperLon = 0;
        int line = 0; // start line subset - default set to 0
        int sample = 0; // start sample subset - default set to 0
        int width = 0; // -99 means no subsetting
        int height = 0; // -99 means no subsetting
        int save_intermediates = 0; // save intermediate files if any were created
        SWIGTYPE_p_double p_range_scale = null; // NULL for no scaling
        SWIGTYPE_p_double p_azimuth_scale = null; // NULL for no scaling
        SWIGTYPE_p_double p_correct_y_pixel_size = null; // NULL for no fixing
        int apply_ers2_gain_fix = 0; // TRUE for correction of ers2 data
        String inMetaNameOption = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000.ldr"; // NULL for normal metadata naming
        // otherwise, this is the meta file name
        String inBaseName = "/media/arwillis/Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000"; // input file
        String ancillary_file = ""; // ancillary file (if needed for input file)
        String colormapName = ""; // colormap file
        String slave_file = ""; // slave metadata file
        String interferogram_file = ""; // interferogram file
        String coherence_file = ""; // coherence image file
        String baseline_file = ""; // baseline file
        int complex_gamma_file = 0; // TRUE for complex GAMMA file
        String uavsar_type = ""; // data type for UAVSAR data
        int metaonly = 0; // flat for generating XML metadata file only
        String outBaseName = "test"; // output file

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
    }
}
