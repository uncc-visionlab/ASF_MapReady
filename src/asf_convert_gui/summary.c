#include "asf_convert_gui.h"

/*
int zone(float lon)
{
return((int)(((lon + 180.0) / 6.0) + 1.0));
}
*/

SIGNAL_CALLBACK
void update_summary()
{
    GtkWidget * summary_label;
    Settings * s;
    gchar text[1024];
    gchar * type;

    s = settings_get_from_gui();

    strcpy(text, "Import: ");

    switch (s->data_type)
    {
        case INPUT_TYPE_SIGMA:
            type = "Sigma";
            break;
            
        case INPUT_TYPE_BETA:
            type = "Beta";
            break;
            
        case INPUT_TYPE_GAMMA:
            type = "Gamma";
            break;
            
        default:
        case INPUT_TYPE_AMP:
            type = "Amplitude";
            break;
            
        case INPUT_TYPE_POWER:
            type = "Power";
            break;
    }
    
    switch (s->input_data_format)
    {
    case INPUT_FORMAT_CEOS_LEVEL0:
    {
        GtkWidget *process_to_level1_checkbutton =
            glade_xml_get_widget(glade_xml, "process_to_level1_checkbutton");

        gboolean process_to_level1_is_checked =
            gtk_toggle_button_get_active(
                GTK_TOGGLE_BUTTON(process_to_level1_checkbutton));

        strcat(text, "CEOS Level Zero");

        if (process_to_level1_is_checked)
        {
            sprintf(text, "%s\n   Process to Level 1\nData Type: %s",
                    text, type);
        }
    }    
    break;

    default:
    case INPUT_FORMAT_CEOS_LEVEL1:
        sprintf(text, "%sCEOS Level One\nData type: %s",
            text, type);

        break;

    case INPUT_FORMAT_STF:
        strcat(text, "STF");
        if (s->latitude_checked)
        {
            sprintf(text, "%s (Lat: %.2f - %.2f)", text, s->latitude_low,
                s->latitude_hi);
        }
        break;

    case INPUT_FORMAT_COMPLEX:
        strcat(text, "Complex");
        break;
    }

    if (s->terrcorr_is_checked || s->refine_geolocation_is_checked)
    {
        GtkWidget *dem_entry;
	char *dem;

        if (s->terrcorr_is_checked)
            strcat(text, "\nTerrain Correction: Yes");
        else
            strcat(text, "\nRefine Geolocation: Yes");

	dem_entry = glade_xml_get_widget(glade_xml, "dem_entry");
	dem = strdup(gtk_entry_get_text(GTK_ENTRY(dem_entry)));

	if (!dem || strlen(dem) == 0)
	{
	    strcat(text, "\n   DEM: <none>");
	}
	else
	{
	    char *p = strrchr(dem, DIR_SEPARATOR);
	    sprintf(text, "%s\n   DEM: %s", text, p ? p+1 : dem);
	}

	free(dem);

        if (s->terrcorr_is_checked)
        {
            if (s->specified_tc_pixel_size)
            {
                sprintf(text, "%s\n   Pixel Size: %.2f m", text,
                        s->tc_pixel_size);
            }
            else if (s->specified_pixel_size)
            {
                sprintf(text, "%s\n   Pixel Size: %.2f m (from geocode)",
                        text, s->pixel_size);
            }

            sprintf(text, "%s\n   Interpolate Layover: %s", text,
                    s->interp ? "Yes" : "No");
        }

        if (s->mask_file_is_checked)
        {
            GtkWidget *mask_entry;
            char *mask;

            mask_entry = glade_xml_get_widget(glade_xml, "mask_entry");
            mask = strdup(gtk_entry_get_text(GTK_ENTRY(mask_entry)));

            if (!mask || strlen(mask) == 0)
            {
                strcat(text, "\n   Mask: <none>");
            }
            else
            {
                char *p = strrchr(mask, DIR_SEPARATOR);
                sprintf(text, "%s\n   Mask: %s", text, p ? p+1 : mask);
            }
        }
        else if (s->auto_water_mask_is_checked)
        {
            strcat(text, "\n   Automatic Mask");
        }

        if (s->generate_dem && s->generate_layover_mask)
        {
            strcat(text, "\n   Save DEM & Layover Mask");
        }
        else if (s->generate_dem)
        {
            strcat(text, "\n   Save DEM");
        }
        else if (s->generate_layover_mask)
        {
            strcat(text, "\n   Save Layover Mask");
        }
    }
    else
    {
        strcat(text, "\nTerrain Correction: No");
    }

    strcat(text, "\nGeocoding: ");

    if (s->geocode_is_checked)
    {
        switch (s->projection)
        {
        case PROJ_UTM:
	    if (s->zone != 0)
	    {
                sprintf(text, "%sUTM\n   Zone: %d\n",
                    text, s->zone);
	    }
	    else
	    {
                sprintf(text, "%sUTM\n   Zone: <from metadata>\n", text);
	    }
            break;

        case PROJ_PS:
            sprintf(text, "%sPolar Stereo\n"
                "   Center: (%.2f, %.2f)\n",
                text, s->lat0, s->lon0);
            break;

        case PROJ_LAMCC:
            sprintf(text, "%sLambert Conformal Conic\n"
                "   Center: (%.2f, %.2f)\n"
                "   Standard Parallels: (%.2f, %.2f)\n",
                text, s->lat0, s->lon0,
                s->plat1, s->plat2);
            break;

        case PROJ_LAMAZ:
            sprintf(text, "%sLambert Azimuthal Equal Area\n"
                "   Center: (%.2f, %.2f)\n",
                text, s->lat0, s->lon0);
            break;

        case PROJ_ALBERS:
            sprintf(text, "%sAlbers Conical Equal Area\n"
                "   Center: (%.2f, %.2f)\n"
                "   Standard Parallels: (%.2f, %.2f)\n",
                text, s->lat0, s->lon0,
                s->plat1, s->plat2);
            break;
        }

        if (s->specified_height)
            sprintf(text, "%s   Height: %.2f\n", text, s->height);

        if (s->specified_pixel_size)
            sprintf(text, "%s   Pixel Size: %.2f m\n", text, s->pixel_size);

        sprintf(text, "%s   Datum: %s\n", text, datum_string(s->datum));

        sprintf(text, "%s   Resampling Method: %s\n", text,
            resample_method_string(s->resample_method));
    }
    else
    {
        strcat(text, "<none>\n");
    }

    strcat(text, "Export: ");
    if (s->export_is_checked)
    {
        switch (s->output_format)
        {
        case OUTPUT_FORMAT_JPEG:
        default:
            strcat(text, "JPEG");
            break;

        case OUTPUT_FORMAT_PPM:
            strcat(text, "PPM");
            break;

        case OUTPUT_FORMAT_GEOTIFF:
            strcat(text, "geoTIFF");
            break;

        case OUTPUT_FORMAT_TIFF:
            strcat(text, "TIFF");
            break;
        }

        if (s->apply_scaling)
        {
            GtkWidget *longest_dimension_spinbutton;
            gdouble d;

            longest_dimension_spinbutton = 
                glade_xml_get_widget(glade_xml, 
                "longest_dimension_spinbutton");

            d = gtk_spin_button_get_value(
                GTK_SPIN_BUTTON(longest_dimension_spinbutton));

            sprintf(text, "%s (size %d)", text, (int)floor(d + 0.5));
        }

        strcat(text, "\n");
    }
    else
    {
        strcat(text, "<none>\n");
    }

    if (s->output_bytes)
    {
        strcat(text, "   Byte Scaling Method: ");
        switch (s->scaling_method)
        {
        default:
        case SCALING_METHOD_SIGMA:
            strcat(text, "Sigma");
            break;

        case SCALING_METHOD_MINMAX:
            strcat(text, "MinMax");
            break;

        case SCALING_METHOD_TRUNCATE:
            strcat(text, "Truncate");
            break;

        case SCALING_METHOD_HISTOGRAM_EQUALIZE:
            strcat(text, "Histogram Equalize");
            break;
        }
        //strcat(text, "\n");
    }

    summary_label =
        glade_xml_get_widget(glade_xml, "summary_label");

    gtk_label_set_text(GTK_LABEL(summary_label), text);

    settings_delete(s);
}
