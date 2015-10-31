#include "shapefil.h"
#include "asf_vector.h"
#include "meta_project.h"
#include "asf.h"
#include <assert.h>
#include <ctype.h>
#include "asf_nan.h"
#include "dateUtil.h"
#include "float_image.h"
#include "ceos_io.h"
#include "libasf_proj.h"
#include <stdio.h>
#include <math.h>

int point2kml(char *inFile, char *outFile)
{
  char id[255], line[1024];
  double lat, lon;
  FILE *fpIn = FOPEN(inFile, "r");
  FILE *fpOut = FOPEN(outFile, "w");
  kml_header(fpOut);
  fgets(line, 1024, fpIn); // header line
  while (fgets(line, 1024, fpIn)) {
    chomp(line);
    sprintf(id, "%s", get_str(line, 0));
    lat = atof(get_str(line, 1));
    lon = atof(get_str(line, 2));

    // Write information in kml file
    fprintf(fpOut, "<Placemark>\n");
    fprintf(fpOut, "  <name>%s</name>\n", id);
    fprintf(fpOut, "  <description>\n");
    fprintf(fpOut, "    <![CDATA[\n");
    fprintf(fpOut, "<!-- Format: POINT (generated by %s) -->",
      TOOL_SUITE_VERSION_STRING);
    fprintf(fpOut, "      <strong>ID</strong>: %s <br>\n", id);
    fprintf(fpOut, "      <strong>Latitude</strong>: %12.5f <br>\n", lat);
    fprintf(fpOut, "      <strong>Longitude</strong>: %12.5f <br>\n", lon);
    fprintf(fpOut, "    ]]>\n");
    fprintf(fpOut, "  </description>\n");
    fprintf(fpOut, "  <LookAt>\n");
    fprintf(fpOut, "    <longitude>%12.5f</longitude>\n", lon);
    fprintf(fpOut, "    <latitude>%12.5f</latitude>\n", lat);
    fprintf(fpOut, "    <range>400000</range>\n");
    fprintf(fpOut, "  </LookAt>\n");
    fprintf(fpOut, "  <Point>\n");
    fprintf(fpOut, "    <coordinates>%f,%f,0</coordinates>\n", lon, lat);
    fprintf(fpOut, "  </Point>\n");
    fprintf(fpOut, "</Placemark>\n");
  }
  kml_footer(fpOut);
  FCLOSE(fpOut);
  FCLOSE(fpIn);

  return TRUE;
}

void kml_header(FILE *kml_file)
{
    fprintf(kml_file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(kml_file,
            "<kml xmlns=\"http://earth.google.com/kml/2.2\">\n");
    fprintf(kml_file, "<Document>\n");
}

void kml_footer(FILE *kml_file)
{
    fprintf(kml_file, "</Document>\n");
    fprintf(kml_file, "</kml>\n");
}

void write_kml_style_keys_ext(FILE *kml_file, c2v_config *cfg)
{
  // so all of our methods use the same "look" for the
  // boxes/lines.
  fprintf(kml_file, "  <Style>\n");
  fprintf(kml_file, "    <LineStyle>\n");
  if (cfg) {
    fprintf(kml_file, "      <color>%s</color>\n", cfg->color);
    fprintf(kml_file, "      <width>%d</width>\n", cfg->width);
  }
  else {
    fprintf(kml_file, "      <color>ffff9900</color>\n");
    fprintf(kml_file, "      <width>5</width>\n");
  }
  fprintf(kml_file, "    </LineStyle>\n");
  fprintf(kml_file, "    <PolyStyle>\n");
  fprintf(kml_file, "      <color>1fff5500</color>\n");
  fprintf(kml_file, "    </PolyStyle>\n");
  fprintf(kml_file, "  </Style>\n");
}

void write_kml_style_keys(FILE *kml_file)
{
  write_kml_style_keys_ext(kml_file, NULL);
}

void write_kml_attributes(FILE *kml_file, int nAttr, dbf_header_t *dbf)
{
  int ii;

  fprintf(kml_file, "  <description><![CDATA[\n");
  fprintf(kml_file, "  <table width=\"400\"><tr><td>\n");
  
  // Write fields into the KML file
  for (ii=0; ii<nAttr; ii++) {
    if (dbf[ii].format == DBF_STRING) {
      //printf("STR ii=%d nAttr=%d def=%s val=%s\n", ii, nAttr, dbf[ii].definition, dbf[ii].sValue);
      fprintf(kml_file, "  <strong>%s</strong>: %s <br>\n", 
        dbf[ii].definition, dbf[ii].sValue ? dbf[ii].sValue : "");
    } else if (dbf[ii].format == DBF_INTEGER) {
      //printf("INT ii=%d nAttr=%d val=%s\n", ii, nAttr, dbf[ii].definition);
      fprintf(kml_file, "  <strong>%s</strong>: %d <br>\n",
        dbf[ii].definition, dbf[ii].nValue);
    } else if (dbf[ii].format == DBF_DOUBLE) {
      //printf("DBL ii=%d nAttr=%d val=%s\n", ii, nAttr, dbf[ii].definition);
      fprintf(kml_file, "  <strong>%s</strong>: %f <br>\n",
        dbf[ii].definition, dbf[ii].fValue);
    }
  }
  
  fprintf(kml_file, "  </td></tr></table>\n");
  fprintf(kml_file, "  ]]></description>\n");
}

void write_kml_object(FILE *kml_file, c2v_config *cfg, char *png_filename,
  int start, int end, double *lat, double *lon)
{
  int ii;
  if (cfg) { 
    if (strcmp_case(cfg->boundary, "POLYGON") != 0 &&
        strcmp_case(cfg->boundary, "LINE") != 0)
      strcpy(cfg->boundary, "POLYGON");
    if (strcmp_case(cfg->altitude, "RELATIVETOGROUND") != 0 &&
        strcmp_case(cfg->altitude, "CLAMPTOGROUND") != 0 &&
        strcmp_case(cfg->altitude, "ABSOLUTE") != 0)
      strcpy(cfg->altitude, "absolute");
  }
  if (cfg && strcmp_case(cfg->boundary, "LINE") == 0) {
    write_kml_style_keys_ext(kml_file, cfg);
    fprintf(kml_file, "  <LineString>\n");
    fprintf(kml_file, "    <altitudeMode>%s</altitudeMode>\n", cfg->altitude);
    fprintf(kml_file, "    <extrude>1</extrude>\n");
    fprintf(kml_file, "    <tesselate>1</tesselate>\n");
  }
  else {
    write_kml_style_keys_ext(kml_file, cfg);
    fprintf(kml_file, "  <Polygon>\n");
    fprintf(kml_file, "    <extrude>%d</extrude>\n",
      png_filename ? 0 : 1);
    fprintf(kml_file, "    <altitudeMode>%s</altitudeMode>\n", 
            cfg ? cfg->altitude : "3000");
    fprintf(kml_file, "    <outerBoundaryIs>\n");
    fprintf(kml_file, "      <LinearRing>\n");
  }
  fprintf(kml_file, "        <coordinates>\n");
  int height = cfg && cfg->height > 0 ? cfg->height : 7000;
  for (ii=start; ii<end; ii++)
    fprintf(kml_file, "          %.6f,%.6f,%d\n", lon[ii], lat[ii], height);
  fprintf(kml_file, "        </coordinates>\n");
  if (cfg && strcmp_case(cfg->boundary, "LINE") == 0)
    fprintf(kml_file, "  </LineString>\n");
  else {
    fprintf(kml_file, "      </LinearRing>\n");
    fprintf(kml_file, "    </outerBoundaryIs>\n");
    fprintf(kml_file, "  </Polygon>\n");
  }
}

void write_kml_placemark(FILE *kml_file, char *name, double center_lat,
  double center_lon, char *png_filename, dbf_header_t *dbf, int nAttr,
  double *lat, double *lon, int nCoords, c2v_config *cfg)
{
  fprintf(kml_file, "<Placemark>\n");

  write_kml_attributes(kml_file, nAttr, dbf);
  
  fprintf(kml_file, "  <name>%s</name>\n", name);
  fprintf(kml_file, "  <LookAt>\n");
  fprintf(kml_file, "    <longitude>%.6f</longitude>\n", center_lon);
  fprintf(kml_file, "    <latitude>%.6f</latitude>\n", center_lat);
  fprintf(kml_file, "    <range>%d</range>\n", cfg ? cfg->range : 400000);
  fprintf(kml_file, "  </LookAt>\n");
  
  // Check whether we need to split up the polygon
  if (crosses_dateline(lon, 0, nCoords)) {
    int *start = (int *) MALLOC(sizeof(int)*2);
    double *mLat = (double *) MALLOC(sizeof(double)*(nCoords+5));
    double *mLon = (double *) MALLOC(sizeof(double)*(nCoords+5));

    split_polygon(lat, lon, 5, start, mLat, mLon, 60);
    fprintf(kml_file, "  <MultiGeometry>\n");
    write_kml_object(kml_file, cfg, png_filename, 0, start[1], mLat, mLon);
    write_kml_object(kml_file, cfg, png_filename, start[1], 10, mLat, mLon);
    fprintf(kml_file, "  </MultiGeometry>\n");
    
    FREE(mLat);
    FREE(mLon);
    FREE(start);
  }
  else
    write_kml_object(kml_file, cfg, png_filename, 0, 5, lat, lon);
  fprintf(kml_file, "</Placemark>\n");
    
  if (png_filename) {
    fprintf(kml_file, "<GroundOverlay>\n");
    fprintf(kml_file, "  <name>%s</name>\n", name);
    if (cfg)
      fprintf(kml_file, "  <color>%xffffff</color>\n", 
          (int)((100 - cfg->transparency) * 2.55 + 0.5));
    fprintf(kml_file, "  <Icon>\n");
    fprintf(kml_file, "      <href>%s</href>\n", png_filename);
    fprintf(kml_file, "      <viewBoundScale>0.75</viewBoundScale>\n");
    fprintf(kml_file, "  </Icon>\n");
    if (cfg) {
      fprintf(kml_file, "  <LatLonBox>\n");
      fprintf(kml_file, "      <north>%.4f</north>\n", cfg->north);
      fprintf(kml_file, "      <south>%.4f</south>\n", cfg->south);
      fprintf(kml_file, "      <east>%.4f</east>\n", cfg->east);
      fprintf(kml_file, "      <west>%.4f</west>\n", cfg->west);
      fprintf(kml_file, "  </LatLonBox>\n");
    }
    fprintf(kml_file, "</GroundOverlay>\n");
  }  
}

void csv2kml(char *inFile, char *outFile, char *format, c2v_config *cfg)
{
  // Read header information
  FILE *fpIn;
  dbf_header_t *header;
  int ii, kk, nCols;
  char shape_type[25], **cols, str[10], line[1024];
  if (strcmp_case(format, "CSV") == 0) {
    fpIn = FOPEN(inFile, "r");
    fgets(line, 1024, fpIn); // header line
    chomp(line);
    split_into_array(line, ',', &nCols, &cols);
    header = (dbf_header_t *) MALLOC(sizeof(dbf_header_t)*nCols);
    strcpy(shape_type, "polygon");
    for (ii=0; ii<nCols; ii++) {
      header[ii].meta = STRDUP(cols[ii]);
      header[ii].shape = STRDUP(cols[ii]);
      header[ii].format = DBF_STRING;
      header[ii].length = 50;
      header[ii].decimals = 0;
      header[ii].definition = STRDUP(cols[ii]);
      header[ii].column = ii; 
    }
    FREE(cols);
  }
  else
    read_header_config(format, &header, &nCols, shape_type);

  // Figure out how many vertices we got
  int nColumns, nVertices = 0, nLat = 0, nLon = 0;
  fpIn = FOPEN(inFile, "r");
  fgets(line, 1024, fpIn); // header line
  chomp(line);
  split_into_array(line, ',', &nColumns, &cols);
  
  for (ii=0; ii<nCols; ii++) {
  
    // Assign the column we need to read from
    for (kk=0; kk<nColumns; kk++) {
      if (strcmp_case(header[ii].meta, cols[kk]) == 0)
        header[ii].column = kk;
    }
    
    // Assuming that we don't have more than 12 vertices
    for (kk=1; kk<13; kk++) {
      sprintf(str, "LAT%d", kk);
      if (strcmp_case(header[ii].shape, str) == 0)
        nLat++;
      sprintf(str, "LON%d", kk);
      if (strcmp_case(header[ii].shape, str) == 0)
        nLon++;
    }

    // Alternative column names - only things that we already know
    if (strcmp_case(header[ii].shape, "NEAR START LAT") == 0 ||
      strcmp_case(header[ii].shape, "NEAR_START_LAT") == 0)
      nLat++;
    else if (strcmp_case(header[ii].shape, "NEAR START LON") == 0 ||
      strcmp_case(header[ii].shape, "NEAR_START_LON") == 0)
      nLon++;
    else if (strcmp_case(header[ii].shape, "FAR START LAT") == 0 ||
      strcmp_case(header[ii].shape, "FAR_START_LAT") == 0)
      nLat++;
    else if (strcmp_case(header[ii].shape, "FAR START LON") == 0 ||
      strcmp_case(header[ii].shape, "FAR_START_LON") == 0)
      nLon++;
    else if (strcmp_case(header[ii].shape, "NEAR END LAT") == 0 ||
      strcmp_case(header[ii].shape, "NEAR_END_LAT") == 0)
      nLat++;
    else if (strcmp_case(header[ii].shape, "NEAR END LON") == 0 ||
      strcmp_case(header[ii].shape, "NEAR_END_LON") == 0)
      nLon++;
    else if (strcmp_case(header[ii].shape, "FAR END LAT") == 0 ||
      strcmp_case(header[ii].shape, "FAR_END_LAT") == 0)
      nLat++;
    else if (strcmp_case(header[ii].shape, "FAR END LON") == 0 ||
      strcmp_case(header[ii].shape, "FAR_END_LON") == 0)
      nLon++;
  }
  if (nLat != nLon)
    asfPrintError("Found %d latitude and %d longitude columns.\n"
      "Can't convert this information properly!\n", nLat, nLon);
  else {
    nVertices = nLat;
    asfPrintStatus("Found %d vertices of a polygon.\n", nVertices);
  }

  // Open KML file  
  FILE *kml_file = FOPEN(outFile, "w");
  kml_header(kml_file);
  fprintf(kml_file, "<!-- Format: %s (generated by %s) -->\n", format,
    TOOL_SUITE_VERSION_STRING);

  // Read polygon information
  double *lat = (double *) MALLOC(sizeof(double)*(nVertices+1));
  double *lon = (double *) MALLOC(sizeof(double)*(nVertices+1)); 
  int column;
  char name[100];
  double center_lat = 0.0, center_lon = 0.0;
  while (fgets(line, 1024, fpIn)) {
    chomp(line);
    split_into_array(line, ',', &nColumns, &cols);
    for (ii=0; ii<nCols; ii++) {
      column = header[ii].column;
      if (header[ii].format == DBF_STRING)
        header[ii].sValue = STRDUP(cols[column]);
      else if (header[ii].format == DBF_INTEGER)
        header[ii].nValue = atoi(cols[column]);
      else if (header[ii].format == DBF_DOUBLE)
        header[ii].fValue = atof(cols[column]);
        
      // Standard LAT/LON columns - mentioned in the documentation
      for (kk=0; kk<nVertices; kk++) {
        sprintf(str, "LAT%d", kk+1);
        if (strcmp_case(header[ii].shape, str) == 0)
          lat[kk] = atof(cols[column]);
        sprintf(str, "LON%d", kk+1);
        if (strcmp_case(header[ii].shape, str) == 0)
          lon[kk] = atof(cols[column]);
      }
      
      // Alternative column names - only things that we already know
      if (strcmp_case(header[ii].shape, "NEAR START LAT") == 0 ||
        strcmp_case(header[ii].shape, "NEAR_START_LAT") == 0)
        lat[0] = atof(cols[column]);
      else if (strcmp_case(header[ii].shape, "NEAR START LON") == 0 ||
        strcmp_case(header[ii].shape, "NEAR_START_LON") == 0)
        lon[0] = atof(cols[column]);
      else if (strcmp_case(header[ii].shape, "FAR START LAT") == 0 ||
        strcmp_case(header[ii].shape, "FAR_START_LAT") == 0)
        lat[1] = atof(cols[column]);
      else if (strcmp_case(header[ii].shape, "FAR START LON") == 0 ||
        strcmp_case(header[ii].shape, "FAR_START_LON") == 0)
        lon[1] = atof(cols[column]);
      else if (strcmp_case(header[ii].shape, "NEAR END LAT") == 0 ||
        strcmp_case(header[ii].shape, "NEAR_END_LAT") == 0)
        lat[3] = atof(cols[column]);
      else if (strcmp_case(header[ii].shape, "NEAR END LON") == 0 ||
        strcmp_case(header[ii].shape, "NEAR_END_LON") == 0)
        lon[3] = atof(cols[column]);
      else if (strcmp_case(header[ii].shape, "FAR END LAT") == 0 ||
        strcmp_case(header[ii].shape, "FAR_END_LAT") == 0)
        lat[2] = atof(cols[column]);
      else if (strcmp_case(header[ii].shape, "FAR END LON") == 0 ||
        strcmp_case(header[ii].shape, "FAR_END_LON") == 0)
        lon[2] = atof(cols[column]);
    }
    lat[nVertices] = lat[0];
    lon[nVertices] = lon[0];
    FREE(cols);

    // Writing KML placemark  
    sprintf(name, "%s", inFile);
    for (ii=0; ii<nVertices; ii++) {
      center_lat += lat[ii];
      center_lon += lon[ii];
    }
    center_lat /= nVertices;
    center_lon /= nVertices;
    write_kml_placemark(kml_file, name, center_lat, center_lon, NULL, 
      header, nCols, lat, lon, nVertices+1, cfg);
  }
  kml_footer(kml_file);
  FCLOSE(kml_file);
}

// Convert to KML file
int convert2kml(char *inFile, char *outFile, char *format, int list, 
  c2v_config *cfg)
{
  FILE *kml_file;
  dbf_header_t *dbf = NULL;
  char line[1024];
  int ii, n = 0, nAttr = 0, nCoords = 0;
  double *lat = NULL, *lon = NULL, center_lat = 0.0, center_lon = 0.0;
  meta_parameters *meta = NULL;
  char name[512];
  
  if (list) {
    FILE *fpIn = FOPEN(inFile, "r");
    kml_file = FOPEN(outFile, "w");
    kml_header(kml_file);
    while (fgets(line, 1024, fpIn)) {
      chomp(line);
      if (strcmp_case(format, "META") == 0) {
        meta = meta2vector(line, &dbf, &nAttr, &lat, &lon, &nCoords);
        sprintf(name, "%s", meta->general->basename);
        center_lat = meta->general->center_latitude;
        center_lon = meta->general->center_longitude;
        fprintf(kml_file, "<!-- Format: META (generated by %s) -->\n", 
          TOOL_SUITE_VERSION_STRING);
        write_kml_placemark(kml_file, name, center_lat, center_lon, NULL, 
          dbf, nAttr, lat, lon, nCoords, cfg);
        meta_free(meta);
      }
      else if (strcmp_case(format, "GEOTIFF") == 0) {
        geotiff2vector(inFile, &dbf, &nAttr, &lat, &lon, &nCoords);
        sprintf(name, "%s", inFile);
        for (ii=0; ii<nCoords; ii++) {
          center_lat += lat[ii];
          center_lon += lon[ii];
        }
        center_lat /= nCoords;
        center_lon /= nCoords;
        fprintf(kml_file, "<!-- Format: GeoTIFF (generated by %s) -->\n", 
          TOOL_SUITE_VERSION_STRING);
        write_kml_placemark(kml_file, name, center_lat, center_lon, NULL, 
          dbf, nAttr, lat, lon, nCoords, cfg);
      }
      else if (strcmp_case(format, "POLYGON") == 0) {
        polygon2vector(inFile, &dbf, &nAttr, &lat, &lon, &nCoords);
        sprintf(name, "%s", inFile);
        for (ii=0; ii<nCoords; ii++) {
          center_lat += lat[ii];
          center_lon += lon[ii];
        }
        center_lat /= nCoords;
        center_lon /= nCoords;
        kml_file = FOPEN(outFile, "w");
        kml_header(kml_file);
        fprintf(kml_file, "<!-- Format: Polygon (generated by %s) -->\n", 
          TOOL_SUITE_VERSION_STRING);
        write_kml_placemark(kml_file, name, center_lat, center_lon, NULL, 
          dbf, nAttr, lat, lon, nCoords, cfg);
        kml_footer(kml_file);
      }
      else
        asfPrintError("List option for %s format not supported\n", format);
      n++;
    }
    FCLOSE(fpIn);
    kml_footer(kml_file);
    FCLOSE(kml_file);
  }
  else {
    if (strcmp_case(format, "META") == 0) {
      meta = meta2vector(inFile, &dbf, &nAttr, &lat, &lon, &nCoords);
      sprintf(name, "%s", meta->general->basename);
      center_lat = meta->general->center_latitude;
      center_lon = meta->general->center_longitude;
      kml_file = FOPEN(outFile, "w");
      kml_header(kml_file);
      fprintf(kml_file, "<!-- Format: META (generated by %s) -->\n", 
        TOOL_SUITE_VERSION_STRING);
      if (strlen(cfg->overlay) > 0)
        write_kml_placemark(kml_file, name, center_lat, center_lon, 
          cfg->overlay, dbf, nAttr, lat, lon, nCoords, cfg);
      else
        write_kml_placemark(kml_file, name, center_lat, center_lon, NULL, 
          dbf, nAttr, lat, lon, nCoords, cfg);
      kml_footer(kml_file);
      meta_free(meta);
      FCLOSE(kml_file);
    }
    else if (strcmp_case(format, "GEOTIFF") == 0) {
      geotiff2vector(inFile, &dbf, &nAttr, &lat, &lon, &nCoords);
      sprintf(name, "%s", inFile);
      for (ii=0; ii<nCoords; ii++) {
        center_lat += lat[ii];
        center_lon += lon[ii];
      }
      center_lat /= nCoords;
      center_lon /= nCoords;
      kml_file = FOPEN(outFile, "w");
      kml_header(kml_file);
      fprintf(kml_file, "<!-- Format: GEOTIFF (generated by %s) -->\n", 
        TOOL_SUITE_VERSION_STRING);
      write_kml_placemark(kml_file, name, center_lat, center_lon, NULL, 
        dbf, nAttr, lat, lon, nCoords, cfg);
      kml_footer(kml_file);
      FCLOSE(kml_file);
    }
    else if (strcmp_case(format, "POLYGON") == 0) {
      polygon2vector(inFile, &dbf, &nAttr, &lat, &lon, &nCoords);
      sprintf(name, "%s", inFile);
      for (ii=0; ii<nCoords; ii++) {
        center_lat += lat[ii];
        center_lon += lon[ii];
      }
      center_lat /= nCoords;
      center_lon /= nCoords;
      kml_file = FOPEN(outFile, "w");
      kml_header(kml_file);
      fprintf(kml_file, "<!-- Format: POLYGON (generated by %s) -->\n", 
        TOOL_SUITE_VERSION_STRING);
      write_kml_placemark(kml_file, name, center_lat, center_lon, NULL, 
        dbf, nAttr, lat, lon, nCoords, cfg);
      kml_footer(kml_file);
      FCLOSE(kml_file);
    }
    else  // must be some CSV type file
      csv2kml(inFile, outFile, format, cfg);
  }
  FREE(lat);
  FREE(lon);

  return TRUE;
}
