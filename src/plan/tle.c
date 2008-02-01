#include <stdio.h>
#include <dateUtil.h>
#include <assert.h>
#include "asf_meta.h"
#include "sgpsdp.h"

double secs_to_jul(double t)
{
  julian_date jd;
  hms_time hms;
  ymd_date ymd;

  sec2date(t, &jd, &hms);
  date_jd2ymd(&jd, &ymd);

  assert(ymd.year>1900);
  assert(ymd.day>=1 && ymd.day<=31);
  assert(ymd.month>=1&&ymd.month<=12);

  double jdoy = Julian_Date_of_Year(ymd.year);
  double doy = DOY(ymd.year,ymd.month-1,ymd.day);
  double fod = hms.sec/(60.*60.*24.) + hms.min/(60.*24.) + hms.hour/24.0;
  //double fod = (hms.hour + (hms.min + hms.sec/60.0)/60.0)/24.0;
  //double fod = Fraction_of_Day(hms.hour,hms.min,(int)(hms.sec+.5));

  double ret = jdoy + doy + fod;

  //double ret = Julian_Date(&the_time);
  //printf("%f -> %02d/%02d/%04d %02d:%02d:%02d -> %f\n", t,
  //       ymd.month, ymd.day, ymd.year, hms.hour, hms.min, (int)hms.sec, ret);

  return ret;
}

stateVector tle_propagate(sat_t *sat, double t)
{
  sat->jul_utc = secs_to_jul(t);
  sat->tsince = (sat->jul_utc - sat->jul_epoch) * xmnpda;

  SGP4(sat, sat->tsince);

  Convert_Sat_State(&sat->pos, &sat->vel);
  Magnitude(&sat->vel);
  sat->velo = sat->vel.w;

  //geodetic_t obs_geodetic;
  //obs_set_t obs_set;
  //obs_geodetic.lat = obs_geodetic.lon = obs_geodetic.alt = 0;
  //obs_geodetic.theta = 0;
  //Calculate_Obs(sat->jul_utc, &sat->pos, &sat->vel, &obs_geodetic, &obs_set);

  geodetic_t sat_geodetic;
  Calculate_LatLonAlt(sat->jul_utc, &sat->pos, &sat_geodetic);

  while (sat_geodetic.lon < -pi)
    sat_geodetic.lon += twopi;
  while (sat_geodetic.lon > pi)
    sat_geodetic.lon -= twopi;

  //sat->az = Degrees (obs_set.az);
  //sat->el = Degrees (obs_set.el);
  //sat->range = obs_set.range;
  //sat->range_rate = obs_set.range_rate;
  sat->ssplat = Degrees (sat_geodetic.lat);
  sat->ssplon = Degrees (sat_geodetic.lon);
  sat->alt = sat_geodetic.alt;
  sat->ma = Degrees (sat->phase);
  sat->ma *= 256.0/360.0;
  sat->phase = Degrees (sat->phase);

  sat->footprint = 12756.33 * acos (xkmper / (xkmper+sat->alt));
  double age = sat->jul_utc - sat->jul_epoch;
  sat->orbit = (long) floor((sat->tle.xno * xmnpda/twopi +
                             age * sat->tle.bstar * ae) * age +
                            sat->tle.xmo/twopi) + sat->tle.revnum - 1;

  stateVector st;

  st.pos.x = sat->pos.x * 1000;
  st.pos.y = sat->pos.y * 1000;
  st.pos.z = sat->pos.z * 1000;
  st.vel.x = sat->vel.x * 1000;
  st.vel.y = sat->vel.y * 1000;
  st.vel.z = sat->vel.z * 1000;

  julian_date jd;
  hms_time hms;
  sec2date(t,&jd,&hms);
  double gha = utc2gha(jd.year,jd.jd,hms.hour,hms.min,hms.sec);
  //printf("t= %f --> gha= %f\n", t, gha);
  //gha-=30.555437; if (gha<0) gha+=360;
  gei2fixed(&st,gha);

  return st;
}

void read_tle(const char *tle_filename, const char *satellite, sat_t *sat)
{
  FILE *fp = fopen(tle_filename, "r");
  char tle_str[3][80];
  if (fp) {

    while (1) {

      if (!fgets(tle_str[0], 80, fp)) break;
      if (!fgets(tle_str[1], 80, fp)) break;
      if (!fgets(tle_str[2], 80, fp)) break;

      if (strncmp(tle_str[0], satellite, strlen(satellite)) == 0) {
        // found
        int ret = Get_Next_Tle_Set (tle_str, &sat->tle);

        if (ret != 1) {
          printf("Error reading '%s' from %s\n", satellite, tle_filename);
          break;
        }
        else {
          // good satellite data obtained
          sat->flags = 0;
          select_ephemeris(sat);

          sat->jul_utc = 0.0;
          sat->tsince = 0.0;
          sat->az = 0.0;
          sat->el = 0.0;
          sat->range = 0.0;
          sat->range_rate = 0.0;
          sat->ra = 0.0;
          sat->dec = 0.0;
          sat->ssplat = 0.0;
          sat->ssplon = 0.0;
          sat->alt = 0.0;
          sat->velo = 0.0;
          sat->ma = 0.0;
          sat->footprint = 0.0;
          sat->phase = 0.0;
          sat->aos = 0.0;
          sat->los = 0.0;

          sat->jul_epoch = Julian_Date_of_Epoch(sat->tle.epoch);

          // calculate the satellite's position at t=0 (epoch time)
          SGP4 (sat, 0.0);

          // scale position and velocity to km and km/sec
          Convert_Sat_State (&sat->pos, &sat->vel);

          // get the velocity of the satellite
          Magnitude (&sat->vel);
          sat->velo = sat->vel.w;

          geodetic_t sat_geodetic;
          Calculate_LatLonAlt (sat->jul_utc, &sat->pos, &sat_geodetic);

          while (sat_geodetic.lon < -pi)
            sat_geodetic.lon += twopi;
	
          while (sat_geodetic.lon > (pi))
            sat_geodetic.lon -= twopi;

          sat->ssplat = Degrees (sat_geodetic.lat);
          sat->ssplon = Degrees (sat_geodetic.lon);
          sat->alt = sat_geodetic.alt;
          sat->ma = Degrees (sat->phase);
          sat->ma *= 256.0/360.0;
          sat->footprint = 2.0 * xkmper * acos (xkmper/sat->pos.w);
          double age = 0.0;
          sat->orbit = (long) floor((sat->tle.xno * xmnpda/twopi +
                                     age * sat->tle.bstar * ae) * age +
                                    sat->tle.xmo/twopi) + sat->tle.revnum - 1;

          // orbit type
          sat->otype = ORBIT_TYPE_UNKNOWN;

          fclose(fp);
          return;
        }
      }
    }

    // reached end of file without finding TLE
    printf("Couldn't find '%s' in %s.\n", satellite, tle_filename);
  }
  else { // failed to open TLE file
    printf("Couldn't read: %s\n", tle_filename);
  }
}
