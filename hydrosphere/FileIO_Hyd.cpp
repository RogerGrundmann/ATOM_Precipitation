
#include "cHydrosphereModel.h"

#include <cmath>

using namespace std;
using namespace AtomUtils;

void cHydrosphereModel::read_Hydrosphere_Surface_Data(int Ma){

cout << endl << "      OGCM: read_Hydrosphere_Surface_Data ......................." << endl;


    if(!is_global_temperature_curve_loaded()) 
        load_global_temperature_curve();

    if(!is_equat_temperature_curve_loaded()) 
        load_equat_temperature_curve();

    if(!is_pole_temperature_curve_loaded()) 
        load_pole_temperature_curve();


    //Prepare the temperature, precipitation and salinity, Name_Transfer_File
    string Name_Transfer_File;
    stringstream ssName_v_w_Transfer_File;
    string Name_SurfaceTemperature_File = temperature_file;
    string Name_SurfaceNASATemperature_File  = temperature_file;
    string Name_SurfacePrecipitation_File = precipitation_file;
    string Name_SurfaceSalinity_File = salinity_file;


    if(Ma != 0 && use_earthbyte_reconstruction){

        Name_SurfaceTemperature_File = output_path
            + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";
        Name_SurfacePrecipitation_File = output_path
            + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";
        Name_SurfaceSalinity_File = output_path + std::to_string(Ma)
            + "Ma_Reconstructed_Salinity.xyz";

        struct stat info;

        if(stat(Name_SurfaceSalinity_File.c_str(), &info) != 0 ||
           stat(Name_SurfacePrecipitation_File.c_str(), &info) != 0 ||
           stat(velocity_v_file.c_str(), &info) != 0 ||
           stat(velocity_w_file.c_str(), &info) != 0){

            std::string cmd_str = "python " + reconstruction_script_path 
                + " " + std::to_string(Ma - time_step) 
                + " " + std::to_string(Ma) 
                + " " + output_path 
                + " " + BathymetrySuffix +" hyd";

            int ret = system(cmd_str.c_str());

            std::cout << " reconstruction script returned: " 
                << ret << std::endl;
        }
    }

    if(!has_printed_welcome_msg)  print_welcome_msg();

    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;

    cout << endl << "      bathymetry given by the x-y-z data set:    "
        << bathymetry_name.c_str();

    // init_bathymetry must come first: it populates h.x[im-1][j][k] which
    // HydrosphereDataTransfer uses to mask Evaporation and Precipitation
    // for land cells.  Calling HydrosphereDataTransfer first left h all-zero
    // (all-water), so the masking was never applied.
    init_bathymetry(bathymetry_path + "/" + bathymetry_name);
    { auto dot = bathymetry_name.rfind('.'); if(dot != string::npos) bathymetry_name.erase(dot); }

    HydrosphereDataTransfer(bathymetry_name);


    read_IC(velocity_v_file, velocity_v_NASA.y, jm, km);
    read_IC(velocity_w_file, velocity_w_NASA.y, jm, km);    


    if((Ma != 0)&&(use_earthbyte_reconstruction))
        read_IC(Name_SurfaceTemperature_File, t.x[im-1], jm, km);       // reconstructed temperature in °C


    read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);
    read_IC(Name_SurfacePrecipitation_File, precipitation_NASA.y, jm, km);


    cout << endl << "      OGCM: read_Hydrosphere_Surface_Data ended ................." << endl;
}
/*
*
*/
void cHydrosphereModel::HydrosphereDataTransfer(const string &Name_Bathymetry_File){
    cout << endl << "      OGCM: HydrosphereDataTransfer" << endl;

    string stem_t = Name_Bathymetry_File;
    { auto dot = stem_t.rfind('.'); if(dot != string::npos) stem_t = stem_t.substr(0, dot); }
    string base_t = output_path;
    if(!base_t.empty() && base_t.back() == '/') base_t.pop_back();
    string Name_Transfer_File = base_t + "/" + stem_t + "_Transfer_Atm.vwtp";

    ifstream Transfer_File(Name_Transfer_File);

    if(!Transfer_File.is_open()){
        cout << "WARNING: transfer file not found, skipping atmosphere data transfer: "
             << Name_Transfer_File << "\n";
        cout << "      OGCM: HydrosphereDataTransfer skipped" << endl;
        return;
    }

    string line;
    bool truncated = false;
    for(int j = 0; j < jm && !truncated; j++) {
        for(int k = 0; k < km; k++) {
            if(!getline(Transfer_File, line)) {
                cerr << "WARNING: HydrosphereDataTransfer: file truncated at j="
                     << j << " k=" << k << "\n";
                truncated = true;
                break;
            }
            double vv = 0.0, wv = 0.0, tv = 0.0, pv = 0.0, ev = 0.0, pr = 0.0;
            {
                istringstream ss(line);
                ss >> vv >> wv >> tv >> pv >> ev >> pr;
            }
            v.x[im-1][j][k]       = isfinite(vv) ? vv : 0.0;            //non-dimensional
            w.x[im-1][j][k]       = isfinite(wv) ? wv : 0.0;            //non-dimensional
            t.x[im-1][j][k]       = isfinite(tv) ? tv : 0.0;            //non-dimensional
            p_dyn.x[im-1][j][k]   = isfinite(pv) ? pv : 0.0;            //non-dimensional
            Evaporation.y[j][k]   = isfinite(ev) ? ev : 0.0;            //dimensional in [mm/d]
            Precipitation_2D.y[j][k] = isfinite(pr) ? pr : 0.0;         //dimensional in [mm/d]

            if(is_land(h, im-1, j, k)) {
                Evaporation.y[j][k]   = 0.0;
                Precipitation_2D.y[j][k] = 0.0;
            }
        }
    }

    Transfer_File.close();
    cout << "      OGCM: HydrosphereDataTransfer ended" << endl;
}
/*
*
*/
void cHydrosphereModel::HydrospherePlotData(const string &Name_Bathymetry_File){
cout << endl << "      OGCM: HydrospherePlotData" << endl;
    string stem = Name_Bathymetry_File;
    { auto dot = stem.rfind('.'); if(dot != string::npos) stem = stem.substr(0, dot); }
    string base = output_path;
    if(!base.empty() && base.back() == '/') base.pop_back();
    string Name_PlotData_File = base + "/" + stem + "_PlotData_Hyd.xyz";
    ofstream PlotData_File(Name_PlotData_File);
    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    PlotData_File << "lons(deg)"
        << ", " << "lats(deg)"
        << ", " << "topography(m)"
        << ", " << "v-velocity(m/s)"
        << ", " << "w-velocity(m/s)"
        << ", " << "velocity-mag(m/s)"
        << ", " << "temperature(Celsius)"
        << ", " << "salinity(g/kg)"
        << ", " << "Ekman_pumping(cm/d)"
        << ", " << "upwelling(cm/d)"
        << ", " << "downwelling(cm/d)"  << endl;
    std::vector<std::string> buffer(km * jm);
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double vel_v   = v.x[im-1][j][k] * u_0;
            double vel_w   = w.x[im-1][j][k] * u_0;
            double vel_mag = sqrt(vel_v * vel_v + vel_w * vel_w);
            std::stringstream ss;
            ss.precision(4);
            ss << std::fixed;
            ss << k << " " << 90-j
                << " " << h.x[im-1][j][k]
                << " " << vel_v
                << " " << vel_w
                << " " << vel_mag
                << " " << t.x[im-1][j][k] * t_0 - t_0
                << " " << c.x[im-1][j][k] * c_35
                << " " << EkmanPumping.y[j][k]
                << " " << Upwelling.y[j][k]
                << " " << Downwelling.y[j][k] << " \n";
            buffer[k * jm + j] = ss.str();
        }
    }
    for(const auto& line : buffer) PlotData_File << line;
    cout << "      OGCM: HydrospherePlotData ended" << endl;
}
/*
*
*/
void cHydrosphereModel::load_global_temperature_curve(){
    load_map_from_file(temperature_global_file, m_global_temperature_curve);
}
/*
*
*/
void cHydrosphereModel::load_equat_temperature_curve(){
    load_map_from_file(temperature_equat_file, m_equat_temperature_curve);
}
/*
*
*/
void cHydrosphereModel::load_pole_temperature_curve(){
    load_map_from_file(temperature_pole_file, m_pole_temperature_curve);
}
/*
*
*/
float cHydrosphereModel::get_temperatures_from_curve(float time, 
    std::map<float, float>& m) const{
    if(time < m.begin()->first 
        || time > (--m.end())->first){
        std::cout << "Input time out of range: " << time << std::endl;    
        return NAN;
    }
    if(m.size() < 2){
        std::cout << "No enough data in map m" << std::endl;
        return NAN;
    }
    map<float, float>::const_iterator upper = m.begin(), 
        bottom = ++m.begin(); 
    for(map<float, float>::const_iterator it = m.begin();
            it != m.end(); ++it){
        if(time < it->first){
            bottom = it;
            break;
        }else{
            upper = it;
        }
    }
    return upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second);
}
/*
*
*/
void cHydrosphereModel::CircumPolarCurrent(){
      cout << endl << "      OGCM: CircumPolarCurrent" << endl;
  // south polar sea                                                                                                                                                                 
  // antarctic circumpolar current(-5000m deep)(from j=147 until j=152 compares to 57°S until 62°S,
  // from k=0 until k=km compares to 0° until 360°)                                                                                                                                  
                  
      double w_cpc = 0.56/u_0;                                            // 0.56 m/s average value

      #pragma omp parallel for collapse(3)
      for(int i = 0; i < im; i++){
          for(int j = 147; j < 153; j++){
              for(int k = 0; k < km; k++){
                  if(is_water(h, i, j, k)){
                      u.x[i][j][k] = 0.0;
                      v.x[i][j][k] = 0.0;
                      w.x[i][j][k] = w_cpc;                             // [./.]
                  }
              }
          }
      }
      cout << "      OGCM: CircumPolarCurrent ended" << endl;
  }
/*
*
*/
void cHydrosphereModel::IC_t_WestEastCoast(){
      cout << endl << "      OGCM: IC_t_WestEastCoast" << endl;
  // initial conditions for temperature at                                                                                                                                           
  // the sea surface close to west coasts
  // transition between coast flows and open sea flows included                                                                                                                      
  // upwelling along west coasts causes a temperature reduction

      int k_mid_temp = 15; //  maximun extension of temperature management from the west coasts
      double temp_red = 0.95;

      #pragma omp parallel for
      for(int j = 0; j < jm; j++){
          for(int k = 1; k < km-1; k++){
              if((is_water(h, im-1, j, k-1))&&(is_land(h, im-1, j, k))){
                  while(k >= 1){
                      if(is_water(h, im-1, j, k))  break;
                      double t_neg = temp_red * t.x[im-1][j][k];
                      int k_temp;
                      if(k <= k_mid_temp) k_temp = k_mid_temp - k;
                      else k_temp = k_mid_temp;
                      for(int l = 0; l <= k_temp; l++){
                          double temp_inter = t.x[im-1][j][k-l];
                          t.x[im-1][j][k-l] = (t.x[im-1][j][k-k_mid_temp] - t_neg)
                             /(double)k_mid_temp * (double)l + t_neg;
                          if(is_land(h, im-1, j, k-l))  t.x[im-1][j][k-l] = temp_inter;
                          if(t.x[im-1][j][k-l] <= t_pole_salt)  t.x[im-1][j][k-l] = t_pole_salt;
                      }
                      break;
                  } // end while
              } // end if
          } // end k
      } // end j

      int j_var_north = 86;
      int j_var_south = 94;
      int k_mid_eq = 100;

      #pragma omp parallel for
      for(int j = j_var_north; j <= j_var_south; j++){
          for(int k = 1; k < km-1; k++){
              if((is_water(h, im-1, j, k-1))&&(is_land(h, im-1, j, k))){
                  int k_temp = k_mid_eq;
                  double t_neg = temp_red * t.x[im-1][j][k];
                  for(int l = 0; l <= k_temp; l++){
                      double temp_inter = t.x[im-1][j][k-l];
                      t.x[im-1][j][k-l] = (t.x[im-1][j][k-k_mid_eq] - t_neg)
                         /(double)k_mid_eq * (double)l + t_neg;
                      if(is_land(h, im-1, j, k-l))  t.x[im-1][j][k-l] = temp_inter;
                      if(t.x[im-1][j][k-l] <= t_pole_salt)  t.x[im-1][j][k-l] = t_pole_salt;
                  }
              } // end if
          } // end k
      } // end j

      cout << "      OGCM: IC_t_WestEastCoast" << endl;
  }
/*
*
*/
void cHydrosphereModel::IC_u_WestEastCoast(){
      cout << endl << "      OGCM: IC_u_WestEastCoast" << endl;
  // initial conditions for the u velocity component at the sea surface close to east or west coasts                                                                                 
      int i_beg = 20;                                                     // == 100m depth
      int i_half = i_beg + 10;                                            // location of u-max                                                                                       
      int j_half =(jm-1)/2;
      int i_max = im-1;
      double d_i_half =(double)i_half;
  // search for east coasts and associated velocity components to close the circulations
  // transition between coast flows and open sea flows included
  // northern hemisphere: east coast
      int k_grad = 20;                                                    // extension of velocity change

      #pragma omp parallel for
      for(int j = 0; j <= j_half; j++){                                   // outer loop: latitude
          int k_water = 0;
          int k_sequel = 1;
          for(int k = 0; k < km; k++){                                    // inner loop: longitude
              if(is_land(h, i_half, j, k)) k_sequel = 0;                  // if solid ground: k_sequel = 0
              if((is_water(h, i_half, j, k))&&(k_sequel == 0))
                  k_water = 0;                                            // if water and and k_sequel = 0 then is water closest to coast
              else k_water = 1;                                           // somewhere on water
              if((is_water(h, i_half, j, k))&&(k_water == 0)){            // if water is closest to coast, change of velocity components begins
                  for(int l = 0; l < k_grad; l++){                        // extension of change, sign change in v-velocity and distribution of u-velocity with depth
                      if(k+l > km-1) break;
                      for(int i = i_beg; i <= i_half; i++){               // loop in radial direction, extension for u-velocity component, downwelling here
                          int m = i_max -(i - i_beg);
                          double d_i =(double)i;
                          u.x[i][j][k+l] = - d_i/d_i_half/(double)(l+1);  // increase with depth, decrease with distance from coast
                          u.x[m][j][k+l] = - d_i/d_i_half/(double)(l+1);  // decrease with depth, decrease with distance from coast
                      }
                  }
                  k_sequel = 1;                                           // looking for another east coast
              }
          }                                                               // end of longitudinal loop
      }                                                                   // end of latitudinal loop

  // southern hemisphere: east coast
      #pragma omp parallel for
      for(int j = j_half+1; j < jm; j++){
          int k_water = 0;
          int k_sequel = 1;
          for(int k = 0; k < km; k++){
              if(is_land(h, i_half, j, k)) k_sequel = 0;
              if((is_water(h, i_half, j, k))&&(k_sequel == 0))
                  k_water = 0;
              else k_water = 1;
              if((is_water(h, i_half, j, k))&&(k_water == 0)){
                  for(int l = 0; l < k_grad; l++){
                      if(k+l > km-1) break;
                      for(int i = i_beg; i <= i_half; i++){
                          int m = i_max -(i - i_beg);
                          double d_i =(double)i;
                          u.x[i][j][k+l] = - d_i/d_i_half/(double)(l+1);
                          u.x[m][j][k+l] = - d_i/d_i_half/(double)(l+1);
                      }
                  }
                  k_sequel = 1;
              }
          }
      }

  // search for east coasts and associated velocity components to close the circulations
  // transition between coast flows and open sea flows included
  // northern hemisphere: west coast
      #pragma omp parallel for
      for(int j = 0; j <= j_half; j++){                                   // outer loop: latitude
          int k_water = 0;
          int flip = 0;
          for(int k = 0; k < km; k++){                                    // inner loop: longitude
              if(is_water(h, i_half, j, k)){                              // if somewhere on water
                  k_water = 0;
                  flip = 0;
              }
              else k_water = 1;                                           // first time on land
              if((flip == 0)&&(k_water == 1)){                            // on water closest to land
                  for(int l = k; l >(k - k_grad + 1); l--){               // backward extention of velocity change: nothing changes
                      if(l < 0) break;
                      for(int i = i_beg; i <= i_half; i++){               // loop in radial direction, extension for u -velocity component, downwelling here
                          int m = i_max -(i - i_beg);
                          double d_i =(double)i;
                          u.x[i][j][l] = + d_i/d_i_half/((double)(k - l + 1));
                          u.x[m][j][l] = + d_i/d_i_half/((double)(k - l + 1));
                      }
                  }
                  flip = 1;
              }
          }
      }

  // southern hemisphere: west coast
      #pragma omp parallel for
      for(int j = j_half+1; j < jm; j++){
          int k_water = 0;
          int flip = 0;
          for(int k = 0; k < km; k++){
              if(is_water(h, i_half, j, k)){
                  k_water = 0;
                  flip = 0;
              }
              else k_water = 1;
              if((flip == 0)&&(k_water == 1)){
                  for(int l = k; l >(k - k_grad + 1); l--){
                      if(l < 0) break;
                      for(int i = i_beg; i <= i_half; i++){
                          int m = i_max -(i - i_beg);
                          double d_i =(double)i;
                          u.x[i][j][l] = + d_i/d_i_half/((double)(k - l + 1));
                          u.x[m][j][l] = + d_i/d_i_half/((double)(k - l + 1));
                      }
                  }
                  flip = 1;
              }
          }
      }

      #pragma omp parallel for collapse(3)
      for(int i = 0; i < i_beg; i++){
          for(int j = 0; j < jm; j++){
              for(int k = 0; k < km; k++){
                  u.x[i][j][k] = 0.0;
              }
          }
      }

      cout << "      OGCM: IC_u_WestEastCoast" << endl;
  }
/*
*
*/
void cHydrosphereModel::IC_Equatorial_Currents(){
      cout << endl << "      OGCM: IC_Equatorial_Currents" << endl;
                                                                                                                                                                                     
      double IC_water = 1.0;
      int i_beg = 0;                                                                                                                                                                 
      int i_EIC_u = i_beg;                                                                                                                                                         
      int i_EIC_o = 0;
      int i_SCC_u = 0;
      int i_SCC_o = 0;
      int i_ECC_u = 0;
      int i_ECC_o = im;
      int i_EUC_u = 0;
      int i_EUC_o = 20;
      int j_half = (jm - 1)/2;

      // Phase 1: collect ocean segments along the equator
      struct Segment { int k_beg; int k_end; };
      std::vector<Segment> segments;

      int k_beg = 0;
      for(int k = 1; k < km-1; k++){
          if((is_land(h, im-1, j_half, k))
              &&(is_water(h, im-1, j_half, k+1))){
              k_beg = k;
          }
          if((is_water(h, im-1, j_half, k))
              &&(is_land(h, im-1, j_half, k+1))){
              segments.push_back({k_beg, k});
          }
          if(((is_water(h, im-1, j_half, k))
              &&(is_water(h, im-1, j_half, k+1)))
              &&(k == km-2)){
              segments.push_back({k_beg, k});
          }
      }

      // Phase 2: apply currents to each segment in parallel
      for(const auto &seg : segments){

          // NECC (j=83..87, full depth)
          #pragma omp parallel for collapse(3)
          for(int i = i_ECC_u; i < i_ECC_o; i++){
              for(int j = 83; j < 88; j++){
                  for(int k = seg.k_beg; k <= seg.k_end; k++){
                      if(is_water(h, i, j, k)){
                          v.x[i][j][k] = 0.;
                          w.x[i][j][k] = + .01 * IC_water
                              *(double)(i - i_ECC_u)
                             /(double)(i_ECC_o - i_ECC_u);
                      }
                  }
              }
          }

          // SECC (j=93..97, full depth)
          #pragma omp parallel for collapse(3)
          for(int i = i_ECC_u; i < i_ECC_o; i++){
              for(int j = 93; j < 98; j++){
                  for(int k = seg.k_beg; k <= seg.k_end; k++){
                      if(is_water(h, i, j, k)){
                          v.x[i][j][k] = 0.;
                          w.x[i][j][k] = + .01 * IC_water
                              *(double)(i - i_ECC_u)
                             /(double)(i_ECC_o - i_ECC_u);
                      }
                  }
              }
          }

          // EUC - Cromwell current (j=87..93)
          #pragma omp parallel for collapse(3)
          for(int i = i_EUC_u; i < i_EUC_o; i++){
              for(int j = 87; j < 94; j++){
                  for(int k = seg.k_beg; k <= seg.k_end; k++){
                      if(is_water(h, i, j, k)){
                          v.x[i][j][k] = 0.;
                          w.x[i][j][k] = + .4 * IC_water;
                      }
                  }
              }
          }

          // EIC (j=88..92)
          #pragma omp parallel for collapse(3)
          for(int i = i_EIC_u; i < i_EIC_o; i++){
              for(int j = 88; j < 93; j++){
                  for(int k = seg.k_beg; k <= seg.k_end; k++){
                      if(is_water(h, i, j, k)){
                          v.x[i][j][k] = 0.;
                          w.x[i][j][k] = - .05 * IC_water
                              *(double)(i - i_EIC_u)
                             /(double)(i_EIC_o - i_EIC_u);
                      }
                  }
              }
          }

          // NSCC (j=86..87)
          #pragma omp parallel for collapse(3)
          for(int i = i_SCC_u; i < i_SCC_o; i++){
              for(int j = 86; j < 88; j++){
                  for(int k = seg.k_beg; k <= seg.k_end; k++){
                      if(is_water(h, i, j, k)){
                          v.x[i][j][k] = 0.;
                          w.x[i][j][k] = .025 * IC_water
                              *(double)(i - i_SCC_u)
                             /(double)(i_SCC_o - i_SCC_u);
                      }
                  }
              }
          }

          // SSCC (j=93..94)
          #pragma omp parallel for collapse(3)
          for(int i = i_SCC_u; i < i_SCC_o; i++){
              for(int j = 93; j < 95; j++){
                  for(int k = seg.k_beg; k <= seg.k_end; k++){
                      if(is_water(h, i, j, k)){
                          v.x[i][j][k] = 0.;
                          w.x[i][j][k] = .025 * IC_water
                              *(double)(i - i_SCC_u)
                             /(double)(i_SCC_o - i_SCC_u);
                      }
                  }
              }
          }
      }

      cout << "      OGCM: IC_Equatorial_Currents ended" << endl;
  }
/*
*
*/
void cHydrosphereModel::init_bathymetry(const string &bathymetry_file){

      cout << endl << "      OGCM: init_bathymetry" << endl;                                                                                                                         
   
      ifstream ifile(bathymetry_file);                                                                                                                                               
                  
      if (!ifile.is_open()){
          cerr << "ERROR: could not open bathymetry file " << bathymetry_file << "\n";
          abort();
      }

      int j, k;
      double lon, lat, depth;

      // file I/O must remain serial
      for(j = 0; j < jm && !ifile.eof(); j++){
          for(k = 0; k < km && !ifile.eof(); k++){

              depth = 999.0;
              ifile >> lon >> lat >> depth;

              if( depth > 0.0){
                  depth = 0.0;
              }

              // Invert stretched coordinate: depth (≤0) → grid index i_top.
              // Uniform fallback (dr_stretch≈0): ξ = 1 + depth/L_hyd  (original formula).
              // Stretched (dr_stretch > 0):       ξ = 1 − asinh(|depth|/L_hyd · sinh(β)) / β
              //
              // Note: the atmosphere model (init_topography) does NOT need this explicit
              // inversion because it calls get_layer_height(i), which reads m_layer_heights[]
              // computed directly from rad.z[i] in init_layer_heights().  The ocean model has
              // no equivalent layer-height table and originally computed i_top arithmetically
              // from uniform spacing, so the inversion must be done explicitly here.
              double xi;
              if (dr_stretch < 1.0e-10) {
                  xi = 1.0 + depth / L_hyd;
              } else {
                  xi = 1.0 - std::asinh(-depth / L_hyd * std::sinh(dr_stretch)) / dr_stretch;
              }
              int i_top_raw = static_cast<int>(std::round((im - 1) * xi));
              int i_top = std::max(0, std::min(im - 1, i_top_raw));

              i_bathymetry[j][k] = i_top;

              Bathymetry.y[j][k] = - depth;

              // When xi < 0 the ocean is deeper than the model domain: i_top was
              // clamped up from a negative value, so i=0 is not a real seafloor.
              // Start from i=1 in that case so no artificial bottom line is drawn.
              int i_floor = (i_top_raw >= 0) ? 0 : 1;
              for(int i = i_floor; i <= i_top; i++){
                  h.x[i][j][k] = 1.0;
              }

              if(ifile.fail()){
                  ifile.clear();
                  std::string tmp;
                  std::getline(ifile, tmp);
                  logger() << "bad data in topography at: " << lon << " " << lat << " " << tmp << std::endl;
              }

          }
      }


  //  reduction and smoothing of peaks and needles in bathymetry
  //  i=0 is excluded: the bottom row is never a spike and must not be opened
      #pragma omp parallel for collapse(3)
      for(int i = 1; i < im; i++){
          for(int j = 1; j < jm-1; j++){
              for(int k = 1; k < km-1; k++){

                  if((is_land(h, i, j, k))
                      &&((is_water(h, i, j-1, k))
                      &&(is_water(h, i, j+1, k))))
                      h.x[i][j][k] = 0.0;

                  if((is_land(h, i, j, k))
                      &&((is_water(h, i, j, k-1))
                      &&(is_water(h, i, j, k+1))))
                      h.x[i][j][k] = 0.0;
              }
          }
      }

  //  rewrite bathymetric data from -180° - 0° - +180° to 0°- 360°
      #pragma omp parallel for
      for(int j = 0; j < jm; j++){

          move_data(Bathymetry.y[j], km);

          for(int i = 0; i < im; i++){
              move_data(h.x[i][j], km);
              move_data(i_bathymetry[j], km);
          }
      }


      cout << "      OGCM: init_bathymetry ended" << endl;
  }
/*
*
*/
// Binary dump of the prognostic ocean fields so a spin-up can resume at a chosen
// iter instead of re-running from the IC.  Mirrors the atmosphere
// (cAtmosphereModel::save_state/load_state, FileIO_Atm.cpp).  Only the genuinely
// prognostic arrays are stored; densities, forces and all diagnostics are
// recomputed each iteration from these.
std::vector<Array*> cHydrosphereModel::restart_arrays(){
    return { &t, &u, &v, &w, &c,
             &tn, &un, &vn, &wn, &cn,
             &p_dyn, &p_hydro,
             &tke, &dis, &tken, &disn, &nue };
}

void cHydrosphereModel::save_state(int iter){
    const string fn = output_path + "/hyd_restart_" + std::to_string(iter) + ".bin";
    std::ofstream f(fn, std::ios::binary);
    if(!f){
        cout << "      OGCM: save_state FAILED to open " << fn << endl;
        return;
    }
    // Header: magic, dimensions, and the iter to resume at.
    const int32_t hdr[5] = { 0x4F434D31 /*"OCM1"*/, im, jm, km, total_iter_count };
    f.write(reinterpret_cast<const char*>(hdr), sizeof(hdr));

    std::vector<Array*> arrs = restart_arrays();
    for(Array* a : arrs)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                f.write(reinterpret_cast<const char*>(a->x[i][j]), km * sizeof(double));

    cout << "      OGCM: save_state wrote " << arrs.size() << " arrays to "
         << fn << " (total_iter_count " << total_iter_count << ")" << endl;
}

bool cHydrosphereModel::load_state(int iter){
    const string fn = output_path + "/hyd_restart_" + std::to_string(iter) + ".bin";
    std::ifstream f(fn, std::ios::binary);
    if(!f){
        cout << "      OGCM: load_state: no file " << fn
             << " — running from scratch" << endl;
        return false;
    }
    int32_t hdr[5];
    f.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
    if(!f || hdr[0] != 0x4F434D31 || hdr[1] != im || hdr[2] != jm || hdr[3] != km){
        cout << "      OGCM: load_state: bad header / grid mismatch in " << fn
             << " — running from scratch" << endl;
        return false;
    }

    std::vector<Array*> arrs = restart_arrays();
    for(Array* a : arrs)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++){
                f.read(reinterpret_cast<char*>(a->x[i][j]), km * sizeof(double));
                if(!f){
                    cout << "      OGCM: load_state: truncated file " << fn
                         << " — running from scratch" << endl;
                    return false;
                }
            }

    total_iter_count = hdr[4];
    cout << "      OGCM: load_state restored " << arrs.size() << " arrays from "
         << fn << " (resuming at total_iter_count " << total_iter_count << ")" << endl;
    return true;
}
/*
*
*/
