
#include "cHydrosphereModel.h"

#include <cmath>
#include <dirent.h>
#include <vector>
#include <sstream>

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
    // Select which atmosphere surface-coupling snapshot to read. The atmosphere writes an
    // iter-stamped <stem>_Transfer_Atm_<iter>.vwtp every checkpoint (see AtmosphereDataTransfer).
    //   atm_transfer_iter >= 0 : read exactly that iteration.
    //   atm_transfer_iter <  0 : (default) read the LATEST — the highest-iter snapshot present in
    //       the output dir, i.e. wherever the atmosphere run happened to stop. Found by scanning
    //       the directory (no fixed-name "latest" file is written by the atmosphere any more).
    const string prefix = stem_t + "_Transfer_Atm_";
    string Name_Transfer_File;
    if(atm_transfer_iter >= 0){
        Name_Transfer_File = base_t + "/" + prefix + std::to_string(atm_transfer_iter) + ".vwtp";
    }else{
        int best = -1;
        if(DIR *dir = opendir(base_t.c_str())){
            while(struct dirent *ent = readdir(dir)){
                const string fn = ent->d_name;
                // match <prefix><digits>.vwtp exactly
                if(fn.size() > prefix.size() + 5
                && fn.compare(0, prefix.size(), prefix) == 0
                && fn.compare(fn.size() - 5, 5, ".vwtp") == 0){
                    const string mid = fn.substr(prefix.size(), fn.size() - prefix.size() - 5);
                    if(!mid.empty() && mid.find_first_not_of("0123456789") == string::npos){
                        const int n = std::stoi(mid);
                        if(n > best) best = n;
                    }
                }
            }
            closedir(dir);
        }
        if(best < 0){
            cout << "WARNING: no iter-stamped transfer file (" << prefix
                 << "<iter>.vwtp) found in " << base_t
                 << ", skipping atmosphere data transfer\n";
            cout << "      OGCM: HydrosphereDataTransfer skipped" << endl;
            return;
        }
        Name_Transfer_File = base_t + "/" + prefix + std::to_string(best) + ".vwtp";
    }
    cout << "      OGCM: atmosphere transfer file = " << Name_Transfer_File
         << (atm_transfer_iter >= 0 ? "  (atm_transfer_iter)" : "  (latest snapshot)") << endl;

    ifstream Transfer_File(Name_Transfer_File);

    if(!Transfer_File.is_open()){
        cout << "WARNING: transfer file not found, skipping atmosphere data transfer: "
             << Name_Transfer_File << "\n";
        cout << "      OGCM: HydrosphereDataTransfer skipped" << endl;
        return;
    }

    string line;

    // The atmosphere writes a one-line headline first (e.g. "# iter_n = 155").
    // Consume it so the jm*km data rows line up; tolerate older header-less
    // transfer files by rewinding to the start when no '#' headline is present.
    if(getline(Transfer_File, line)){
        if(!line.empty() && line[0] == '#'){
            cout << "      OGCM: transfer file headline: " << line << endl;
        }else{
            Transfer_File.clear();
            Transfer_File.seekg(0, ios::beg);
        }
    }

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
// Reverse coupling channel (hydrosphere -> atmosphere). Writes the ocean surface
// SST (t.x[im-1], non-dimensional) so a subsequent atmosphere run can blend it into
// its ocean surface temperature (read_Hydrosphere_SST / sst_coupling_alpha), closing
// the atm<->hyd Picard loop. Mirrors AtmosphereDataTransfer:
//   - iter-stamped file <stem>_Transfer_Hyd_SST_<total_iter_count>.vwtp, never clobbered
//     (total_iter_count, like the VTK stamp, so a restart cannot overwrite an earlier
//     round's snapshot; the atmosphere reader picks the highest-iter file = latest);
//   - one SST value per (j,k) row, land cells written 0, a "# iter_n = ..." headline;
//   - NaN/Inf guard: if any surface value is non-finite the file is NOT written, so the
//     last clean snapshot stays on disk (a corrupt SST must never seed the atmosphere).
// This is a pure OUTPUT of the current ocean state; it does not touch any field.
void cHydrosphereModel::HydrosphereSSTTransfer(const string &Name_Bathymetry_File){
    cout << endl << "      OGCM: HydrosphereSSTTransfer   total_iter_count = "
         << total_iter_count << endl;

    string stem_t = Name_Bathymetry_File;
    { auto dot = stem_t.rfind('.'); if(dot != string::npos) stem_t = stem_t.substr(0, dot); }
    string base_t = output_path;
    if(!base_t.empty() && base_t.back() == '/') base_t.pop_back();
    const string stamped_name = base_t + "/" + stem_t + "_Transfer_Hyd_SST_"
                              + std::to_string(total_iter_count) + ".vwtp";

    // Build the buffer first, checking every value for NaN/Inf before touching disk.
    std::vector<std::string> line_buffer(jm * km);
    bool has_nonfinite = false;

    #pragma omp parallel for collapse(2) schedule(static) reduction(||: has_nonfinite)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            std::stringstream ss;
            ss.precision(6);
            ss.setf(ios::fixed);
            if(is_land(h, im-1, j, k)){
                ss << "0.000000";
            }else{
                double sst = t.x[im-1][j][k];                          // non-dimensional
                if(!is_finite_safe(sst)) has_nonfinite = true;
                ss << sst;
            }
            line_buffer[j * km + k] = ss.str();
        }
    }

    if(has_nonfinite){
        cout << "      OGCM: HydrosphereSSTTransfer SKIPPED — non-finite surface SST present;"
             << " previous snapshot left intact" << endl;
        return;
    }

    ofstream SST_File(stamped_name);
    if(!SST_File.is_open()){
        cerr << "ERROR: could not open hydrosphere SST transfer file: " << stamped_name << "\n";
        abort();
    }
    SST_File << "# iter_n = " << total_iter_count << "\n";
    for(int i = 0; i < jm * km; i++)
        SST_File << line_buffer[i] << "\n";
    SST_File.close();
    cout << "      OGCM: HydrosphereSSTTransfer ended (wrote " << stamped_name << ")" << endl;
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
    // THE SIZE TEST MUST COME FIRST. It used to sit BELOW the range test, which
    // dereferences m.begin() and decrements m.end() -- both undefined behaviour on an
    // empty map, and (--m.end()) is UB whether or not the map is empty when begin()==end().
    // The guard that was written to catch a too-small map could not run until after the
    // code it was guarding. Found in ATHAD, where the curve machinery was deleted outright
    // (one epoch, no time slices); here the curves are real, so the fix is the order.
    if(m.size() < 2){
        std::cout << "No enough data in map m" << std::endl;
        return NAN;
    }
    if(time < m.begin()->first 
        || time > (--m.end())->first){
        std::cout << "Input time out of range: " << time << std::endl;    
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

void cHydrosphereModel::save_state(int iter, int Ma){
    const string fn = output_path + "/hyd_restart_" + std::to_string(Ma) + "Ma_"
                      + std::to_string(iter) + ".bin";
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

bool cHydrosphereModel::load_state(int iter, int Ma){
    const string fn = output_path + "/hyd_restart_" + std::to_string(Ma) + "Ma_"
                      + std::to_string(iter) + ".bin";
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
