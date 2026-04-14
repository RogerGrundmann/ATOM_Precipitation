#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;
 

void cAtmosphereModel::read_Atmosphere_Surface_Data(int Ma){

cout << endl << endl << endl << "      AGCM: read_Atmosphere_Surface_Data ......................." << endl;

    auto begin = std::chrono::high_resolution_clock::now();


    if(!has_printed_welcome_msg)  print_welcome_msg();


// loading topography
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;

    std::cout << endl << "      topography given by the x-y-z data set:    " 
        << bathymetry_name.c_str() << endl;

    init_topography(bathymetry_path + "/" + bathymetry_name);
    { auto dot = bathymetry_name.rfind('.'); if(dot != string::npos) bathymetry_name.erase(dot); }




// temperature file loading for around 500 Mio years
    #pragma omp parallel sections
    {
        #pragma omp section
        { if(!is_global_temperature_curve_loaded()) load_global_temperature_curve(); }// loading of global temperature files from Scotese et al., ESR 2021 in °C
        #pragma omp section
        { if(!is_equat_temperature_curve_loaded()) load_equat_temperature_curve(); }  // loading of equatorial temperature files from Scotese et al., ESR 2021 in °C
        #pragma omp section
        { if(!is_pole_temperature_curve_loaded()) load_pole_temperature_curve(); }    // loading of polar temperature files from Scotese et al., ESR 2021 in °C
    }


// preparing temperature, precipitation and velocity data file constructed by reconstruct_atom_data.py
    string Name_SurfaceTemperature_File  = temperature_file;            // if  Ma > 0, temperature file name
    string Name_SurfaceNASATemperature_File  = temperature_file;        // if  Ma = 0, NASA temperature file name
    string Name_SurfaceNASAPrecipitation_File = precipitation_file;     // precipitation file name

    if((*get_current_time() > 0)&&(use_earthbyte_reconstruction)){      // reading temperature, precipitation, v and w in case Ma > 0
        Name_SurfaceTemperature_File = output_path
            + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";  // reconstructed temperature in °C

        Name_SurfaceNASAPrecipitation_File = output_path
            + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";// reconstructed precipitation in mm/d

        velocity_v_file = output_path + std::to_string(Ma)
            + "Ma_Reconstructed_wind_v.xyz";                            // reconstructed v-velocity in m/s

        velocity_w_file = output_path + std::to_string(Ma)
            + "Ma_Reconstructed_wind_w.xyz";                            // reconstructed w-velocity in m/s


        struct stat info;                                               // initiates the search of corrupted reconstruction files

        if(stat(output_path.c_str(), &info) != 0){                      // if output path does not exist mkdir is applied
             mkdir(output_path.c_str(), 0777);
        }


        if(stat(Name_SurfaceTemperature_File.c_str(), &info) != 0 ||    // if any file is corrupt an error message is posted,  0 == executed successfully
           stat(Name_SurfaceNASAPrecipitation_File.c_str(), &info) != 0 ||
           stat(velocity_v_file.c_str(), &info) != 0 ||
           stat(velocity_w_file.c_str(), &info) != 0){

            std::string cmd_str = "python " + reconstruction_script_path// python ../reconstruction/reconstruct_atom_data.py, runs the python code
                + " " + std::to_string(Ma - time_step)                  // sys.argv[1], preceding timestep
                + " " + std::to_string(Ma)                              // sys.argv[2], recent timestep
                + " " + output_path                                     // sys.argv[3], output path where to store reconstructed files
                + " " + BathymetrySuffix                                // sys.argv[4], BathymetrySuffix here 'Ma_smooth.xyz'
                + " atm";                                               // sys.argv[5], here Atmosphere code (atm or hyd) applied

            int ret = system(cmd_str.c_str());                          // executes an external operating system command stored in std::string cmd_str, here reconstruct_atom_data.py

            std::cout << " reconstruction script returned: "            // return value ret = 0 means successful ended, ret = -1 means the command processor for python cannot be created
                << ret << std::endl;
        } 
    }


    #pragma omp parallel sections
    {
        #pragma omp section
        {
            if((*get_current_time() > 0) && (use_earthbyte_reconstruction)){
                read_IC(Name_SurfaceTemperature_File, t.x[0], jm, km);
            }
        }
        #pragma omp section
        {
            read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);
            read_IC(Name_SurfaceNASAPrecipitation_File, precipitation_NASA.y, jm, km);

//            read_IC(velocity_v_file, velocity_v_NASA.y, jm, km);                // reconstructed v-velocity in m/s

//            read_IC(velocity_w_file, velocity_w_NASA.y, jm, km);                // reconstructed v-velocity in m/s
        }
    }


    cout << "      AGCM: read_Atmosphere_Surface_Data ended ................." << endl << endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for read_Atmosphere_Surface_Data\n", elapsed.count() * 1e-9);

}
/*
*
*/
void cAtmosphereModel::AtmosphereDataTransfer(const string &Name_Bathymetry_File){
    cout << endl << "      AGCM: AtmosphereDataTransfer" << endl;

    string stem_t = Name_Bathymetry_File;
    { auto dot = stem_t.rfind('.'); if(dot != string::npos) stem_t = stem_t.substr(0, dot); }
    string base_t = output_path;
    if(!base_t.empty() && base_t.back() == '/') base_t.pop_back();
    string Name_Transfer_File = base_t + "/" + stem_t + "_Transfer_Atm.vwtp";
    ofstream Transfer_File(Name_Transfer_File);

    if(!Transfer_File.is_open()){
        cerr << "ERROR: could not open transfer file: " << Name_Transfer_File << "\n";
        abort();
    }

    // prepare buffer for parallel line formatting
    std::vector<std::string> line_buffer(jm * km);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            std::stringstream ss;
            ss.precision(4);
            ss.setf(ios::fixed);

            if(is_land(h, 0, j, k)){
                ss << "0.0000 0.0000 0.0000 0.0000 0.0000 0.0000";
            }else{
                ss << v.x[0][j][k] << " "                               //non-dimensional
                   << w.x[0][j][k] << " "                               //non-dimensional
                   << t.x[0][j][k] << " "                               //non-dimensional
                   << p_dyn.x[0][j][k] << " "                           //non-dimensional
                   << Evaporation.y[j][k] << " "                        //dimensional in [mm/d]
                   << Precipitation.x[0][j][k] * 86400.0;               //dimensional in [mm/d]
            }

            // write to buffer (index must be unique per thread)
            line_buffer[j * km + k] = ss.str();
        }
    }

    // serial write of finished buffer to file
    for(int i = 0; i < jm * km; i++){
        Transfer_File << line_buffer[i] << "\n";
    }

    Transfer_File.close();
    cout << "      AGCM: AtmosphereDataTransfer ended" << endl;
}
/*
*
*/
void cAtmosphereModel::AtmospherePlotData(const string &Name_Bathymetry_File){
    cout << endl << "      AGCM: AtmospherePlotData" << endl;

    string stem = Name_Bathymetry_File;
    { auto dot = stem.rfind('.'); if(dot != string::npos) stem = stem.substr(0, dot); }
    string base = output_path;
    if(!base.empty() && base.back() == '/') base.pop_back();
    string Name_PlotData_File = base + "/" + stem + "_PlotData_Atm.xyz";
    ofstream PlotData_File(Name_PlotData_File);

    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);

    // Header schreiben
    PlotData_File << "lons(deg), lats(deg), topography(m), v-velocity(m/s), w-velocity(m/s), "
                  << "velocity-mag(m/s), temperature(Celsius), water_vapour(g/kg), "
                  << "precipitation(mm/d), precipitable water(mm), pressure-static (hPa), "
                  << "temp_landscape (Celsius), p_stat_landscape (hPa)" << endl;

    // Puffer für alle Zeilen (km * jm)
    std::vector<std::string> buffer(km * jm);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            // Berechnungen parallel durchführen
            double vel_v = v.x[0][j][k] * u_0;
            double vel_w = w.x[0][j][k] * u_0;
            double vel_mag = sqrt(vel_v * vel_v + vel_w * vel_w); // schneller als pow()

            std::stringstream ss;
            ss.precision(4);
            ss.setf(ios::fixed);

            ss << k << " " << 90-j 
               << " " << h.x[0][j][k] 
               << " " << vel_v
               << " " << vel_w
               << " " << vel_mag
               << " " << t.x[0][j][k] * t_0 - t_0
               << " " << c.x[0][j][k] * 1000.0
               << " " << Precipitation.x[0][j][k] * 86400.0
               << " " << precipitable_water.y[j][k] 
               << " " << p_stat.x[0][j][k]
               << " " << temp_landscape.y[j][k]
               << " " << p_stat_landscape.y[j][k];

            buffer[k * jm + j] = ss.str();
        }
    }

    // Serielles Schreiben des Puffers
    for(const auto& line : buffer) {
        PlotData_File << line << "\n";
    }

    PlotData_File.close();
    cout << "      AGCM: AtmospherePlotData ended" << endl;
}
/*
*
*/
void cAtmosphereModel::init_topography(const string &topo_filename){

    cout << endl << endl << endl << "      AGCM: init_topography" << endl;

    ifstream ifile(topo_filename);

    if(! ifile.is_open()){
        std::cerr << "ERROR: could not open Name_Bathymetry_File file: " 
            <<  topo_filename << std::endl;
        abort();
    }

    double lon, lat, height;

    for(int j = 0; j < jm && !ifile.eof(); j++){
        for(int k = 0; k < km && !ifile.eof(); k++){

            height = -999.0;                                            // in case the height is NaN

            ifile >> lon >> lat >> height;

            if(!(height > 0.0)){
                h.x[0][j][k] = Topography.y[j][k] = 0.0;
           }else{
                Topography.y[j][k] = height;

                for(int i = 0; i < im; i++){

                    if(height > get_layer_height(i)){
                        h.x[i][j][k] = 1.0;
                    }else{
                        i_topography[j][k] = i-1;
                        break;
                    }
                }
            }

            if(ifile.fail()){
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " 
                    << lat << " " << tmp << std::endl;
            }   

        }   
    }
//  reduction and smoothing of peaks and needles in the topography
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < im; i++){
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){

                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j-1, k))
                     &&(is_air(h, i, j+1, k)))){
                    h.x[i][j][k] = 0.0;
                }
                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j, k-1))
                     &&(is_air(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.0;
                }

            }
        }
    }


    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im; i++){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                if((is_air(h, i, j, k))&&(is_land(h, i-1, j, k))){
                    i_topography[j][k] = i;
                    i_landscape[j][k] = get_layer_height(i);
                    Landscape.y[j][k] = get_layer_height(i);
                }
            }
        }
    }



// rewriting bathymetrical data from -180° _ 0° _ +180° coordinate system to 0°- 360°
    #pragma omp parallel for schedule(dynamic)
    for(int j = 0; j < jm; j++){
        move_data(Topography.y[j], km);
        move_data(i_topography[j], km);
        move_data(i_landscape[j], km);
        move_data(Landscape.y[j], km);

        for(int i = 0; i < im; i++){
            move_data(h.x[i][j], km);
        }
    }


    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){

            int i_mount = i_topography[j][k];

            for(int i = 0; i < im; i++){

                if(i == i_mount)
                    Landscape.y[j][k] = get_layer_height(i);
                break;
            }
        }
    }

    cout << "      AGCM: init_topography ended" << endl;
}
