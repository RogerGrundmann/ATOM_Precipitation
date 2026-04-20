/*
 * Cube General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#include "cCubeModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

namespace ParaViewAtm{
    inline double safe_val(double v){ return std::isfinite(v) ? v : 0.0; }

    void dump_array(const string &name, Array &a, double multiplier, ofstream &f){
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(int k = 0; k < a.km; k++){
            for(int j = 0; j < a.jm; j++){
                for(int i = 0; i < a.im; i++){
                    f << safe_val(a.x[i][j][k] * multiplier) << "\n";
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
        f.clear();
    }
/*
 * 
*/
    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << "\n";
        f << "LOOKUP_TABLE default" << "\n";
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.x[i][j][k] * multiplier) << "\n";
            }
        }
    }
/*
 * 
*/
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << "\n";
        f << "LOOKUP_TABLE default" << "\n";
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.y[j][k] * multiplier) << "\n";
            }
        }
    }
/*
 * 
*/
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << "\n";
        f <<  "LOOKUP_TABLE default" << "\n";
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << safe_val(a.x[i][j][k] * multiplier) << "\n";
            }
        }
    }
/*
 * 
*/
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << "\n";
        f << "LOOKUP_TABLE default" << "\n";
        for(int i = 0; i < a.im; i++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.x[i][j][k] * multiplier) << "\n";
            }
        }
    }
}
/*
 * 
*/
void cCubeModel::paraview_panorama_vts(string &Name_Bathymetry_File, int n){

    using namespace ParaViewAtm;

    string Cube_panorama_vts_File_Name = output_path + "/"
        + Name_Bathymetry_File + "Cube_turbulent_panorama_" + std::to_string(n) + ".vts";
    ofstream Cube_panorama_vts_File;

    Cube_panorama_vts_File.precision(4);
    Cube_panorama_vts_File.setf(ios::fixed);

    Cube_panorama_vts_File.open(Cube_panorama_vts_File_Name);

    if(!Cube_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open panorama_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    ostringstream buf;
    buf.precision(4);
    buf.setf(ios::fixed);

    buf << "<?xml version=\"1.0\"?>\n\n"
        << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\n"
        << " <StructuredGrid WholeExtent=\"" << 1 << " " << im << " " << 1 << " " << jm << " " << 1 << " " << km << "\">\n\n"
        << "  <Piece Extent=\"" << 1 << " " << im << " " << 1 << " " << jm << " " << 1 << " " << km << "\">\n\n"
        << "   <PointData Vectors=\"Velocity MagneticField\" Scalars=\"Topography u-component v-component w-component Temperature CondensationTemp EvaporationTemp Epsilon_3D PressureDynamic PressureStatic WaterVapour CloudWater CloudIce CO2-Concentration Q_Latent Rain RainSuper Ice PrecipitationRain PrecipitationSnow PrecipitationConv Updraft Downdraft\">\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n\n";

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                buf << u.x[i][j][k] << " "
                    << v.x[i][j][k] << " " << w.x[i][j][k] << '\n';
            }
            buf << "\n\n";
        }
        buf << "\n\n";
    }

    buf << "\n\n    </DataArray>\n\n";
    Cube_panorama_vts_File << buf.str();
    buf.str("");
    buf.clear();



    dump_array("Topography", h, 1.0, Cube_panorama_vts_File);

    dump_array("u-component", u, u_0, Cube_panorama_vts_File);
    dump_array("v-component", v, u_0, Cube_panorama_vts_File);
    dump_array("w-component", w, u_0, Cube_panorama_vts_File);

    dump_array("p_dyn", p_dyn, p_0, Cube_panorama_vts_File);


    if(turb_model != "none"){
        dump_array("TKE",           tke,        1e3, Cube_panorama_vts_File);
        dump_array("Dissipation",   dis,        1.0, Cube_panorama_vts_File);
        dump_array("EddyViscosity", nue,        1e3, Cube_panorama_vts_File);
        dump_array("TurbProd",      prod,       1e3, Cube_panorama_vts_File);
        dump_array("TKE_Source",    tke_source, 1e3, Cube_panorama_vts_File);
        dump_array("Dis_Source",    dis_source, 1.0, Cube_panorama_vts_File);
    }

    buf << "   </PointData>\n\n"
        << "   <Points>\n\n"
        << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n\n";

    double dx = 0.1;
    double dy = 0.1;
    double dz = 0.1;

    for(int k = 0; k < km; k++){
        double z = k * dz;
        for(int j = 0; j < jm; j++){
            double y = j * dy;
            for(int i = 0; i < im; i++){
                buf << i * dx << " " << y << " " << z << '\n';
            }
            buf << "\n\n";
        }
        buf << "\n\n";
    }

    buf << "    </DataArray>\n\n"
        << "   </Points>\n\n"
        << "  </Piece>\n\n"
        << " </StructuredGrid>\n\n"
        << "</VTKFile>\n\n";
    Cube_panorama_vts_File << buf.str();
    Cube_panorama_vts_File.close();

    cout << "   File:  " << Cube_panorama_vts_File_Name
        << "  has been written\n";
}
/*
*
*/
void cCubeModel::paraview_vtk_radial(string &Name_Bathymetry_File, 
    int i_radial, int n){
    using namespace ParaViewAtm;
    string Cube_radial_File_Name = output_path + "/"
        + Name_Bathymetry_File + "Cube_turbulent_radial_" + std::to_string(i_radial)
        + "_" + std::to_string(n) + ".vtk";
    ofstream Cube_vtk_radial_File;
    Cube_vtk_radial_File.precision(4);
    Cube_vtk_radial_File.setf(ios::fixed);
    Cube_vtk_radial_File.open(Cube_radial_File_Name);
    if(!Cube_vtk_radial_File.is_open()){
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    ostringstream buf;
    buf.precision(4);
    buf.setf(ios::fixed);

    buf << "# vtk DataFile Version 3.0\n"
        << "Radial_Data_Cube_Circulation\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_GRID\n"
        << "DIMENSIONS " << km << " " << jm << " " << 1 << '\n'
        << "POINTS " << jm * km << " float\n";

    double dx = 0.1;
    double dy = 0.1;

    for(int j = 0; j < jm; j++){
        double x = j * dx;
        for(int k = 0; k < km; k++){
            buf << x << " " << k * dy << " " << 0.0 << '\n';
        }
    }

    buf << "POINT_DATA " << jm * km << '\n';
    Cube_vtk_radial_File << buf.str();
    Cube_vtk_radial_File.flush();
    if(!Cube_vtk_radial_File){
        cerr << "ERROR: write failed for " << Cube_radial_File_Name
             << " after header (" << jm * km << " points, ~"
             << (buf.str().size() / 1024) << " KB). Stream state: "
             << "fail=" << Cube_vtk_radial_File.fail()
             << " bad=" << Cube_vtk_radial_File.bad() << "\n";
        abort();
    }
    buf.str("");
    buf.clear();
    dump_radial("u-Component", u, u_0, i_radial, Cube_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Cube_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Cube_vtk_radial_File);

    dump_radial("Topography", h, 1.0, i_radial, Cube_vtk_radial_File);

    dump_radial("p_dyn", p_dyn, p_0, i_radial, Cube_vtk_radial_File);


    if(turb_model != "none"){
        dump_radial("TKE",        tke,         1e3, i_radial, Cube_vtk_radial_File);
        dump_radial("Dissipation", dis,        1.0, i_radial, Cube_vtk_radial_File);
        dump_radial("EddyViscosity", nue,      1e3, i_radial, Cube_vtk_radial_File);
        dump_radial("TurbProd",   prod,        1e3, i_radial, Cube_vtk_radial_File);
        dump_radial("TKE_Source", tke_source,  1e3, i_radial, Cube_vtk_radial_File);
        dump_radial("Dis_Source", dis_source,  1.0, i_radial, Cube_vtk_radial_File);
    }


    buf << "VECTORS v-w-ATOM float\n";
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            buf << v.x[i_radial][j][k]
                << " " << w.x[i_radial][j][k] << " " << 0.0 << '\n';
        }
    }

    Cube_vtk_radial_File << buf.str();
    buf.str("");
    buf.clear();


    Cube_vtk_radial_File << buf.str();


    Cube_vtk_radial_File.close();
    cout << "   File:  " << Cube_radial_File_Name
        << "  has been written\n";
}
/*
 * 
*/
void cCubeModel::paraview_vtk_zonal(string &Name_Bathymetry_File, 
    int k_zonal, int n){
    using namespace ParaViewAtm;
    string Cube_zonal_File_Name = output_path + "/" + Name_Bathymetry_File
        + "Cube_turbulent_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    ofstream Cube_vtk_zonal_File;
    Cube_vtk_zonal_File.precision(4);
    Cube_vtk_zonal_File.setf(ios::fixed);
    Cube_vtk_zonal_File.open(Cube_zonal_File_Name);
    if(!Cube_vtk_zonal_File.is_open()){
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    // Buffer all output to reduce I/O syscalls
    ostringstream buf;
    buf.precision(4);
    buf.setf(ios::fixed);

    buf << "# vtk DataFile Version 3.0\n"
        << "Zonal_Data_Cube_Circulation\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_GRID\n"
        << "DIMENSIONS " << jm << " " << im << " " << 1 << '\n'
        << "POINTS " << im * jm << " float\n";

//    double dx = 0.1;
    double dx = 0.05;
    double dy = 0.05;

    for(int i = 0; i < im; i++){
        double x = i * dx;
        for(int j = 0; j < jm; j++){
            buf << x << " " << j * dy << " " << 0.0 << '\n';
        }
    }

    buf << "POINT_DATA " << im * jm << '\n';
    Cube_vtk_zonal_File << buf.str();
    buf.str("");
    buf.clear();

    dump_zonal("u-Component", u, u_0, k_zonal, Cube_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Cube_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Cube_vtk_zonal_File);

    double inv_u_0 = 1.0 / u_0;

    buf << "VECTORS u-v-Cell float\n";
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            buf << u.x[i][j][k_zonal] * inv_u_0 << " "
                << v.x[i][j][k_zonal] * inv_u_0 << " " << 0.0 << '\n';
        }
    }
    Cube_vtk_zonal_File << buf.str();
    buf.str("");
    buf.clear();

    // height is invariant across j, hoist get_layer_height out of j-loop
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        double height = get_layer_height(i);
        for(int j = 0; j < jm; j++){
            aux_t.x[i][j][k_zonal] = height;
        }
    }


    dump_zonal("Topography", h, 1.0, k_zonal, Cube_vtk_zonal_File);
    dump_zonal("height", aux_t, 1e-3, k_zonal, Cube_vtk_zonal_File);

    dump_zonal("p_dyn", p_dyn, p_0, k_zonal, Cube_vtk_zonal_File);

    if(turb_model != "none"){
        dump_zonal("TKE",          tke,        1e3, k_zonal, Cube_vtk_zonal_File);
        dump_zonal("Dissipation",  dis,        1.0, k_zonal, Cube_vtk_zonal_File);
        dump_zonal("EddyViscosity", nue,       1e3, k_zonal, Cube_vtk_zonal_File);
        dump_zonal("TurbProd",     prod,       1e3, k_zonal, Cube_vtk_zonal_File);
        dump_zonal("TKE_Source",   tke_source, 1e3, k_zonal, Cube_vtk_zonal_File);
        dump_zonal("Dis_Source",   dis_source, 1.0, k_zonal, Cube_vtk_zonal_File);
    }

    Cube_vtk_zonal_File << buf.str();


    Cube_vtk_zonal_File.close();
    cout << "   File:  " << Cube_zonal_File_Name
        << "  has been written\n";
}
/*
 * 
*/
void cCubeModel::paraview_vtk_longal(string &Name_Bathymetry_File, 
    int j_longal, int n){
    using namespace ParaViewAtm;
    string Cube_longal_File_Name = output_path + "/" + Name_Bathymetry_File
        + "Cube_turbulent_longal_" + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    ofstream Cube_vtk_longal_File;
    Cube_vtk_longal_File.precision(4);
    Cube_vtk_longal_File.setf(ios::fixed);
    Cube_vtk_longal_File.open(Cube_longal_File_Name);
    if(!Cube_vtk_longal_File.is_open()){
        cerr << "ERROR: could not open vtk_longal file " << __FILE__
            << " at line " << __LINE__ << "\n";
        abort();
    }

    ostringstream buf;
    buf.precision(4);
    buf.setf(ios::fixed);

    buf << "# vtk DataFile Version 3.0\n"
        << "Longitudinal_Data_Cube_Circulation\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_GRID\n"
        << "DIMENSIONS " << km << " " << im << " " << 1 << '\n'
        << "POINTS " << im * km << " float\n";

//    double dx = 0.1;
    double dx = 0.05;
    double dz = 0.025;

    for(int i = 0; i < im; i++){
        double x = i * dx;
        for(int k = 0; k < km; k++){
            buf << x << " " << 0.0 << " " << k * dz << '\n';
        }
    }

    buf << "POINT_DATA " << im * km << '\n';
    Cube_vtk_longal_File << buf.str();
    buf.str("");
    buf.clear();
    dump_longal("u-Component", u, u_0, j_longal, Cube_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Cube_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Cube_vtk_longal_File);

    Cube_vtk_longal_File << buf.str();
    buf.str("");
    buf.clear();

    for(int i = 0; i < im; i++){
        double height = get_layer_height(i);
        for(int k = 0; k < km; k++){
            aux_t.x[i][j_longal][k] = height;
        }
    }

    dump_longal("Topography", h, 1.0, j_longal, Cube_vtk_longal_File);
    dump_longal("height", aux_t, 1e-3, j_longal, Cube_vtk_longal_File);

    dump_longal("p_dyn", p_dyn, p_0, j_longal, Cube_vtk_longal_File);


    if(turb_model != "none"){
        dump_longal("TKE",           tke,        1e3, j_longal, Cube_vtk_longal_File);
        dump_longal("Dissipation",   dis,        1.0, j_longal, Cube_vtk_longal_File);
        dump_longal("EddyViscosity", nue,        1e3, j_longal, Cube_vtk_longal_File);
        dump_longal("TurbProd",      prod,       1e3, j_longal, Cube_vtk_longal_File);
        dump_longal("TKE_Source",    tke_source, 1e3, j_longal, Cube_vtk_longal_File);
        dump_longal("Dis_Source",    dis_source, 1.0, j_longal, Cube_vtk_longal_File);
    }


    buf << "VECTORS u-w-Cell float\n";
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            buf << u.x[i][j_longal][k]
                << " " << 0.0 << " " << w.x[i][j_longal][k] << '\n';
        }
    }
    Cube_vtk_longal_File << buf.str();
    buf.str("");
    buf.clear();


    Cube_vtk_longal_File << buf.str();


    Cube_vtk_longal_File.close();
    cout << "   File:  " << Cube_longal_File_Name
        << "  has been written\n";
}
