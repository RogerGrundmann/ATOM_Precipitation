/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/


#include <string>
#include <fstream>
#include <cmath>

#include "cHydrosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

namespace ParaViewHyd{
    inline double safe_val(double v){ return std::isfinite(v) ? v : 0.0; }

    void dump_array(const string &name, Array &a, double multiplier, ofstream &f){
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(int k = 0; k < a.km; k++){
            for(int j = 0; j < a.jm; j++){
                for(int i = 0; i < a.im; i++){
                    f << safe_val(a.x[i][j][k] * multiplier) << endl;
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
    }
/*
*
*/
    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
*
*/
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.y[j][k] * multiplier) << endl;
            }
        }
    }
/*
*
*/
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << endl;
        f <<  "LOOKUP_TABLE default" << endl;
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << safe_val(a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
*
*/
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int i = 0; i < a.im; i++){
            for(int k = 0; k < a.km; k++){
                f << safe_val(a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::paraview_sphere_vts(const string &Name_Bathymetry_File, int n){

    using namespace ParaViewHyd;

    double x, y, z, sinthe, sinphi, costhe, cosphi;

    string Hydrosphere_sphere_vts_File_Name = output_path + "/" + Name_Bathymetry_File + "_Hyd_shere_" + std::to_string(n) + ".vts";

    ofstream Hydrosphere_sphere_vts_File;

    Hydrosphere_sphere_vts_File.precision(4);
    Hydrosphere_sphere_vts_File.setf(ios::fixed);
    
    Hydrosphere_sphere_vts_File.open(Hydrosphere_sphere_vts_File_Name);

    if(!Hydrosphere_sphere_vts_File.is_open()){
        cerr << "ERROR: could not open paraview_sphere_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    Hydrosphere_sphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Hydrosphere_sphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Hydrosphere_sphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_sphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_sphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Bathymetry Temperature PressureDyn Salinity\">\n"  << endl;
    Hydrosphere_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;

    for(int k = 0; k < km; k++){
        sinphi = sin(phi.z[k]);
        cosphi = cos(phi.z[k]);

        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);

            for(int i = 0; i < im; i++){
                aux_u.x[i][j][k] = sinthe * cosphi * u.x[i][j][k] + costhe * cosphi *
                    v.x[i][j][k] - sinphi * w.x[i][j][k];
                aux_v.x[i][j][k] = sinthe * sinphi * u.x[i][j][k] + sinphi * costhe *
                    v.x[i][j][k] + cosphi * w.x[i][j][k];
                aux_w.x[i][j][k] = costhe * u.x[i][j][k] - sinthe * v.x[i][j][k];

                Hydrosphere_sphere_vts_File << aux_u.x[i][j][k] << " " << aux_v.x[i][j][k]
                    << " " << aux_w.x[i][j][k]  << endl;
            }
            Hydrosphere_sphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_sphere_vts_File <<  "\n"  << endl;
    }

    Hydrosphere_sphere_vts_File <<  "\n"  << endl;
    Hydrosphere_sphere_vts_File <<  "    </DataArray>\n" << endl;



    dump_array("Bathymetry", h, 1.0, Hydrosphere_sphere_vts_File);

    dump_array("u-component", aux_u, u_0 * 1e2, Hydrosphere_sphere_vts_File);
    dump_array("v-component", aux_v, u_0 * 1e2, Hydrosphere_sphere_vts_File);
    dump_array("w-component", aux_w, u_0 * 1e2, Hydrosphere_sphere_vts_File);

    dump_array("Temperature", t, 1.0, Hydrosphere_sphere_vts_File);
    dump_array("PressureDynamic", p_dyn, p_0, Hydrosphere_sphere_vts_File);


    Hydrosphere_sphere_vts_File <<  "   </PointData>\n" << endl;


    Hydrosphere_sphere_vts_File <<  "   <Points>\n"  << endl;
    Hydrosphere_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n" << endl;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){

                x = (rad.z[i] + 10.0) * sin(the.z[j]) * cos(phi.z[k]);
                y = (rad.z[i] + 10.0) * sin(the.z[j]) * sin(phi.z[k]);
                z = (rad.z[i] + 10.0) * cos(the.z[j]);

                Hydrosphere_sphere_vts_File << x << " " << y << " " << z << endl;
            }
            Hydrosphere_sphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_sphere_vts_File <<  "\n"  << endl;
    }

    Hydrosphere_sphere_vts_File <<  "    </DataArray>\n" << endl;

    Hydrosphere_sphere_vts_File <<  "   </Points>\n" << endl;

    Hydrosphere_sphere_vts_File <<  "  </Piece>\n" << endl;
    Hydrosphere_sphere_vts_File <<  " </StructuredGrid>\n" << endl;
    Hydrosphere_sphere_vts_File <<  "</VTKFile>\n" << endl;
    Hydrosphere_sphere_vts_File.close();


    cout << "   File:  " << Hydrosphere_sphere_vts_File_Name
        << "  has been written\n";
}
/*
*
*/
void cHydrosphereModel::paraview_panorama_vts(const string &Name_Bathymetry_File, int n){

    using namespace ParaViewHyd;

    double x, y, z, dx, dy, dz;

    string Hydrosphere_panorama_vts_File_Name = output_path + "/"
        + Name_Bathymetry_File + "_Hyd_panorama_" + std::to_string(n) + ".vts";

    ofstream Hydrosphere_panorama_vts_File;

    Hydrosphere_panorama_vts_File.precision(4);
    Hydrosphere_panorama_vts_File.setf(ios::fixed);
    Hydrosphere_panorama_vts_File.open(Hydrosphere_panorama_vts_File_Name);

    if(!Hydrosphere_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open panorama_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    Hydrosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Hydrosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Bathymetry Temperature u-velocity v-velocity w-velocity PressureDynamic Salinity DensitySaltWater BuoyancyForce CoriolisForce PresGradForce\">\n"  << endl;

    Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Hydrosphere_panorama_vts_File << u.x[i][j][k] << " " << v.x[i][j][k] << " " << w.x[i][j][k] << endl;
            }
            Hydrosphere_panorama_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_panorama_vts_File <<  "\n"  << endl;
    }

    Hydrosphere_panorama_vts_File <<  "\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


    dump_array("Bathymetry", h, 1.0, Hydrosphere_panorama_vts_File);

    dump_array("u-velocity", u, u_0 * 1e2, Hydrosphere_panorama_vts_File);
    dump_array("v-velocity", v, u_0 * 1e2, Hydrosphere_panorama_vts_File);
    dump_array("w-velocity", w, u_0 * 1e2, Hydrosphere_panorama_vts_File);

    dump_array("Temperature", t, 1.0, Hydrosphere_panorama_vts_File);

    dump_array("PressureDynamic", p_dyn, p_0, Hydrosphere_panorama_vts_File);
//    dump_array("PressureHydro", p_hydro, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("Salinity", c, c_35, Hydrosphere_panorama_vts_File);
//    dump_array("DensityWater", r_water, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("DensitySaltWater", r_salt_water, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("CoriolisForce", CoriolisForce, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("PresGradForce", PresGradForce, 1.0, Hydrosphere_panorama_vts_File);
    if (turb_model != "laminar") {
        dump_array("TKE",           tke,        1.0, Hydrosphere_panorama_vts_File);
        dump_array("Dissipation",   dis,        1.0, Hydrosphere_panorama_vts_File);
        dump_array("EddyViscosity", nue,        1.0, Hydrosphere_panorama_vts_File);
        dump_array("TurbProd",      prod,       1.0, Hydrosphere_panorama_vts_File);
        dump_array("TKE_Source",    tke_source, 1.0, Hydrosphere_panorama_vts_File);
        dump_array("Dis_Source",    dis_source, 1.0, Hydrosphere_panorama_vts_File);
    }

    Hydrosphere_panorama_vts_File <<  "   </PointData>\n" << endl;
    Hydrosphere_panorama_vts_File <<  "   <Points>\n"  << endl;


    Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n" << endl;

    x = 0.0;
    y = 0.0;
    z = 0.0;

//    dx = 0.025;
    dx = 0.1;
    dy = 0.1;
    dz = 0.1;

    // Regular structured grid: x = i*dx, y = j*dy, z = k*dz. The previous incremental
    // form forced x = 0 for the ENTIRE k=0 plane and every j=0 row (via the in-loop
    // "if(k==0||j==0) x=0.0"), collapsing those cells to zero width -> degenerate
    // (zero-volume) cells, which made VTK's gradient/contour filters fail with
    // "Unable to factor linear system / Cannot compute gradient of grid".
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                x = i * dx;  y = j * dy;  z = k * dz;
                Hydrosphere_panorama_vts_File << x << " " << y << " " << z  << endl;
            }
            Hydrosphere_panorama_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_panorama_vts_File <<  "\n"  << endl;
    }

    Hydrosphere_panorama_vts_File <<  "    </DataArray>\n"  << endl;

    Hydrosphere_panorama_vts_File <<  "   </Points>\n"  << endl;

    Hydrosphere_panorama_vts_File <<  "  </Piece>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  " </StructuredGrid>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "</VTKFile>\n"  << endl;
    Hydrosphere_panorama_vts_File.close();


    cout << "   File:  " << Hydrosphere_panorama_vts_File_Name
        << "  has been written\n";
}
/*
*
*/
void cHydrosphereModel::paraview_vtk_longal(const string &Name_Bathymetry_File, 
    int j_longal, int n){
    using namespace ParaViewHyd;
    double x, y, z, dx, dz;
    string Hydrosphere_longal_File_Name = output_path + "/" + Name_Bathymetry_File
        + "_Hyd_longal_" + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    ofstream Hydrosphere_vtk_longal_File;
    Hydrosphere_vtk_longal_File.precision(4);
    Hydrosphere_vtk_longal_File.setf(ios::fixed);
    Hydrosphere_vtk_longal_File.open(Hydrosphere_longal_File_Name);
    if(!Hydrosphere_vtk_longal_File.is_open()){
        cerr << "ERROR: could not open vtk_longal file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
    Hydrosphere_vtk_longal_File <<  "Longitudinal_Data_Hydrosphere_Circulation\n";
    Hydrosphere_vtk_longal_File <<  "ASCII" << endl;
    Hydrosphere_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Hydrosphere_vtk_longal_File <<  "DIMENSIONS " << km << " "<< im << " " << 1 << endl;
    Hydrosphere_vtk_longal_File <<  "POINTS " << im * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dz = 0.025;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(k == 0) z = 0.0;
            else z = z + dz;

            Hydrosphere_vtk_longal_File << x << " " << y << " "<< z << endl;
        }
        z = 0.0;
        x = x + dx;
    }
    Hydrosphere_vtk_longal_File <<  "POINT_DATA " << im * km << endl;

    Hydrosphere_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
    Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_longal_File << t.x[i][j_longal][k] 
                * t_0 - t_0 << endl;
        }
    }

/*
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            aux_u.x[i][j_longal][k] = (r_salt_water.x[i][j_longal][k] * cp_w 
                * t.x[i][j_longal][k] * t_0 
                + 0.5 * r_salt_water.x[i][j_longal][k] 
                * pow(sqrt((u.x[i][j_longal][k] * u.x[i][j_longal][k] 
                          + v.x[i][j_longal][k] * v.x[i][j_longal][k] 
                          + w.x[i][j_longal][k] * w.x[i][j_longal][k]) 
                          * u_0 * u_0/3.0), 2.0)) * 1.0e-2;
            aux_v.x[i][j_longal][k] = (0.5 * r_salt_water.x[i][j_longal][k] 
                * pow(sqrt((u.x[i][j_longal][k] * u.x[i][j_longal][k] 
                          + v.x[i][j_longal][k] * v.x[i][j_longal][k] 
                          + w.x[i][j_longal][k] * w.x[i][j_longal][k]) 
                          * u_0 * u_0/3.0), 2.0)) * 1.0e-2;
            aux_w.x[i][j_longal][k] = r_salt_water.x[i][j_longal][k] * cp_w 
                * t.x[i][j_longal][k] * t_0 * 1.0e-2;
        }
    }
*/
    dump_longal("Bathymetry", h, 1.0, j_longal, Hydrosphere_vtk_longal_File);

    // Layer height [m]: signed vertical coordinate, 0 at the sea surface (i=im-1),
    // negative downward to -depth at the bottom (i=0). Mirrors the atmosphere's
    // "height" field. height = (rad.z[i] - rad.z[im-1]) * L_hyd.
    Array height_l(im, jm, km, 0.0);
    for(int i = 0; i < im; i++){
        double height = (rad.z[i] - rad.z[im-1]) * L_hyd;
        for(int k = 0; k < km; k++)
            height_l.x[i][j_longal][k] = height;
    }
    dump_longal("height", height_l, 1.0, j_longal, Hydrosphere_vtk_longal_File);

    dump_longal("u-Component", u, u_0 * 1e2, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("v-Component", v, u_0 * 1e2, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("w-Component", w, u_0 * 1e2, j_longal, Hydrosphere_vtk_longal_File);
/*
    dump_longal("rhs_u", rhs_u, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("rhs_v", rhs_v, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("rhs_w", rhs_w, 1.0, j_longal, Hydrosphere_vtk_longal_File);
*/
    dump_longal("PressureDynamic", p_dyn, p_0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("PressureHydro", p_hydro, 1.0, j_longal, Hydrosphere_vtk_longal_File);
//    dump_longal("PressureTotal", aux_u, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("PressureDynamic", p_dyn, p_0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("Salinity", c, c_35, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("DensityWater", r_water, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("DensSaltWater", r_salt_water, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("CoriolisForce", CoriolisForce, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("CentrifugalForce", CentrifugalForce, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("PresGradForce", PresGradForce, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    if (turb_model != "laminar") {
        dump_longal("TKE",           tke,        1.0, j_longal, Hydrosphere_vtk_longal_File);
        dump_longal("Dissipation",   dis,        1.0, j_longal, Hydrosphere_vtk_longal_File);
        dump_longal("EddyViscosity", nue,        1.0, j_longal, Hydrosphere_vtk_longal_File);
        dump_longal("TurbProd",      prod,       1.0, j_longal, Hydrosphere_vtk_longal_File);
        dump_longal("TKE_Source",    tke_source, 1.0, j_longal, Hydrosphere_vtk_longal_File);
        dump_longal("Dis_Source",    dis_source, 1.0, j_longal, Hydrosphere_vtk_longal_File);
    }
    Hydrosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_longal_File << u.x[i][j_longal][k] << " " << y << " " << w.x[i][j_longal][k] << endl;
        }
    }

    Hydrosphere_vtk_longal_File.close();
    cout << "   File:  " << Hydrosphere_longal_File_Name
        << "  has been written\n";
}
/*
*
*/
void cHydrosphereModel::paraview_vtk_radial(const string &Name_Bathymetry_File, 
    int i_radial, int n){
    using namespace ParaViewHyd;
    double x, y, z, dx, dy;
    string Hydrosphere_radial_File_Name = output_path + "/" + Name_Bathymetry_File
        + "_Hyd_radial_" + std::to_string(i_radial) + "_" + std::to_string(n) + ".vtk";
    ofstream Hydrosphere_vtk_radial_File;
    Hydrosphere_vtk_radial_File.precision(4);
    Hydrosphere_vtk_radial_File.setf(ios::fixed);
    Hydrosphere_vtk_radial_File.open(Hydrosphere_radial_File_Name);
    if(!Hydrosphere_vtk_radial_File.is_open()){
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
    Hydrosphere_vtk_radial_File <<  "Radial_Data_Hydrosphere_Circulation\n";
    Hydrosphere_vtk_radial_File <<  "ASCII" << endl;
    Hydrosphere_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
    Hydrosphere_vtk_radial_File <<  "DIMENSIONS " << km << " "<< jm << " " << 1 << endl;
    Hydrosphere_vtk_radial_File <<  "POINTS " << jm * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.1;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(k == 0) y = 0.0;
            else y = y + dy;

            Hydrosphere_vtk_radial_File << x << " " << y << " "<< z << endl;
        }
        y = 0.0;
        x = x + dx;
    }
    Hydrosphere_vtk_radial_File <<  "POINT_DATA " << jm * km << endl;
    Hydrosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_radial_File << t.x[i_radial][j][k] * t_0 - t_0 << endl;
            aux_w.x[i_radial][j][k] = Evaporation.y[j][k] - Precipitation_2D.y[j][k];
//            if(is_land(h, i_radial, j, k))  aux_w.x[i_radial][j][k] = 0.;
        }
    }
    Hydrosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_radial_File << t.x[i_radial][j][k] * t_0 
                - t_0 << endl;
/*
            aux_u.x[i_radial][j][k] = (r_salt_water.x[i_radial][j][k] * cp_w 
                * t.x[i_radial][j][k] * t_0 
                + 0.5 * r_salt_water.x[i_radial][j][k] 
                * pow(sqrt((u.x[i_radial][j][k] * u.x[i_radial][j][k] 
                          + v.x[i_radial][j][k] * v.x[i_radial][j][k] 
                          + w.x[i_radial][j][k] * w.x[i_radial][j][k]) 
                          * u_0 * u_0/3.0), 2.0)) * 1.0e-2;
            aux_v.x[i_radial][j][k] = (0.5 * r_salt_water.x[i_radial][j][k] 
                * pow(sqrt((u.x[i_radial][j][k] * u.x[i_radial][j][k] 
                          + v.x[i_radial][j][k] * v.x[i_radial][j][k] 
                          + w.x[i_radial][j][k] * w.x[i_radial][j][k]) 
                          * u_0 * u_0/3.0), 2.0)) * 1.0e-2;
*/
        }
    }
    dump_radial("Bathymetry", h, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial_2d("Bathymetry_m", Bathymetry, 1.0, Hydrosphere_vtk_radial_File);

    dump_radial_2d("v-velocity_NASA", velocity_v_NASA, 1.0, Hydrosphere_vtk_radial_File);
    dump_radial_2d("w-velocity_NASA", velocity_w_NASA, 1.0, Hydrosphere_vtk_radial_File);

    dump_radial("u-Component", u, u_0 * 1e2, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("v-Component", v, u_0 * 1e2, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("w-Component", w, u_0 * 1e2, i_radial, Hydrosphere_vtk_radial_File);
/*
    dump_radial("rhs_u", rhs_u, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("rhs_v", rhs_v, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("rhs_w", rhs_w, 1.0, i_radial, Hydrosphere_vtk_radial_File);
*/
    dump_radial("PressureDynamic", p_dyn, p_0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("PressureHydro", p_hydro, 1.0, i_radial, Hydrosphere_vtk_radial_File);
//    dump_radial("PressureTotal", aux_u, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("PressureDynamic", p_dyn, p_0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("Salinity", c, c_35, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("DensityWater", r_water, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("DensSaltWater", r_salt_water, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("CoriolisForce", CoriolisForce, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("CentrifugalForce", CentrifugalForce, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("CentrifugalForce", CentrifugalForce, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("PresGradForce", PresGradForce, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    // EkmanPumping/Up/Downwelling are already a physical vertical velocity in
    // cm/d (coeff_pumping in ThermoHyd.h). Dump with multiplier 1.0 like the
    // other rate fields (Evaporation/Precipitation) — the old u_0*1e2 double-
    // scaled them, so ParaView disagreed with the MinMax log by ~24x.
    dump_radial_2d("EkmanPumping_cm_d", EkmanPumping, 1.0, Hydrosphere_vtk_radial_File);
    dump_radial_2d("Upwelling_cm_d", Upwelling, 1.0, Hydrosphere_vtk_radial_File);
    dump_radial_2d("Downwelling_cm_d", Downwelling, 1.0, Hydrosphere_vtk_radial_File);
    dump_radial_2d("Evaporation", Evaporation, 1.0, Hydrosphere_vtk_radial_File);
    dump_radial_2d("Precipitation", Precipitation_2D, 1.0, Hydrosphere_vtk_radial_File);
    dump_radial_2d("SalinityEvap", salinity_evaporation, c_35, Hydrosphere_vtk_radial_File);
    dump_radial("Evap-Precip", aux_w, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    if (turb_model != "laminar") {
        dump_radial("TKE",           tke,        1.0, i_radial, Hydrosphere_vtk_radial_File);
        dump_radial("Dissipation",   dis,        1.0, i_radial, Hydrosphere_vtk_radial_File);
        dump_radial("EddyViscosity", nue,        1.0, i_radial, Hydrosphere_vtk_radial_File);
        dump_radial("TurbProd",      prod,       1.0, i_radial, Hydrosphere_vtk_radial_File);
        dump_radial("TKE_Source",    tke_source, 1.0, i_radial, Hydrosphere_vtk_radial_File);
        dump_radial("Dis_Source",    dis_source, 1.0, i_radial, Hydrosphere_vtk_radial_File);
    }

    Hydrosphere_vtk_radial_File <<  "VECTORS v-w-ATOM float" << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_radial_File << v.x[i_radial][j][k] << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }

    Hydrosphere_vtk_radial_File <<  "VECTORS v-w-NASA float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_radial_File << velocity_v_NASA.y[j][k] 
                << " " << velocity_w_NASA.y[j][k] << " " << z << endl;
        }
    }

    Hydrosphere_vtk_radial_File.close();
    cout << "   File:  " << Hydrosphere_radial_File_Name
        << "  has been written\n";
}
/*
*
*/
void cHydrosphereModel::paraview_vtk_zonal(const string &Name_Bathymetry_File, 
    int k_zonal, int n){
    using namespace ParaViewHyd;
    double x, y, z, dx, dy;
    string Hydrosphere_zongal_File_Name = output_path + "/" + Name_Bathymetry_File
        + "_Hyd_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    ofstream Hydrosphere_vtk_zonal_File;
    Hydrosphere_vtk_zonal_File.precision(4);
    Hydrosphere_vtk_zonal_File.setf(ios::fixed);
    Hydrosphere_vtk_zonal_File.open(Hydrosphere_zongal_File_Name);
    if(!Hydrosphere_vtk_zonal_File.is_open()){
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
    Hydrosphere_vtk_zonal_File <<  "Zonal_Data_Hydrosphere_Circulation\n";
    Hydrosphere_vtk_zonal_File <<  "ASCII" << endl;
    Hydrosphere_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Hydrosphere_vtk_zonal_File <<  "DIMENSIONS " << jm << " "<< im << " " << 1 << endl;
    Hydrosphere_vtk_zonal_File <<  "POINTS " << im * jm << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.05;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if(j == 0) y = 0.0;
            else y = y + dy;
            Hydrosphere_vtk_zonal_File << x << " " << y << " "<< z << endl;
        }
        y = 0.0;
        x = x + dx;
    }
    Hydrosphere_vtk_zonal_File <<  "POINT_DATA " << im * jm << endl;

    Hydrosphere_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Hydrosphere_vtk_zonal_File << t.x[i][j][k_zonal] * t_0 
                - t_0 << endl;
        }
    }
    for(int i = 0; i <= im-1; i++){
        for(int j = 0; j < jm; j++){
            aux_w.x[i][j][k_zonal] = r_salt_water.x[i][j][k_zonal] * cp_w 
                * t.x[i][j][k_zonal] * t_0 * 1.0e-2;
        }
    }
    dump_zonal("Bathymetry", h, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);

    // Layer height [m]: signed vertical coordinate, 0 at the sea surface (i=im-1),
    // negative downward to -depth at the bottom (i=0). Mirrors the atmosphere's
    // "height" field; lets ParaView place/scale layers at their true depth on the
    // stretched radial grid. height = (rad.z[i] - rad.z[im-1]) * L_hyd.
    Array height_z(im, jm, km, 0.0);
    for(int i = 0; i < im; i++){
        double height = (rad.z[i] - rad.z[im-1]) * L_hyd;
        for(int j = 0; j < jm; j++)
            height_z.x[i][j][k_zonal] = height;
    }
    dump_zonal("height", height_z, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);

    dump_zonal("u-Component", u, u_0 * 1e2, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0 * 1e2, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0 * 1e2, k_zonal, Hydrosphere_vtk_zonal_File);

    dump_zonal("PressureDynamic", p_dyn, p_0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("PressureHydro", p_hydro, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
//    dump_zonal("PressureTotal", aux_u, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("PressureDynamic", p_dyn, p_0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("Energy_inner", aux_w, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("Salinity", c, c_35, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("DensityWater", r_water, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("DensSaltWater", r_salt_water, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("CoriolisForce", CoriolisForce, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("CentrifugalForce", CentrifugalForce, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("PresGradForce", PresGradForce, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    if (turb_model != "laminar") {
        dump_zonal("TKE",           tke,        1.0, k_zonal, Hydrosphere_vtk_zonal_File);
        dump_zonal("Dissipation",   dis,        1.0, k_zonal, Hydrosphere_vtk_zonal_File);
        dump_zonal("EddyViscosity", nue,        1.0, k_zonal, Hydrosphere_vtk_zonal_File);
        dump_zonal("TurbProd",      prod,       1.0, k_zonal, Hydrosphere_vtk_zonal_File);
        dump_zonal("TKE_Source",    tke_source, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
        dump_zonal("Dis_Source",    dis_source, 1.0, k_zonal, Hydrosphere_vtk_zonal_File);
    }

    Hydrosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Hydrosphere_vtk_zonal_File << u.x[i][j][k_zonal] << " " << v.x[i][j][k_zonal] << " " << z << endl;
        }
    }

    Hydrosphere_vtk_zonal_File.close();
    cout << "   File:  " << Hydrosphere_zongal_File_Name
        << "  has been written\n";
}
