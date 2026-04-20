#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

#include "cCubeModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;
 

void cCubeModel::init_topography(const string &topo_filename){

//    cout << endl << endl << endl << "      AGCM: init_topography" << endl;

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
                }
            }
        }
    }



// rewriting bathymetrical data from -180° _ 0° _ +180° coordinate system to 0°- 360°
    #pragma omp parallel for schedule(dynamic)
    for(int j = 0; j < jm; j++){
        move_data(Topography.y[j], km);
        move_data(i_topography[j], km);

        for(int i = 0; i < im; i++){
            move_data(h.x[i][j], km);
        }
    }

//    cout << "      AGCM: init_topography ended" << endl;
}
