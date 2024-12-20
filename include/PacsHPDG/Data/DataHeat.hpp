/**
 * @file DataHeat.hpp
 * @author Lorenzo Citterio (github.com/Citte00)
 * @brief Data for the Heat equation problem.
 * @date 2024-12-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef INCLUDE_PACSHPDG_DATA_DATAHEAT_HPP
#define INCLUDE_PACSHPDG_DATA_DATAHEAT_HPP

#include "../Base.hpp"
#include "../Fem.hpp"
#include "../Algebra.hpp"
#include "../Geometry.hpp"

#include <iostream>
#include <array>

namespace pacs {

    struct DataHeat {

        // Geometrical properties
        std::vector<Point> domain = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
        int N = 300;
        bool meshFromFile = true;
        std::string VTKMeshFileName = "Mesh.vtk";
        std::string meshFileSeq = "meshes/square/square_300.poly";

        // Material properties
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real> D_ext = [](Vector<Real> x, Vector<Real> y, Real t) { return 1.0 + 0.0 * x; };
        
        // Forcing Term
        bool homog_source_f = false;
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real, Vector<Real>> source_f = [](Vector<Real> x, Vector<Real> y, Real t, Vector<Real> D) {
            Vector<Real> result{x.size()};

            for (size_t i = 0; i < x.size(); i++)
                result[i] = - (std::cos(2.0*M_PI*x[i])*std::cos(2*M_PI*y[i])+2) + 8*M_PI*M_PI*D[i]*std::cos(2*M_PI*x[i])*std::cos(2*M_PI*y[i])*(1-t);
            
            return result;
        };

        // Boundary Conditions
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real> DirBC = [](Vector<Real> x, Vector<Real> y, Real t) {
                
            Vector<Real> result{x.size()};

            for (size_t i = 0; i < x.size(); i++)
                result[i] = (std::cos(2.0*M_PI*x[i]) * std::cos(2.0*M_PI*y[i]) + 2.0) * (1-t);

            return result;
        };

        // Gradients of the Boundary Conditions
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real> DirBC_dx = [](Vector<Real> x, Vector<Real> y, Real t) {
                
            Vector<Real> result{x.size()};

            for (size_t i = 0; i < x.size(); i++)
                result[i] = -2.0L * M_PI * std::sin(2.0*M_PI*x[i]) * std::cos(2.0*M_PI*y[i]) * (1-t);

            return result;
        };

        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real> DirBC_dy = [](Vector<Real> x, Vector<Real> y, Real t) {
                
            Vector<Real> result{x.size()};

            for (size_t i = 0; i < x.size(); i++)
                result[i] = -2.0L * M_PI * std::cos(2.0*M_PI*x[i]) * std::sin(2.0*M_PI*y[i]) * (1-t);

            return result;
        };
        
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>> DirBC_dt = [](Vector<Real> x, Vector<Real> y) {
                
            Vector<Real> result{x.size()};

            for (size_t i = 0; i < x.size(); i++)
                result[i] = - (std::cos(2.0*M_PI*x[i]) * std::cos(2.0*M_PI*y[i]) + 2.0);

            return result;
        };

        // Exact Solution
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real> c_ex = DirBC;

        // Gradients of the Exact Solution
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real> dc_dx_ex = DirBC_dx;
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>, Real> dc_dy_ex = DirBC_dy;
        GenFunc<Vector<Real>, Vector<Real>, Vector<Real>> dc_dt_ex = DirBC_dt;

        // Time discretization
        Real t_0 = 0.0;
        Real t_f = 2.0;
        Real dt = 1e-2;
        Real theta = 0.5;

        // Space discretization
        size_t degree = 4;
        Real penalty_coeff = 10.0;

        // Visualization settings
        bool PlotExact = true;
        bool PlotGridSol = true;
        bool PlotIniCond = true;
        int VisualizationStep = 10;
        int NqnVisualization = 5;

        // p-adaptivity
        bool isAdaptive = true;
        std::string p_adaptive = "eta";
        double grad_threshold = -INFINITY;
        double eta_threshold = -INFINITY;
        double error_threshold = -INFINITY;

        // Save solution settings
        double SaveSolutionStep = 0.05;
        std::string VTKpFileName = "Distribution_p_" + std::to_string(degree) + "_t_";

    };

}


#endif