#include <igl/triangle/scaf.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOBJ.h>
#include <igl/Timer.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/MappingEnergyType.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/flipped_triangles.h>
#include <igl/topological_hole_fill.h>
#include <igl/winding_number.h>

#include "../include/readDAT.h"
#include "../include/readSTLbnd.h"
#include "../include/normalprismatic.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;
igl::Timer timer;
igl::triangle::SCAFData scaf_data;

bool show_uv = false;
float uv_scale = 0.2f;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    if (key == '1')
        show_uv = false;
    else if (key == '2')
        show_uv = true;

    if (key == ' ')
    {
        timer.start();

        igl::triangle::scaf_solve(scaf_data, 10);
        std::cout << "time = " << timer.getElapsedTime() << std::endl;
    }

    const auto& V_uv = uv_scale * scaf_data.w_uv;// .topRows(V.rows());
    if (show_uv)
    {
        Eigen::MatrixXi Topo;
        igl::cat(1,F, scaf_data.s_T, Topo);
        viewer.data().clear();
        viewer.data().set_mesh(V_uv, Topo);
        viewer.data().set_uv(V_uv);
        viewer.core().align_camera_center(V_uv, Topo);
        Eigen::MatrixXd color_map(V_uv.rows(),3);
        for (int i = 0; i < V_uv.rows(); i++) {
            color_map(i, 0) = 1;
            color_map(i, 1) = 1;
            color_map(i, 2) = 1;
        }
        viewer.data().set_colormap(color_map);
    }
    else
    {
        //viewer.data().set_mesh(V, F);
        //viewer.data().set_uv(V_uv);
        //viewer.core().align_camera_center(V, F);
    }

    viewer.data().compute_normals();

    return false;
}
#define PA "C:/code/libigl-example-project/build/_deps/libigl_tutorial_tata-src"
int main(int argc, char* argv[])
{
    using namespace std;
    // Load a mesh in OFF format

    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    Eigen::MatrixXd V_uvp;
    Eigen::MatrixXi F_uvp;
    expan::readSTLbnd("ushape.wrl",V1,F1);
   // expan::readNACADAT("naca4.dat", V1, F1);
    std::map<int, std::vector<int>> next;
    std::vector<std::array<int, 2>> loop;
    for (int j = 0; j < F1.rows(); j++) {
        next[F1(j, 0)].push_back(F1(j, 1));
        next[F1(j, 1)].push_back(F1(j,0));
    }
    loop.push_back(std::array<int, 2>{ F1(0, 0),F1(0, 1) });
    for (int j = 1; j < F1.rows(); j++) {
        int p1=loop.back()[1];
        int p2= loop.back()[0];
        int p3 = next[p1][0] + next[p1][1] - p2;
        loop.push_back({ p1,p3 });
    }

    igl::list_to_matrix(loop,F1);

    
    //orient loop
    //Eigen::MatrixXd p(1,2);
    //p(0, 0) = 1e10;
    //p(0, 1) = 1e10;
    //
    //if (igl::winding_number(V1, F1, p) > 0) {
        //for (int j = 0; j < F1.rows(); j++) {
        //    swap(F1(j, 0), F1(j, 1));
        //}
    //}


    //igl::readOBJ(string(PA)+"/camel_b.obj", V, F);
    NormalPrismaticMesh mesh(V1, F1, 0.01, 1.05);
    mesh.getParaMesh(V,F);
    mesh.getCylinderMesh(V_uvp, F_uvp);
   
    
    Eigen::MatrixXd bnd_uv, uv_init;

    Eigen::VectorXd M;
    igl::doublearea(V, F, M);
  //  std::vector<std::vector<int>> all_bnds;
   // igl::boundary_loop(F, all_bnds);

   // // Heuristic primary boundary choice: longest
   // auto primary_bnd = std::max_element(all_bnds.begin(), all_bnds.end(), [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); });

   // Eigen::VectorXi bnd = Eigen::Map<Eigen::VectorXi>(primary_bnd->data(), primary_bnd->size());

   // igl::map_vertices_to_circle(V, bnd, bnd_uv);
   //// bnd_uv *= sqrt(M.sum() / (2 * igl::PI));
   // if (all_bnds.size() == 1)
   // {
   //     if (bnd.rows() == V.rows()) // case: all vertex on boundary
   //     {
   //         uv_init.resize(V.rows(), 2);
   //         for (int i = 0; i < bnd.rows(); i++)
   //             uv_init.row(bnd(i)) = bnd_uv.row(i);
   //     }
   //     else
   //     {
   //         igl::harmonic(V, F, bnd, bnd_uv, 1, uv_init);
   //         if (igl::flipped_triangles(uv_init, F).size() != 0)
   //             igl::harmonic(F, bnd, bnd_uv, 1, uv_init); // fallback uniform laplacian
   //     }
   // }
   // else
   // {
   //     // if there is a hole, fill it and erase additional vertices.
   //     all_bnds.erase(primary_bnd);
   //     Eigen::MatrixXi F_filled;
   //     igl::topological_hole_fill(F, bnd, all_bnds, F_filled);
   //     igl::harmonic(F_filled, bnd, bnd_uv, 1, uv_init);
   //     uv_init.conservativeResize(V.rows(), 2);
   // }
    scaf_data.all_bnds = mesh.getAllBound();
    V_uvp.conservativeResize(V.rows(), 2);
    uv_init = V_uvp;
    std::vector<int> fix_id(V1.size());
    for (int i = 0; i < fix_id.size(); i++) { fix_id[i] = i; }
    igl::list_to_matrix(fix_id, scaf_data.fixed_ids);
    
    


    Eigen::VectorXi b; Eigen::MatrixXd bc;
    igl::triangle::scaf_precompute(V, F, uv_init, scaf_data, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    const auto& V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
    viewer.data().set_uv(V_uv);
    viewer.callback_key_down = &key_down;

    // Enable wireframe
    viewer.data().show_lines = true;

    // Draw checkerboard texture
    viewer.data().show_texture = true;


    std::cerr << "Press space for running an iteration." << std::endl;
    std::cerr << "Press 1 for Mesh 2 for UV" << std::endl;

    // Launch the viewer
    viewer.launch();
}
