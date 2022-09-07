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
#include "../include/loops.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;
igl::Timer timer;
igl::triangle::SCAFData hybrid_data;

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

        igl::triangle::scaf_solve(hybrid_data, 1);
        std::cout << "time = " << timer.getElapsedTime() << std::endl;
    }

    const auto& V_uv = uv_scale * hybrid_data.w_uv;// .topRows(V.rows());
    if (show_uv)
    {
        Eigen::MatrixXi Topo;
        igl::cat(1,F, hybrid_data.s_T, Topo);
        viewer.data().clear();
        viewer.data().set_mesh(V_uv, Topo);
        viewer.data().set_uv(V_uv);
        viewer.core().align_camera_center(V_uv, Topo);
        //Eigen::MatrixXd color_map(Topo.rows(),3);
        //for (int i = 0; i < Topo.rows(); i++) {
        //    color_map(i, 0) = 1;
        //    color_map(i, 1) = 1;
        //    color_map(i, 2) = 1;
        //}
        //viewer.data().set_colormap(color_map);
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
    RIGIDT::readSTLbnd("IMR.wrl",V1,F1);

    
    auto L=RIGIDT::findLoop(F1);
    auto loop_group=RIGIDT::orientLoop(L, V1);
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F2;
    RIGIDT::reoriganize(loop_group,V1,V2,F2);



    //igl::readOBJ(string(PA)+"/camel_b.obj", V, F);
    hybrid_data.mesh=std::shared_ptr<NormalPrismaticMesh>(new NormalPrismaticMesh(V2, F2, 0.01, 1.1));


    hybrid_data.mesh->getCylinderMesh(V,F);
    hybrid_data.mesh->getUVMesh(V_uvp, F_uvp);
   
    
    Eigen::MatrixXd bnd_uv, uv_init;

    Eigen::VectorXd M;
    igl::doublearea(V, F, M);
    hybrid_data.all_bnds = hybrid_data.mesh->getAllBound();

   

    V_uvp.conservativeResize(V.rows(), 2);
    uv_init = V_uvp;
    std::vector<int> fix_id(V1.size());
    for (int i = 0; i < fix_id.size(); i++) { fix_id[i] = i; }
    igl::list_to_matrix(fix_id, hybrid_data.fixed_ids);
    
    


    Eigen::VectorXi b; Eigen::MatrixXd bc;
    igl::triangle::scaf_precompute(V, F, uv_init, hybrid_data, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);



    igl::triangle::scaf_solve(hybrid_data, 2);

    Eigen::MatrixXd V_final = hybrid_data.w_uv.topRows(V.rows());
    hybrid_data.mesh->saveVTK("test_quad.vtk", V_final);


    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    const auto& V_uv = uv_scale * hybrid_data.w_uv.topRows(V.rows());
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
