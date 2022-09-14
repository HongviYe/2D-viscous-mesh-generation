#ifndef _MULTIPLE_NORMAL_H_
#define _MULTIPLE_NORMAL_H_

#include <Eigen/Core>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <cstdio>

#include <igl/boundary_loop.h>
#include <igl/triangle/triangulate.h>
#include <igl/combine.h>
#include <igl/list_to_matrix.h>

/**
 * 2D virtual mesh generation.
 */
namespace RIGIDT {
	constexpr double VIRTURALSCALE = 1e-6;
	void VirtualMeshGeneration(const Eigen::MatrixXd  V_input, const Eigen::MatrixXi  F_input, Eigen::MatrixXd & V_output, Eigen::MatrixXi & F_output);
}
#endif //!_MULTIPLE_NORMAL_H_



