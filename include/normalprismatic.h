#ifndef _NORMAL_PRISM_H_
#define _NORMAL_PRISM_H_

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

#include "loops.h"




class NormalPrismaticMesh {
public:
	NormalPrismaticMesh(Eigen::MatrixXd& v, Eigen::MatrixXi& f, double fll = 0.001, double rt = 1.2) {
		F_ = f;
		V_ = (v);
		first_layer_length_ = fll;
		ratio_ = rt;
		calculateLayer(); setNormal();
	}

	std::vector<std::vector<int>> getAllBound();

	/**
	 * Get the initial mesh.
	 * 
	 * \param V_c
	 * \param F_c
	 */
	void getUVMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c);

	/**
	 * get target mesh.
	 * 
	 * \param V_c
	 * \param F_c
	 * \param scale
	 * \param discrete_map_2_nondiscrete
	 */
	void getDiscreteCylinderMesh(Eigen::MatrixXd& V_c, Eigen::MatrixXi& F_c, Eigen::VectorXd& scale = Eigen::VectorXd(), std::vector<int>& discrete_map_2_nondiscrete = std::vector<int>());
	
	void generateOuterTriMesh(const RIGIDT::LoopGroup& loops, const Eigen::MatrixXd& V_input, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

	void getTopLoop(RIGIDT::LoopGroup lg,
		const Eigen::MatrixXd& V_boundary,
		const Eigen::MatrixXd& uv,
		Eigen::MatrixXd& newV,
		RIGIDT::LoopGroup& newl);


	void getHoles(const RIGIDT::LoopGroup& loops, const Eigen::MatrixXd& V_input, Eigen::MatrixXd& H);


	void getCylinderMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c, Eigen::VectorXd& scale= Eigen::VectorXd());


	void saveVTK(std::string filename, Eigen::MatrixXd& coordinate);
	void getSimpleCylinderMesh(Eigen::MatrixXd& V_c, Eigen::MatrixXi& F_c);

	void setPreservingLayer(const double& val);

	void setMultipleNormal();
protected:
	void setNormal();
	
public:
	
	int calculateLayer();
	int dim();
public:
	Eigen::MatrixXd V_;
	Eigen::MatrixXi F_;
	double first_layer_length_;

	double ratio_;
	
	int nlayer_;

	RIGIDT::LoopGroup loops_;

	double preserving_layer_height_;
	double target_length_;

protected:
	Eigen::MatrixXd point_normals_;
	Eigen::MatrixXd front_normals_;

};



#endif //!_NORMAL_PRISM_H_

