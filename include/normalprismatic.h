#ifndef _NORMAL_PRISM_H_
#define _NORMAL_PRISM_H_

#include <Eigen/Core>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <cstdio>


class NormalPrismaticMesh {
public:
	NormalPrismaticMesh(Eigen::MatrixXd& v, Eigen::MatrixXi& f, double fll = 0.001, double rt = 1.2) {
		F_ = f;
		V_ = (v);
		first_layer_length_ = fll;
		ratio_ = rt;
		calculateLayer();
	}

	std::vector<std::vector<int>> getAllBound();
	void getCylinderMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c);
	void getParaMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c);
	
public:
	
	int calculateLayer();
	int dim();
public:
	Eigen::MatrixXd V_;
	Eigen::MatrixXi F_;
	double first_layer_length_;

	double ratio_;
	
	int nlayer_;

};


std::vector<std::vector<int>> NormalPrismaticMesh::getAllBound() {
	std::vector<std::vector<int>> ans(2);

	

	double height = 1e-5;
	std::vector<std::array<int, 3>> f_c;
	std::vector<std::array<double, 3>> v_c;

	std::vector<std::array<double, 2>> normals;

	for (int j = 0; j < V_.rows(); j++) {
		v_c.push_back({ V_(j,0),V_(j,1) });
	}

	double sum = 0;
	int i = nlayer_-1;
	double buff = sum + pow(ratio_, i) * first_layer_length_;
	sum = buff;
	int lower_offset = i * V_.rows();
	int index_offset = (i + 1) * V_.rows();

	for (int j = 0; j < F_.rows(); j++) {
		ans[0].push_back(F_(j, 0));
		ans[1].push_back(F_(j, 0) + index_offset);

	}
	return ans;
		


}


void NormalPrismaticMesh::getCylinderMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c)
{
	double height = 1e-5;
	std::vector<std::array<int, 3>> f_c;
	std::vector<std::array<double, 3>> v_c;

	std::vector<std::array<double, 2>> normals;

	std::vector<std::pair<int, int>> pre_next(V_.rows());
	for (int i = 0; i < F_.rows(); i++) {
		pre_next[F_(i, 0)].first = F_(i, 1);
		pre_next[F_(i, 1)].second = F_(i, 0);
	}

	for (int j = 0; j < V_.rows(); j++) {
		int pre = pre_next[j].first;
		int next= pre_next[j].second;


		double x1 = V_(j,0) - V_(next,0);
		double y1 = V_(j, 1) - V_(next, 1);
		double x2 = V_(pre, 0) - V_(j, 0);
		double y2 = V_(pre, 1) - V_(j, 1);

		double x3 = -y1;
		double x4=- y2;
		double y3 = x1; 
		double y4=x2; 
		double s3 = sqrt(x3 * x3 + y3 * y3);
		double s4 = sqrt(x4 * x4 + y4 * y4);
		x3 /= s3; y3 /= s3; x4 /= s4; y4 /= s4;


		normals.push_back({ x3+x4,y3+y4 });
	}

	for (int j = 0; j < V_.rows(); j++) {
		v_c.push_back({ V_(j,0),V_(j,1) });
	}

	double sum = 0;
	for (int i = 0; i < nlayer_; i++) {
		double buff = sum + pow(ratio_, i) * first_layer_length_;
		sum = buff;
		int lower_offset = i * V_.rows();
		int index_offset = (i + 1) * V_.rows();
		for (int j = 0; j < V_.rows(); j++) {
			v_c.push_back({ V_(j,0) + buff * height * normals[j][0],V_(j,1) + buff * height * normals[j][1],0 });
		}
		for (int j = 0; j < V_.rows(); j++) {
			f_c.push_back({ F_(j,0)+ lower_offset,F_(j,1)+ lower_offset,F_(j,1) + index_offset });
			f_c.push_back({ F_(j,0)+ lower_offset,F_(j,1)+ index_offset,F_(j,0) + index_offset });

			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,0) + index_offset });
			f_c.push_back({ F_(j,1) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });
		}
	}

	bool V_rect = igl::list_to_matrix(v_c, V_c);


	bool F_rect = igl::list_to_matrix(f_c, F_c);


}


void NormalPrismaticMesh::getParaMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c)
{

	std::vector<std::array<int, 3>> f_c;
	std::vector<std::array<double, 3>> v_c;

	for (int j = 0; j < V_.rows(); j++) {
		v_c.push_back({ V_(j,0),V_(j,1),0 });
	}

	double sum = 0;
	for (int i = 0; i < nlayer_; i++) {
		double buff = sum + pow(ratio_, i) * first_layer_length_;
		sum = buff;
		int index_offset = (i + 1) * V_.rows();
		int lower_offset = i * V_.rows();
		for (int j = 0; j < V_.rows(); j++) {
			v_c.push_back({ V_(j,0),V_(j,1),buff });
		}
		for (int j = 0; j < V_.rows(); j++) {
			f_c.push_back({ F_(j,0)+ lower_offset,F_(j,1)+ lower_offset,F_(j,1) + index_offset });
			f_c.push_back({ F_(j,0)+ lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });

			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,0) + index_offset });
			f_c.push_back({ F_(j,1) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });
		}
	}

	bool V_rect = igl::list_to_matrix(v_c, V_c);


	bool F_rect = igl::list_to_matrix(f_c, F_c);

}


int NormalPrismaticMesh::calculateLayer()
{
	double length = 0;
	for (int i = 1; i < V_.rows(); i++) {
		length += sqrt(pow(V_(i,0) - V_(i - 1,0),2)+ pow(V_(i, 1) - V_(i - 1, 1), 2));;
	}
	length /= V_.rows() - 1;
	nlayer_ = log(length / first_layer_length_) / log(ratio_);
	if (!nlayer_)
		nlayer_ = 2;
	return 0;
}


int NormalPrismaticMesh::dim()
{
	return F_.rows();
}

#endif //!_NORMAL_PRISM_H_

