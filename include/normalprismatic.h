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
	void getUVMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c);

	void getDiscreteCylinderMesh(Eigen::MatrixXd& V_c, Eigen::MatrixXi& F_c, Eigen::VectorXd& scale, std::vector<int>& discrete_map_2_nondiscrete);
	



	void getCylinderMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c, Eigen::VectorXd& scale= Eigen::VectorXd());

	void saveVTK(std::string filename, Eigen::MatrixXd& coordinate);
	void getSimpleCylinderMesh(Eigen::MatrixXd& V_c, Eigen::MatrixXi& F_c);
	
public:
	
	int calculateLayer();
	int dim();
public:
	Eigen::MatrixXd V_;
	Eigen::MatrixXi F_;
	double first_layer_length_;

	double ratio_;
	
	int nlayer_;

	double target_length_;

};
inline void NormalPrismaticMesh::saveVTK(std::string filename, Eigen::MatrixXd& coordinate)
{
	std::ofstream fout(filename);
	fout << "# vtk DataFile Version 2.0" << std::endl;
	fout << "boundary layer mesh" << std::endl;
	fout << "ASCII" << std::endl;
	fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fout << "POINTS " << coordinate.rows() * 3 << " double" << std::endl;
	int target_layer = nlayer_ - 1;
	for (int i = 0; i < coordinate.rows(); i++) {
		fout << coordinate(i, 0) << " " << coordinate(i, 1) << " " << 0 << std::endl;
	}
	fout << "CELLS " << target_layer * V_.rows() << " " << 5 * target_layer * V_.rows() << std::endl;

	std::vector<std::array<int, 4>> f_c;
	for (int i = 0; i < target_layer; i++) {

		int index_offset = (i + 1) * V_.rows();
		int lower_offset = i * V_.rows();

		for (int j = 0; j < V_.rows(); j++) {
			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });
		}
	}

	for (auto i : f_c) {
		fout << 4 << " " << i[0] << " " << i[1] << " " << i[2] << " " << i[3] << std::endl;
	}
	fout << "CELL_TYPES " << nlayer_ * V_.rows() << std::endl;
	for (int i = 0; i < nlayer_ * V_.rows(); i++) {
		fout << 9 << std::endl;
	}
	fout.close();

}

std::vector<std::vector<int>> NormalPrismaticMesh::getAllBound() {
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
		int next = pre_next[j].second;


		double x1 = V_(j, 0) - V_(next, 0);
		double y1 = V_(j, 1) - V_(next, 1);
		double x2 = V_(pre, 0) - V_(j, 0);
		double y2 = V_(pre, 1) - V_(j, 1);

		double x3 = -y1;
		double x4 = -y2;
		double y3 = x1;
		double y4 = x2;
		double s3 = sqrt(x3 * x3 + y3 * y3);
		double s4 = sqrt(x4 * x4 + y4 * y4);
		x3 /= s3; y3 /= s3; x4 /= s4; y4 /= s4;


		normals.push_back({ x3 + x4,y3 + y4 });
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
			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,1) + index_offset });
			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });
		}
	}
	Eigen::MatrixXi F_c;


	bool F_rect = igl::list_to_matrix(f_c, F_c);
	std::vector<std::vector<int>> ans;
	igl::boundary_loop(F_c, ans);
	return ans;


}


void NormalPrismaticMesh::getUVMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c)
{
	double height = 1e-3;
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


		normals.push_back({ -x3-x4,-y3-y4 });
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

	for (auto& i : f_c) {
		std::swap(i[0], i[1]);
	}

	bool F_rect = igl::list_to_matrix(f_c, F_c);


}
void NormalPrismaticMesh::getDiscreteCylinderMesh(Eigen::MatrixXd& V_c,
	Eigen::MatrixXi& F_c, 
	Eigen::VectorXd& scale,
	std::vector<int>& discrete_map_2_nondiscrete) {


	std::vector<int> &d2nd= discrete_map_2_nondiscrete;
	std::vector<std::array<int, 3>> f_c;
	std::vector<std::array<double, 3>> v_c;




	std::vector<double> sum(F_.rows(),0);
	for (int j = 0; j < F_.rows(); j++) {
		d2nd.push_back(F_(j, 0)); d2nd.push_back(F_(j, 1));

		v_c.push_back({ V_(F_(j,0),0),V_(F_(j,0),1),0 });		
		v_c.push_back({ V_(F_(j,1),0),V_(F_(j,1),1),0 });
	}
	for (int i = 0; i < nlayer_; i++) {
	for (int j = 0; j < F_.rows(); j++) {
		double q = 1;
		if (scale.rows())
			q = scale(j);	
			double buff = sum[j] + q*pow(ratio_, i) * first_layer_length_;
			sum[j] = buff;

			int curr_s = v_c.size();
			int lower_s = v_c.size() - 2 * F_.rows();

			f_c.push_back({ lower_s,lower_s+1,curr_s+1 });
			f_c.push_back({ lower_s,curr_s +1,curr_s });

			f_c.push_back({ lower_s,lower_s+1,curr_s });
			f_c.push_back({ lower_s+1,curr_s+1,curr_s });


			int index_offset = (i + 1) * V_.rows();


			d2nd.push_back(F_(j, 0)+ index_offset);
			d2nd.push_back(F_(j, 1) + index_offset);

			v_c.push_back({ V_(F_(j,0),0),V_(F_(j,0),1),buff });
			v_c.push_back({ V_(F_(j,1),0),V_(F_(j,1),1),buff });



		}
	}

	bool V_rect = igl::list_to_matrix(v_c, V_c);

	for (auto& i : f_c) {
		std::swap(i[0], i[1]);
	}

	bool F_rect = igl::list_to_matrix(f_c, F_c);
}

void NormalPrismaticMesh::getCylinderMesh(Eigen::MatrixXd & V_c, Eigen::MatrixXi& F_c,Eigen::VectorXd& scale)
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
			double q = 1;
			if (scale.rows())
				q = scale(j);
			v_c.push_back({ V_(j,0),V_(j,1),q*buff });
		}
		for (int j = 0; j < V_.rows(); j++) {
			f_c.push_back({ F_(j,0)+ lower_offset,F_(j,1)+ lower_offset,F_(j,1) + index_offset });
			f_c.push_back({ F_(j,0)+ lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });

			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,0) + index_offset });
			f_c.push_back({ F_(j,1) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });
		}
	}
	target_length_ = sum;
	bool V_rect = igl::list_to_matrix(v_c, V_c);

	for (auto& i : f_c) {
		std::swap(i[0], i[1]);
	}

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
	nlayer_++;// add a preserving layer
	return 0;
}


int NormalPrismaticMesh::dim()
{
	return F_.rows();
}

#endif //!_NORMAL_PRISM_H_

