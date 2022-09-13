#include <igl/cat.h>
#include <igl/segment_segment_intersect.h>

#include "../include/dtiso2deigen.h"
#include "../include/normalprismatic.h"
void NormalPrismaticMesh::saveVTK(std::string filename, Eigen::MatrixXd& coordinate)
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

void NormalPrismaticMesh::setPreservingLayer(const double& val)
{
	preserving_layer_height_ = val;
}

void NormalPrismaticMesh::setMultipleNormal()
{
	//first we calculate the angle
	//the numbe of extra normal=(angle-180)/89
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


void NormalPrismaticMesh::getUVMesh(Eigen::MatrixXd& V_c, Eigen::MatrixXi& F_c)
{
	
	
	std::vector<std::array<int, 3>> f_c;
	std::vector<std::array<double, 3>> v_c_start;
	std::vector<std::array<double, 3>> v_c_dir;


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

		double s = sqrt((-x3 - x4) * (-x3 - x4) + (-y3 - y4) * (-y3 - y4));
		normals.push_back({( - x3 - x4)/s,( - y3 - y4)/s});
	}

	for (int j = 0; j < V_.rows(); j++) {
		v_c_start.push_back({ V_(j,0),V_(j,1) });
		v_c_dir.push_back({ 0,0,0 });
	}

	double sum = 0;
	for (int i = 0; i < nlayer_; i++) {
		double buff = sum + pow(ratio_, i) * first_layer_length_;
		sum = buff;
		int lower_offset = i * V_.rows();
		int index_offset = (i + 1) * V_.rows();
		for (int j = 0; j < V_.rows(); j++) {
			v_c_start.push_back({ V_(j,0) ,V_(j,1),0 });
			v_c_dir.push_back({ buff  * normals[j][0], buff  * normals[j][1],0 });
		}
		for (int j = 0; j < V_.rows(); j++) {
			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,1) + index_offset });
			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });

			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,0) + index_offset });
			f_c.push_back({ F_(j,1) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });
		}
	}
	Eigen::MatrixXd V_c_start;
	Eigen::MatrixXd V_c_dir;

	igl::list_to_matrix(v_c_start, V_c_start);
	igl::list_to_matrix(v_c_dir, V_c_dir);
	// devided
	double height_start = 10;
	while (true) {
		bool inter = false;
		for (int i = 0; i < F_.rows(); i++) {
			for (int j = i+1; j < F_.rows(); j++) {
				Eigen::Vector3d p1; Eigen::Vector3d dir1; Eigen::Vector3d p2; Eigen::Vector3d dir2;

				if (F_(i, 1) == F_(j, 0) || F_(i, 0) == F_(j, 1))// if share the same point
					continue;

				for (int k = 0; k < 2; k++) {
					p1(k) = v_c_start[v_c_start.size() - V_.size() + F_(i, 0)][k] + height_start*v_c_dir[v_c_dir.size() - V_.size() + F_(i, 0)][k];
					dir1(k) = v_c_start[v_c_start.size() - V_.size() + F_(i, 1)][k] + height_start*v_c_dir[v_c_dir.size() - V_.size() + F_(i, 1)][k];
					p2(k) = v_c_start[v_c_start.size() - V_.size() + F_(j, 0)][k] + height_start*v_c_dir[v_c_dir.size() - V_.size() + F_(j, 0)][k];
					dir2(k) = v_c_start[v_c_start.size() - V_.size() + F_(j, 1)][k] + height_start*v_c_dir[v_c_dir.size() - V_.size() + F_(j, 1)][k];
				}
				dir1 = dir1 - p1; dir2 = dir2 - p2;
				
				double t, u, eps = 1e-9;
				inter |=igl::segment_segment_intersect(p1,dir1,p2,dir2,t,u,eps);
			}
		}
		if (!inter)
			break;
		std::cout << height_start << std::endl;
		height_start /= 2;
	}



	V_c = V_c_start + height_start * V_c_dir;
	for (auto& i : f_c) {
		std::swap(i[0], i[1]);
	}

	bool F_rect = igl::list_to_matrix(f_c, F_c);


}
void NormalPrismaticMesh::getDiscreteCylinderMesh(Eigen::MatrixXd& V_c,
	Eigen::MatrixXi& F_c,
	Eigen::VectorXd& scale,
	std::vector<int>& discrete_map_2_nondiscrete) {


	std::vector<int>& d2nd = discrete_map_2_nondiscrete;
	std::vector<std::array<int, 3>> f_c;
	std::vector<std::array<double, 3>> v_c;




	std::vector<double> sum(F_.rows(), 0);
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
			double buff = sum[j] + q * pow(ratio_, i) * first_layer_length_;
			if (i == nlayer_ - 1 && preserving_layer_height_ > 0) {
				buff = sum[j] + q * preserving_layer_height_*pow(ratio_, i) * first_layer_length_;
			}
			sum[j] = buff;

			int curr_s = v_c.size();
			int lower_s = v_c.size() - 2 * F_.rows();

			f_c.push_back({ lower_s,lower_s + 1,curr_s + 1 });
			f_c.push_back({ lower_s,curr_s + 1,curr_s });

			f_c.push_back({ lower_s,lower_s + 1,curr_s });
			f_c.push_back({ lower_s + 1,curr_s + 1,curr_s });


			int index_offset = (i + 1) * V_.rows();


			d2nd.push_back(F_(j, 0) + index_offset);
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

void NormalPrismaticMesh::getCylinderMesh(Eigen::MatrixXd& V_c, Eigen::MatrixXi& F_c, Eigen::VectorXd& scale)
{

	std::vector<std::array<int, 3>> f_c;
	std::vector<std::array<double, 3>> v_c;

	for (int j = 0; j < V_.rows(); j++) {
		v_c.push_back({ V_(j,0),V_(j,1),0 });
	}


	double sum = 0;
	for (int i = 0; i < nlayer_; i++) {
		double buff = sum + pow(ratio_, i) * first_layer_length_;
		if (preserving_layer_height_ > 0 && i == nlayer_ - 1) {
			buff= sum + preserving_layer_height_* pow(ratio_, i) * first_layer_length_;
		}
		sum = buff;

		int index_offset = (i + 1) * V_.rows();
		int lower_offset = i * V_.rows();
		for (int j = 0; j < V_.rows(); j++) {
			double q = 1;
			if (scale.rows())
				q = scale(j);
			v_c.push_back({ V_(j,0),V_(j,1),q * buff });
		}
		for (int j = 0; j < V_.rows(); j++) {
			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + lower_offset,F_(j,1) + index_offset });
			f_c.push_back({ F_(j,0) + lower_offset,F_(j,1) + index_offset,F_(j,0) + index_offset });

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
		length += sqrt(pow(V_(i, 0) - V_(i - 1, 0), 2) + pow(V_(i, 1) - V_(i - 1, 1), 2));;
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

void NormalPrismaticMesh::getTopLoop(RIGIDT::LoopGroup lg,
	const Eigen::MatrixXd& V_boundary,
	const Eigen::MatrixXd& uv,
	Eigen::MatrixXd& newV,
	RIGIDT::LoopGroup& newl) {


	Eigen::MatrixXd V(V_.rows() + lg.boundary.size(), V_.cols()); V.topRows(V_.rows()) = V_;
	Eigen::MatrixXi F(F_.rows() + lg.boundary.size(), F_.cols()); F.topRows(F_.rows()) = F_;

	for (int i = 0; i < V_.rows(); i++) {
		V(i, 0) = uv((nlayer_ - 1) * V_.rows() + i, 0); V(i, 1) = uv((nlayer_ - 1) * V_.rows() + i, 1);
	}
	int count = V_.rows();

	for (int i = 0; i < lg.boundary.size(); i++) {
		V(count + i, 0) = V_boundary(lg.boundary[i], 0);
		V(count + i, 1) = V_boundary(lg.boundary[i], 1);
		F(count + i, 0) = count + i;
		F(count + i, 1) = count + ((i + 1) % lg.boundary.size());
	}
	auto L = RIGIDT::findLoop(F);
	newl = RIGIDT::orientLoop(L, V);
	newV = V;
}
void NormalPrismaticMesh::getHoles(const RIGIDT::LoopGroup& loops, const Eigen::MatrixXd& V_input, Eigen::MatrixXd& H)
{
	std::vector<RIGIDT::loop> ls; int total_length = 0;
	std::vector<Eigen::MatrixXd> VS;
	std::vector<Eigen::MatrixXi> FS;

	VS.resize(VS.size() + 1); FS.resize(FS.size() + 1);

	for (auto i : loops.bloops) {
		double a = RIGIDT::calArea(i.outer_loop, V_input);
		if (a < 0) {

		}
		else {
			auto p = i.outer_loop;
			//std::reverse(p.begin(), p.end());
			ls.push_back(p);
			total_length += i.outer_loop.size();
		}
	}


	H.resize(ls.size()*2, 2);
	for (int i = 0; i < ls.size(); i++) {
		auto l = ls[i];

		// find the normal
		double x1 = V_input(l[1], 0) - V_input(l[2], 0);
		double y1 = V_input(l[1], 1) - V_input(l[2], 1);
		double x2 = V_input(l[0], 0) - V_input(l[1], 0);
		double y2 = V_input(l[0], 1) - V_input(l[1], 1);
		double x3 = -y1;
		double x4 = -y2;
		double y3 = x1;
		double y4 = x2;
		double s3 = sqrt(x3 * x3 + y3 * y3);
		double s4 = sqrt(x4 * x4 + y4 * y4);
		x3 /= s3; y3 /= s3; x4 /= s4; y4 /= s4;

		H(2*i, 0) = -(x3 + x4) * 1e-5 + V_input(l[1], 0);
		H(2*i, 1) = -(y3 + y4) * 1e-5 + V_input(l[1], 1);

		H(2*i+1, 0) = (x3 + x4) * 1e-5 + V_input(l[1], 0);
		H(2*i+1, 1) = (y3 + y4) * 1e-5 + V_input(l[1], 1);

	}

}
void NormalPrismaticMesh::generateOuterTriMesh(const RIGIDT::LoopGroup& loops, const Eigen::MatrixXd& V_input, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
	std::vector<RIGIDT::loop> ls; int total_length = 0;
	std::vector<Eigen::MatrixXd> VS;
	std::vector<Eigen::MatrixXi> FS;
	std::string command = "qYYQa4";
	for (auto i : loops.bloops) {
		double a = RIGIDT::calArea(i.outer_loop, V_input);
		if (a < 0) {
			//is inner loop
			Eigen::MatrixXi E(i.outer_loop.size(), 2);
			auto p = i.outer_loop;
			std::reverse(p.begin(), p.end());
			for (int j = 0; j < i.outer_loop.size(); j++) {
				int k = (j + 1) % i.outer_loop.size();
				E(j, 0) = p[j];
				E(j, 1) = p[k];
			}
			Eigen::MatrixXd H;
			VS.resize(VS.size() + 1); FS.resize(FS.size() + 1);
			//(V_input, E, H, command, VS.back(), FS.back());


			DTISO2D::triangulate(V_input, E, VS.back(), FS.back());
		}
		else {
			auto p = i.outer_loop;
			std::reverse(p.begin(), p.end());
			ls.push_back(p);
			total_length += i.outer_loop.size();
		}
	}
	VS.resize(VS.size() + 1); FS.resize(FS.size() + 1);
	auto p = loops.boundary;
	Eigen::MatrixXi E(loops.boundary.size() + total_length, 2);
	int count = 0;
	std::reverse(p.begin(), p.end());
	for (int j = 0; j < loops.boundary.size(); j++) {
		int k = (j + 1) % loops.boundary.size();
		E(count, 0) = p[j];
		E(count, 1) = p[k];
		count++;
	}
	Eigen::MatrixXd H(ls.size(), 2);
	for (int i = 0; i < ls.size(); i++) {
		auto l = ls[i];
		for (int j = 0; j < l.size(); j++) {
			int k = (j + 1) % l.size();
			E(count, 0) = l[j];
			E(count, 1) = l[k];
			count++;
		}
		// find the normal
		double x1 = V_input(l[1], 0) - V_input(l[2], 0);
		double y1 = V_input(l[1], 1) - V_input(l[2], 1);
		double x2 = V_input(l[0], 0) - V_input(l[1], 0);
		double y2 = V_input(l[0], 1) - V_input(l[1], 1);
		double x3 = -y1;
		double x4 = -y2;
		double y3 = x1;
		double y4 = x2;
		double s3 = sqrt(x3 * x3 + y3 * y3);
		double s4 = sqrt(x4 * x4 + y4 * y4);
		x3 /= s3; y3 /= s3; x4 /= s4; y4 /= s4;

		H(i, 0) = -(x3 + x4) * 1e-5 + V_input(l[1], 0);
		H(i, 1) = -(y3 + y4) * 1e-5 + V_input(l[1], 1);

	}
	//igl::triangle::triangulate(V_input, E, H, command, VS.back(), FS.back());
	DTISO2D::triangulate(V_input, E, VS.back(), FS.back());

	Eigen::MatrixXd V_init;
	igl::combine(VS, FS, V_init, F);
	Eigen::MatrixXd VO = Eigen::MatrixXd::Zero(V_init.rows(), 1);
	igl::cat(2, V_init, VO, V);
}