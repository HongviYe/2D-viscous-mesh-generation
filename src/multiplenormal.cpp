#include <igl/matrix_to_list.h>
#include <igl/PI.h>

#include "../include/multiplenormal.h"

namespace RIGIDT {
	void VirtualMeshGeneration(const Eigen::MatrixXd  V_input, const Eigen::MatrixXi  F_input, Eigen::MatrixXd & V_output, Eigen::MatrixXi & F_output)
	{
		//first we calculate the angle
	//the numbe of extra normal=(angle-180)/91


		const auto& F_ = F_input; const auto& V_ = V_input;
		std::vector<std::array<double, 2>> point_normals;
		std::vector<std::array<double, 2>> front_normals(F_input.rows());


		std::vector<std::pair<int, int>> pre_next(V_output.rows());
		std::vector<std::pair<int, int>> front_pre_next(V_output.rows());
		for (int i = 0; i < F_.rows(); i++) {
			pre_next[F_(i, 0)].first = F_(i, 1);
			front_pre_next[F_(i, 0)].first = i;
			pre_next[F_(i, 1)].second = F_(i, 0);
			front_pre_next[F_(i, 1)].second = i;
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
			front_normals[front_pre_next[j].second] = { -x3,-y3 };	front_normals[front_pre_next[j].first] = { -x4,-y4 };

			double s = sqrt((-x3 - x4) * (-x3 - x4) + (-y3 - y4) * (-y3 - y4));
			point_normals.push_back({ (-x3 - x4) / s,(-y3 - y4) / s });
		}







		std::vector<std::vector<double>> v;
		std::vector<std::vector<int>> f;


		igl::matrix_to_list(V_input, v);
		igl::matrix_to_list(F_input, f);


		int add_count = 0;

		for (int j = 0; j < V_.rows(); j++) {
			int pre = front_pre_next[j].first;
			int next = front_pre_next[j].second;

			Eigen::Vector3d v1 = { front_normals[pre][0],  front_normals[pre][1],0 };
			Eigen::Vector3d v2 = { front_normals[next][0],front_normals[next][1],0 };
			double d = -asin(v1.cross(v2)(2));
			int num_extra_normal = d / (igl::PI / 2.9);

			if (num_extra_normal == 1) {
			
				int new_point_id = v.size();
				double scale = VIRTURALSCALE;
				v.push_back({ V_(j,0) + scale *(v1(0)), V_(j,1) + scale *(v1(1)) });
				v[j] = { V_(j,0) + scale *( v2(0)), V_(j,1) + scale *( v2(1)) };
				f[pre][0] = new_point_id;
				f.push_back({j,new_point_id});
			}
			else if (num_extra_normal == 2) {

			}
			else {
				// to be continue; 
			}
		}
		std::map<int, int> o2n;
		std::vector<int> next(v.size());
		for (int i = 0; i < f.size(); i++) {
			next[f[i][0]] = f[i][1];
		}
		while (o2n.size() < v.size()) {
			for (int i = 0; i < v.size(); i++) {
				int start = i;
				while (true) {
					if (o2n.find(start) == o2n.end()) {
						int s = o2n.size();
						o2n[start] = s;
						start = next[start];
					}
					else
						break;
				}
			}
		}

		std::vector<std::vector<double>> v_copy = v;
		for (int i = 0; i < v_copy.size();i++) {
			v[o2n[i]] = v_copy[i];
		}
		for (auto &i:f) {
			for (auto &j : i) {
				j = o2n[j];
			}
		}
	
		



		igl::list_to_matrix(v, V_output);
		igl::list_to_matrix(f, F_output);
	}
}
