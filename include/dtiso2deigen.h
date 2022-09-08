#ifndef _DTISO_EIGEN_H_
#define _DTISO_EIGEN_H_
#include "genmesh2d.h"

#include<map>
#include<set>
#include <vector>
#include <Eigen/Dense>
namespace DTISO2D {
	void triangulate(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& V_input, Eigen::MatrixXi& F_input) {
		
		int count = 0;
		std::map<int, int> m;
		for (int i = 0; i < F.rows(); i++) {

			if (m.find(F(i, 0)) == m.end())
				m[F(i, 0)] = count++;
			if (m.find(F(i, 1)) == m.end())
				m[F(i, 1)] = count++;
		}
		std::vector<int> index(F.size()); std::vector<double> coord(m.size() * 2);
		for (int i = 0; i < F.rows(); i++) {
			index[2 * i] = m[F(i, 0)]; index[2 * i+1] = m[F(i, 1)];
		}

		for (auto i : m) {

			coord[2 * i.second] = V(i.first, 0);
			coord[2 * i.second+1] = V(i.first, 1);
		}




		//reorder 




		double* coord_out; int* index_out; int nump, nume;

		GenMesh2D(coord.size()/2, index.size()/2, coord.data(), index.data(),
			&nump, &coord_out, &nume, &index_out);
		V_input.resize(nump, 2); F_input.resize(nume, 3);
		for (int i = 0; i < nump; i++) {
			V_input(i, 0) = coord_out[2 * i];
			V_input(i, 1) = coord_out[2 * i+1];
		}
		for (int i = 0; i < nume; i++) {
			F_input(i, 0) = index_out[3 * i];
			F_input(i, 1) = index_out[3 * i + 1];
			F_input(i, 2) = index_out[3 * i + 2];
		}

	}
}
#endif //£¡ _DTISO_EIGEN_H_