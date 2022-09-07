/* This head-only file defined the algorithm of loopping*/
#ifndef _LOOP_H_
#define _LOOP_H_
#include<igl/winding_number.h>

#include<vector>
#include <set>
namespace RIGIDT {
	typedef std::vector<int> loop;

	struct BLoop {
		loop outer_loop;
		std::vector<loop> inner_loop;

	};

	struct LoopGroup {
		loop boundary;// the box
		std::vector<BLoop> bloops;// the face, every one means a domain

	};

	static std::vector<loop> findLoop(const Eigen::MatrixXi& F) {
		std::vector<loop> ans;
		std::vector<std::vector<int>> next(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			next[F(i, 0)].push_back(F(i, 1));
			next[F(i, 1)].push_back(F(i, 0));
		}
		std::set<int> isread;
		int s = 0;
		while (isread.size() < F.rows()) {
			for (; s < F.rows(); s++) {
				if (isread.find(s) == isread.end())
					break;
			}
			loop l;
			l.push_back(s);
			l.push_back(next[s][0]);
			isread.insert(s);
			isread.insert(next[s][0]);

			while (true) {
				int n = next[l.back()][0] + next[l.back()][1] - l[l.size() - 2];
				if (isread.find(n) == isread.end()) {
					l.push_back(n);
					isread.insert(n);
				}
				else
					break;
			}
			ans.push_back(l);

		}
		return ans;

	}




	static double calArea(const loop& l, const Eigen::MatrixXd& V) {
		double area = 0;;
		for (int i = 0; i < l.size(); i++)
		{
			int j = (i + 1) % l.size();
			area += (double)(V(i, 0) * V(j, 1) - V(i, 1) * V(j, 0)) / 2.0;
		}
		return area;
	}
	static bool isClockwise(const loop& l, const Eigen::MatrixXd& V) {
		return calArea(l,V) < 0;
	}
	
	

	static LoopGroup orientLoop(const std::vector<loop>& input, const Eigen::MatrixXd& V) {
		LoopGroup ans;
		std::vector<double> A; 
		auto input_copy = input;
		for (int i = 0; i < input.size(); i++) {
			double a = calArea(input[i], V);
			if (a > 0)
				std::reverse(input_copy[i].begin(), input_copy[i].end());
			A.push_back(abs(a));
		}
		int index=std::max_element(A.begin(), A.end())-A.begin();
		// boundary with the biggest area
		ans.boundary = input_copy[index];
		std::reverse(ans.boundary.begin(), ans.boundary.end());
		input_copy.erase(input_copy.begin() + index);

		std::vector<std::vector<int>> subindex(input_copy.size());
		for (int i = 0; i < input_copy.size(); i++) {
			for (int j = 0; j < input_copy.size(); j++) {
				Eigen::Vector2d p(V(input[i][0],0), V(input[i][0], 1));

				Eigen::MatrixXi F;
				std::vector<std::array<int, 2>> l;
				for (int k = 0; k < input_copy[j].size(); k++) {
					l.push_back({ input_copy[j][k],input_copy[j][(k + 1) % input_copy[j].size()] });
				}
				igl::list_to_matrix(l, F);
				double v=igl::winding_number(V,F,p);
				if (v != 0) {
					subindex[j].push_back(i);
				}
			}
		}
		// If you want to implement the multiple nested ring,please implement the DAG,
		// but in here, we only onsider the situation with no more than 1-level hole
		for (int i = 0; i < subindex.size(); i++) {
			for (int j = 0; j < subindex[i].size(); j++) {
				std::reverse(input_copy[subindex[i][j]].begin(), input_copy[subindex[i][j]].end());			
			}
		}
		for (int i = 0; i < subindex.size(); i++) {
			BLoop bl;
			bl.outer_loop = input_copy[i];
			for (int j = 0; j < subindex[i].size(); j++) {
				bl.inner_loop.push_back(input_copy[subindex[i][j]]);

			}
			ans.bloops.push_back(bl);
		}
		return ans;
	}
	static void reoriganize(const LoopGroup lp, const Eigen::MatrixXd& V_input, Eigen::MatrixXd& V_output, Eigen::MatrixXi& F) {
		std::map<int, int> m; int all_size = 0;
		for (int i = 0; i < lp.bloops.size(); i++) {
			for (auto j : lp.bloops[i].outer_loop) {
				if (m.find(j) == m.end()) {
					int index = m.size();
					m[j] = index;
				}				
			}
			all_size += lp.bloops[i].outer_loop.size();
		}
		V_output.resize(m.size(), V_input.cols());
		for (auto i:m) {
			for(int j=0;j< V_input.cols();j++)
				V_output(i.second, j) = V_input(i.first, j);
		}
		F.resize(all_size, 2); int count = 0;
		for (int i = 0; i < lp.bloops.size(); i++) {
			for (auto j = 0; j < lp.bloops[i].outer_loop.size();j++) {
				F(count, 0) = m[lp.bloops[i].outer_loop[j]];
				F(count, 1) = m[lp.bloops[i].outer_loop[(j+1)% lp.bloops[i].outer_loop.size()]];
				count++;
			}
		}
	
	}
}


#endif //!_LOOP_H_