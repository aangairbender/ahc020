#pragma once

#include <algorithm>
#include <cassert>
#include <set>
#include <array>
#include<chrono>
#include "model.h"
#include "dsu.h"

struct Matching {
	int sid;
	double dist;
};

class Solver {
	const Input& m_input;

	std::vector<Edge> m_sorted_edges;
	std::vector<std::vector<Matching>> m_resident_matchings;
	std::array<std::vector<std::pair<int, int>>, N> m_station_res;
	std::vector<std::vector<Edge>> g;
	std::vector<int> assignment;
public:
	Solver(const Input& input) : m_input(input) {
		m_sorted_edges = input.edges;
		std::sort(m_sorted_edges.begin(), m_sorted_edges.end(), [](const Edge& e1, const Edge& e2) {
			return e1.w < e2.w;
		});

		m_resident_matchings.resize(input.residents.size());
		for (int i = 0; i < input.residents.size(); ++i) {
			auto& r = input.residents[i];
			for (int j = 0; j < input.stations.size(); ++j) {
				auto& s = input.stations[j];
				double dist2 = calc_dist2(r, s);
				if (dist2 <= P_MAX * P_MAX) {
					m_resident_matchings[i].push_back(Matching{j, sqrt(dist2)});
				}
			}
			std::sort(m_resident_matchings[i].begin(), m_resident_matchings[i].end(), [](auto a, auto b) {
				return a.dist < b.dist;
			});
		}

		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < input.residents.size(); ++j) {
				int d2 = calc_dist2(input.stations[i], input.residents[j]);
				if (d2 <= P_MAX * P_MAX) {
					m_station_res[i].push_back({d2, j});
				}
			}
			std::sort(m_station_res[i].begin(), m_station_res[i].end());
		}

		g.resize(N);
		for (auto e : input.edges) {
			g[e.a].push_back(e);
			g[e.b].push_back(e);
			std::swap(g[e.b].back().a, g[e.b].back().b);
		}
	}

	void solve(Solution& output) {
		init_strengths(output);
		optimize_strengths(output);
		optimize_strengths2(output);
		mst(output.station_strengths, output.edge_activations);
	}

	void init_strengths(Solution& output) {
		assignment.resize(m_input.residents.size());
		for (int i = 0; i < m_input.residents.size(); ++i) {
			double best_penalty = 1e9;
			int best_x = 0;
			int best_sid = -1;
			for (auto& m : m_resident_matchings[i]) {
				double r = output.station_strengths[m.sid];
				int x = std::max(0, (int)ceil(m.dist - r));
				double penalty = (2.0 * r + x) * x;

				if (best_penalty > penalty) {
					best_penalty = penalty;
					best_x = x;
					best_sid = m.sid;
				}
			}
			assert(best_sid != -1);
			if (best_penalty > 0) {
				output.station_strengths[best_sid] += best_x;
			}
			assignment[i] = best_sid;
		}
	}

	void recalc_rm(Solution& output, int sid, int rid, const std::vector<int>& assignment) {
		int max_r2 = 0;
		for (auto [d2, i] : m_station_res[sid]) {
			if (assignment[i] == sid &&  d2 > max_r2) max_r2 = d2;
		}
		output.station_strengths[sid] = (int)ceil(sqrt(max_r2));
	}

	void recalc_add(Solution& output, int sid, int rid, const std::vector<int>& assignment) {
		output.station_strengths[sid] = std::max(output.station_strengths[sid], (int)ceil(sqrt(calc_dist2(m_input.stations[sid], m_input.residents[rid]))));
	}

	void optimize_strengths(Solution& output) {
		std::vector<int> cov_cnt(m_input.residents.size(), 0);
		for (int i = 0; i < N; ++i) {
			auto s = m_input.stations[i];
			auto cur_r = output.station_strengths[i];
			for (int j = 0; j < m_input.residents.size(); ++j) {
				auto r = m_input.residents[j];
				auto dist2 = calc_dist2(r, s);
				if (dist2 <= cur_r * cur_r) cov_cnt[j]++;
			}
		}
		
		std::vector<int> order;
		for (int i = 0; i < N; ++i) {
			order.push_back(i);
		}
		std::sort(order.begin(), order.end(), [&](int a, int b) {
			return output.station_strengths[a] > output.station_strengths[b];
		});

		for (int i : order) {
			auto s = m_input.stations[i];
			long long max_dist2 = 0;
			int max_j = -1;
			auto old_r = output.station_strengths[i];
			for (int j = 0; j < m_input.residents.size(); ++j) {
				auto r = m_input.residents[j];
				auto dist2 = calc_dist2(r, s);
				if (dist2 <= old_r * old_r && cov_cnt[j] == 1) {
					long long dist2 = calc_dist2(s, r);
					if (max_dist2 < dist2) {
						max_dist2 = dist2;
						max_j = j;
					}
				}
			}

			output.station_strengths[i] = (int)std::ceil(std::sqrt(max_dist2));
			auto cur_r = output.station_strengths[i];

			for (int j = 0; j < m_input.residents.size(); ++j) {
				auto r = m_input.residents[j];
				auto dist2 = calc_dist2(r, s);
				if (dist2 <= old_r * old_r && dist2 > cur_r * cur_r) cov_cnt[j]--;
			}
		}
	}

	void optimize_strengths2(Solution& output) {
		auto best = output;
		auto best_score = eval(m_input, output);
		std::chrono::time_point start = std::chrono::steady_clock::now();
		int iters = 0;
		while (std::chrono::steady_clock::now() - start < std::chrono::milliseconds(1500)) {
			int i = rand() % m_input.residents.size();
			if (m_resident_matchings[i].size() == 1) continue;
			int b = m_resident_matchings[i][rand() % m_resident_matchings[i].size()].sid;
			if (assignment[i] == b) continue;
			int a = assignment[i];
			assignment[i] = b;
			recalc_rm(output, a, i, assignment);
			recalc_add(output, b, i, assignment);

			auto score = eval(m_input, output);
			if (score < best_score) {
				best_score = score;
				best = output;
			}
			iters++;
		}
		std::cerr << iters << std::endl;
		output = best;
	}

	void mst(const std::vector<int>& stations_used, std::vector<bool>& edge_activations) {
		std::fill(edge_activations.begin(), edge_activations.end(), true);
		Dsu dsu(N);
		for (auto& e : m_sorted_edges) {
			if (dsu.same(e.a, e.b)) {
				edge_activations[e.id] = false;
			}
			else {
				dsu.merge(e.a, e.b);
			}
		}

		optimize_mst(0, -1, stations_used, edge_activations);
	}

	bool optimize_mst(int v, int p, const std::vector<int>& stations_used, std::vector<bool>& edge_activations) {
		bool needed = stations_used[v] != 0;
		for (auto& e : g[v]) {
			if (e.b == p) continue;
			if (edge_activations[e.id] == false) continue;
			needed |= optimize_mst(e.b, v, stations_used, edge_activations);
		}

		if (!needed) {
			for (auto& e : g[v]) {
				edge_activations[e.id] = false;
			}
		}
		return needed;
	}
};