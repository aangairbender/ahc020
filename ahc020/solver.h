#pragma once

#include <algorithm>
#include <cassert>
#include <set>
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
	std::vector<std::vector<Edge>> g;
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

		g.resize(N);
		for (auto e : input.edges) {
			g[e.a].push_back(e);
			g[e.b].push_back(e);
			std::swap(g[e.b].back().a, g[e.b].back().b);
		}
	}

	void solve(Solution& output) {
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
		}

		auto output1 = output;
		mst1(output1.station_strengths, output1.edge_activations);
		auto output2 = output;
		mst2(output2.station_strengths, output2.edge_activations);

		if (eval(m_input, output1) < eval(m_input, output2)) {
			output = output1;
		}
		else {
			output = output2;
		}
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

	void mst1(const std::vector<int>& stations_used, std::vector<bool>& edge_activations) {
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

	void mst2(const std::vector<int>& stations_used, std::vector<bool>& edge_activations) {
		const int INF = 1000000000;

		std::vector<int> min_e(N, INF), sel_e(N, -1);
		min_e[0] = 0;
		std::set<std::pair<int, int> > q;
		q.insert({ 0, 0 });
		while (!q.empty()) {
			int v = q.begin()->second;
			q.erase(q.begin());

			if (sel_e[v] != -1) {
				edge_activations[sel_e[v]] = true;
			}

			for (size_t j = 0; j < g[v].size(); ++j) {
				int to = g[v][j].b;
				if (stations_used[to] == 0) continue;
				int cost = g[v][j].w;
				int eid = g[v][j].id;
				if (cost < min_e[to]) {
					q.erase({ min_e[to], to });
					min_e[to] = cost;
					sel_e[to] = eid;
					q.insert({ min_e[to], to });
				}
			}
		}
	}
};