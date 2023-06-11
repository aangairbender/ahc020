#pragma once

#include <vector>

constexpr int N = 100;
constexpr int M_MAX = 300;
constexpr int K_MAX = 5000;
constexpr int P_MAX = 5000;

struct Vertex {
	int x, y;
};

long long calc_dist2(Vertex a, Vertex b) {
	long long dx = a.x - b.x;
	long long dy = a.y - b.y;
	return dx * dx + dy * dy;
}

struct Edge {
	int a, b, w, id;
};

struct Input {
	std::vector<Vertex> stations;
	std::vector<Edge> edges;
	std::vector<Vertex> residents;
};

struct Solution {
	std::vector<int> station_strengths;
	std::vector<bool> edge_activations;
};

long long eval(const Input& input, const Solution& solution) {
	long long res = 0;
	for (int i = 0; i < N; ++i) {
		res += solution.station_strengths[i] * solution.station_strengths[i];
	}
	for (int i = 0; i < input.edges.size(); ++i) {
		res += solution.edge_activations[i] ? input.edges[i].w : 0;
	}
	return res;
}
