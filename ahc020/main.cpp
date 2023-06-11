#include <iostream>
#include <cmath>
#include "model.h"
#include "solver.h"

int main() {
	int n, m, k;
	std::cin >> n >> m >> k;

	Input input;
	input.stations.resize(n);
	input.edges.resize(m);
	input.residents.resize(k);

	for (auto& s : input.stations) {
		std::cin >> s.x >> s.y;
	}
	for (int i = 0; i < input.edges.size(); ++i) {
		auto& e = input.edges[i];
		e.id = i;
		std::cin >> e.a >> e.b >> e.w;
		// fixing 1-indexing
		e.a--;
		e.b--;
	}
	for (auto& r : input.residents) {
		std::cin >> r.x >> r.y;
	}

	Solution output;
	output.station_strengths.resize(n, 0);
	output.edge_activations.resize(m, false);

	Solver solution(input);
	solution.solve(output);

	for (auto p : output.station_strengths) {
		std::cout << p << ' ';
	}
	std::cout << '\n';
	for (auto b : output.edge_activations) {
		std::cout << (int)b << ' ';
	}
	std::cout << '\n';

	std::cout.flush();
}