#pragma once

#include <vector>

class Dsu {
	std::vector<int> d;
public:
	Dsu(int n) : d(n) {
		for (int i = 0; i < n; ++i) d[i] = i;
	}

	void merge(int a, int b) {
		a = color(a);
		b = color(b);
		if (rand() & 1) std::swap(a, b);
		d[a] = b;
	}

	int color(int x) {
		if (d[x] == x) return x;
		else return d[x] = color(d[x]);
	}

	inline bool same(int a, int b) {
		return color(a) == color(b);
	}
};
