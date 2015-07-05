#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cassert>
#include <string>

int LoadText_Fscanf(const char* path) {
	int nrows = 0;
	std::vector<double> line(5);
	FILE* fd = fopen(path, "r");
	assert(fd);
 	while (!feof(fd)) {
		fscanf(fd, "%lf %lf %lf %lf %lf\n",
					 &line[0], &line[1], &line[2], &line[3], &line[4]);
		++nrows;
	}
	fclose(fd);
	return nrows;
}

int main(int argc, char **argv) {

	int nrows;
    nrows = LoadText_Fscanf(argv[2]);

	return 0;
}
