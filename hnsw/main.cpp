#include <string.h>
#include <stdlib.h>

void sift_test1B(char*, char*, int, int, int, int, float, int);
int main(int argc, char **argv) {
	char  data_set[200];
	int vecsize_;
    int vecdim_;
    int level_;	
	float delta_;
	int qsize_;
	int efsearch_;
	
	//strncpy(data_set, argv[1], sizeof(data_set));
	vecsize_ = atoi(argv[3]);
	vecdim_ = atoi(argv[4]);
	level_ = atoi(argv[5]);
	qsize_ = atoi(argv[6]);
	delta_ = atof(argv[7]);
	efsearch_ = atoi(argv[8]);
	
    sift_test1B(argv[1], argv[2], vecsize_, vecdim_, level_, qsize_, delta_, efsearch_);

    return 0;
};
