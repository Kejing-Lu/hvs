#include <string.h>
#include <stdlib.h>

void sift_test1B(char*, int, int, int, float);
int main(int argc, char **argv) {
	char  data_set[200];
	int vecsize_;
    int vecdim_;
    int level_;	
	float delta_;
	
	//strncpy(data_set, argv[1], sizeof(data_set));
	vecsize_ = atoi(argv[2]);
	vecdim_ = atoi(argv[3]);
	level_ = atoi(argv[4]);
	delta_ = atof(argv[5]);
	
    sift_test1B(argv[1], vecsize_, vecdim_, level_, delta_);

    return 0;
};
