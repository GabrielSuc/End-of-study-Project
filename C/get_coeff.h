struct getFilter {
        double* filter;
        unsigned int flen;

};

struct getFilter* get_matlab_filter(int stage_no);


void getFilter_dump(struct getFilter* getfilter);
