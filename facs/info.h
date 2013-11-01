#ifndef INFO
#define INFO

enum info_errors {
  info_no_error=0,
  info_bad_call,
  info_error_reading_bloom_filter
};

int info_main(int argc, char** argv);

#endif
