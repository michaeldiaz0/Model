void initialize_globals();

void initialize_perturbation();
void initialize_basic_state();
void setup_memory_allocation();
void initialize_serial();
void initialize_subarray(int,int,int);

void reinitialize();
void reinitialize_perturbation(int,double);
void init_friction();
void init_topography();
void initialize_coriolis_parameter(double *lats);


void init_basic_state_vertical_velocity();
void init_basic_state_meridional_velocity(int);
void initialize_vertical_basic_state(int,int);
void initialize_vertical_basic_state2(double*,double*);
double get_QV_Sat(double temperature,double pressure);
void change_humidity();
void stretched_grid(double*,double*,double*,double*,double,int);
void reverse_yz_coord(double*, int, int, int);
void flip_array(double *, int, int, int);

void interpolate_terrain(const char * infile);
void setup_coordinates(double,double);
void initialize_basic_state_from_output_file(const char *myfilename);
void print_vertical_basic_state();
void setup_vertical_height_levels();
void init_topography_to_zero();
void initialize_basic_state_from_reanalysis();

