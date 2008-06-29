// The parse class contains the methods necessary to parse
// the command line, print help, and do error checking

#ifdef MUSE_MPI
  #include <mpi.h>
#endif

namespace drl {
	
class parse {

public:

    // Methods
	
	parse ( int argc, char **argv );
	~parse () {}
	
	// user parameters
	string sim_file;		// .sim file
	string coord_file;		// .coord file
	string parms_file;		// .parms file
	string real_file;	    // .real file
	
	int rand_seed;		// random seed int >= 0
	float edge_cut;			// edge cutting real [0,1]
	int int_out;			// intermediate output, int >= 1
	int edges_out;                  // true if .edges file is requested
	int parms_in;		    // true if .parms file is to be read
	float real_in;		    // true if .real file is to be read
	
private:

	void print_syntax ( const char *error_string );
	
};

} // namespace drl
