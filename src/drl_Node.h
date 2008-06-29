#ifndef __NODE_H__
#define __NODE_H__

// The node class contains information about a given node for
// use by the density server process.

// structure coord used to pass position information between
// density server and graph class

namespace drl {

class Node {

 public:
  
  bool fixed;	// if true do not change the
				// position of this node
  int id;
  
  float x,y;
  float sub_x,sub_y;
  float energy;

 public:
  
  Node( int node_id ) { x = y = 0.0; fixed = false; 
						id = node_id; }
  ~Node() { }
  
};

} // namespace drl

#endif //__NODE_H__
