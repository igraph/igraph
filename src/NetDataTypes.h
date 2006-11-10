/***************************************************************************
                          NetDataTypes.h  -  description
                             -------------------
    begin                : Mon Oct 6 2003
    copyright            : (C) 2003 by Joerg Reichardt
    email                : reichardt@mitte
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef NETDATATYPES_H
#define NETDATATYPES_H

#include <string.h>

//###########################################################################################

struct HUGE_INDEX
{
    unsigned int field_index;
    unsigned long in_field_index;
};

template <class DATA> class HugeArray
{
  private:
  unsigned long int size;
  unsigned int highest_field_index;
  unsigned long max_bit_left;
  unsigned long max_index;
  DATA *data;
  DATA *fields[32];
  public:
  HUGE_INDEX get_huge_index(unsigned long);
  DATA &Set(unsigned long);
  DATA Get(unsigned long);
  HugeArray(void);
  ~HugeArray(void);
  DATA &operator[](unsigned long);
  unsigned long Size(void) {return max_index;}
} ;
//###############################################################################################
template <class L_DATA > class DLList;
template <class L_DATA > class DL_Indexed_List;
template <class L_DATA > class ClusterList;
template <class L_DATA > class DLList_Iter;

template <class L_DATA>
class DLItem
{
      friend class DLList<L_DATA> ;
      friend class DL_Indexed_List<L_DATA>;
      friend class DLList_Iter<L_DATA>;
      private:
      L_DATA  item;
      unsigned long index;
      DLItem *previous;
      DLItem *next;
      DLItem(L_DATA i, unsigned long ind);
      DLItem(L_DATA i, unsigned long ind, DLItem<L_DATA> *p, DLItem<L_DATA> *n); 
      ~DLItem();
};

template <class L_DATA >
class DLList
{
  friend class DLList_Iter<L_DATA>;
  protected:
    DLItem<L_DATA>  *head;
    DLItem<L_DATA>  *tail;
    unsigned long number_of_items;
    DLItem<L_DATA> *pInsert(L_DATA, DLItem<L_DATA>*);
    L_DATA pDelete(DLItem<L_DATA>*);
  public:
    DLList(void);
    ~DLList();
    unsigned long Size(void) { return number_of_items; }
    int Insert(L_DATA, unsigned long);
    int Delete(unsigned long);
    int fDelete(L_DATA);
    L_DATA Push(L_DATA);
    L_DATA Pop(void);
    L_DATA Get(unsigned long);
    int Enqueue(L_DATA);
    L_DATA Dequeue(void);
    unsigned long Is_In_List(L_DATA);
};

template <class L_DATA>
class DL_Indexed_List : virtual public DLList<L_DATA>
{
  friend class DLList_Iter<L_DATA>;
  private:
   DLItem<L_DATA> *pInsert(L_DATA, DLItem<L_DATA>*);
   L_DATA pDelete(DLItem<L_DATA>*);
   HugeArray<DLItem<L_DATA>*> array;
   unsigned long last_index;
  public:
      DL_Indexed_List(void);
    ~DL_Indexed_List();
    L_DATA Push(L_DATA);
    L_DATA Pop(void);
    L_DATA Get(unsigned long);
};

//#####################################################################################################

template <class L_DATA> class DLList_Iter
{
    private:
     DLList<L_DATA>  *list;
     DLItem<L_DATA> *current;
     bool end_reached;
    public:
      DLList_Iter(void);
      ~DLList_Iter() {end_reached=true;};
      L_DATA Next(void);
      L_DATA Previous(void);
      L_DATA First(DLList<L_DATA> *l);
      L_DATA Last(DLList<L_DATA> *l);
      bool End(void) {return end_reached;}
      DLItem<L_DATA> *Get_Current(void) {return current;}
      L_DATA Get_Current_Item(void) {return current->item;}
      void Set_Current(DLItem<L_DATA> *c) {current=c;}
      void Set_Status(bool s) {end_reached=s;}
      bool Swap(DLList_Iter<L_DATA>);  //swapt die beiden Elemente, wenn sie in der gleichen Liste stehen!!
 
};

//#####################################################################################################
 struct RGBcolor {
    unsigned int red;
    unsigned int green;
    unsigned int blue;
    char pajek_c[20];
};
//-------------------------------------------------------------------------------

class NLink;

class NNode
{
  friend class NLink;
  private :
    unsigned long index;
    unsigned long cluster_index;
    unsigned long marker, affiliations;
    unsigned long *state_history;
    unsigned int max_states;
    long distance;
    double clustering;
    double weight;
    double affinity;
//    double old_weight;
    
    DLList<NNode*> *neighbours;    //list with pointers to neighbours
    DLList<NLink*> *n_links;
    DLList<NLink*> *global_link_list;
    char name[255];
    RGBcolor color;
  public :
    NNode(unsigned long, unsigned long, DLList<NLink*>*, char*, int);
    ~NNode();
    unsigned long Get_Index(void)  { return(index); }
    unsigned long Get_ClusterIndex(void) { return(cluster_index);}
    unsigned long Get_Marker(void) { return marker;}
    void Set_Marker(unsigned long m) {marker=m;}
    unsigned long Get_Affiliations(void) { return affiliations;}
    void Set_Affiliations(unsigned long m) {affiliations=m;}
    void Set_ClusterIndex(unsigned long ci) { cluster_index=ci; return;}
    void Set_Index(unsigned long i) {index=i; return;}
    unsigned long Get_Degree(void) { return(neighbours->Size());}
    char *Get_Name(void) {return name;}
    void Set_Name(char* n) {strcpy(name,n);}
    double Get_Links_Among_Neigbours(void);
    double Get_Clustering(void);
    double Get_Weight(void) {return weight;}
    double Get_Affinity(void) {return affinity;}
    unsigned long *Get_StateHistory(void) {return state_history;}
    void Add_StateHistory(unsigned int q);
  //  double Get_OldWeight(void) {return old_weight;}
    void Set_Weight(double w) {weight=w;}
    void Set_Affinity(double w) {affinity=w;}

  //  void Set_OldWeight(double w) {old_weight=w;}
    long Get_Distance(void) {return distance;}
    void Set_Distance(long d) {distance=d;}
    int  Connect_To(NNode*, double);
    DLList<NNode*> *Get_Neighbours(void) {return neighbours;}
    DLList<NLink*> *Get_Links(void) {return n_links;}
    int  Disconnect_From(NNode*);
    int  Disconnect_From_All(void);
    bool Is_Linked_To(NNode*);
    RGBcolor Get_Color(void) {return color;}
    void Set_Color(RGBcolor c);
    NLink *Get_LinkToNeighbour(NNode *neighbour);
};

//#####################################################################################################

class NLink
{
  friend class NNode;
  private :
    NNode *start;
    NNode *end;
    double weight;
    double old_weight;
    unsigned long index;
    unsigned long marker;
  public :
    NLink( NNode*, NNode*, double);
    ~NLink();
    unsigned long Get_Start_Index(void)  { return(start->Get_Index()); }
    unsigned long Get_End_Index(void)    { return(end->Get_Index()); }
    NNode *Get_Start(void) {return(start);}
    NNode *Get_End(void) {return(end);}
    double Get_Weight(void) {return weight;}
    void Set_Weight(double w) {weight=w;}
    double Get_OldWeight(void) {return old_weight;}
    void Set_OldWeight(double w) {old_weight=w;}
    unsigned long Get_Marker(void) {return marker;}
    void Set_Marker(unsigned long m) {marker=m;}
    unsigned long Get_Index() {return index;}
    void Set_Index(unsigned long i) {index=i;}
};

//#####################################################################################################

template <class L_DATA>  class ClusterList : public DLList<L_DATA>
{
  friend class DLList_Iter<L_DATA>;
  private:
  long links_out_of_cluster;
  unsigned long links_inside_cluster;
  unsigned long frequency;
  double cluster_energy;
  DLList<L_DATA> *candidates;
  long marker;
  public:
  ClusterList(void);
   ~ClusterList();
  long Get_Links_OOC(void) {return(links_out_of_cluster);}
  void Set_Links_OOC(long looc) {links_out_of_cluster=looc;}
  unsigned long Get_Links_IC(void) {return(links_inside_cluster);}
  unsigned long Get_Frequency(void) {return(frequency);}
  void IncreaseFrequency(void) {frequency++;}
  void Set_Links_IC(unsigned long lic) {links_inside_cluster=lic;}
  double Get_Energy(void) {return (cluster_energy);}
  void Set_Energy(double e) {cluster_energy=e;}
  DLList<L_DATA> *Get_Candidates(void) {return candidates;}
  bool operator<(ClusterList<L_DATA> &b);
  bool operator==(ClusterList <L_DATA> &b);
  long Get_Marker(void) {return marker;}
  void Set_Marker(long m) {marker=m;}  
};
//#####################################################################################################
template <class L_DATA>
class DL_Node_List : virtual public DL_Indexed_List<NNode*>
{
  friend class DLList_Iter<L_DATA>;
  private:
   DLItem<L_DATA> *pInsert(NNode*, DLItem<NNode*>*);
   NNode* pDelete(DLItem<NNode*>*);
   HugeArray<DLItem<NNode*>*> array;
   unsigned long last_index;
  public:
    DL_Node_List(void);
    ~DL_Node_List();
    NNode* Push(NNode*);
    NNode* Pop(void);
    NNode* Get(unsigned long);
    int Delete(unsigned long);

};
//#####################################################################################################



struct cluster_join_move
{
  ClusterList<NNode*> *c1;
  ClusterList<NNode*> *c2;
  double joint_energy;
  long joint_looc;
  unsigned long joint_lic;
} ;

struct network
{
  DL_Indexed_List<NNode*> *node_list;
  DL_Indexed_List<NLink*> *link_list;
  DL_Indexed_List<ClusterList<NNode*>*> *cluster_list;
  DL_Indexed_List<cluster_join_move*> *moveset;
  unsigned long max_k;
  unsigned long min_k;
  unsigned long diameter;
  double av_weight;
  double max_weight;
  double min_weight;
  double sum_weights;
  double av_k;
  double av_bids;
  unsigned long max_bids;
  unsigned long min_bids;
  unsigned long sum_bids;
} ;

/*
struct network
{
  DLList<NNode*> *node_list;
  DLList<NLink*> *link_list;
  DLList<ClusterList<NNode*>*> *cluster_list;
  DLList<cluster_join_move*> *moveset;
} ;
*/

void test2(void);

#endif

