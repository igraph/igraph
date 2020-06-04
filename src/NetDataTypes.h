/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

/* The original version of this file was written by Jörg Reichardt
   The original copyright notice follows here */

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

#include <cstring>

//###########################################################################################

struct HUGE_INDEX {
    unsigned int field_index;
    unsigned long in_field_index;
};

template <class DATA> class HugeArray {
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
    unsigned long Size(void) {
        return max_index;
    }
} ;
//###############################################################################################
template <class L_DATA > class DLList;
template <class L_DATA > class DL_Indexed_List;
template <class L_DATA > class ClusterList;
template <class L_DATA > class DLList_Iter;

template <class L_DATA>
class DLItem {
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
public:
    void del() {
        delete item;
    }
};

template <class L_DATA >
class DLList {
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
    unsigned long Size(void) {
        return number_of_items;
    }
    int Insert(L_DATA, unsigned long);
    int Delete(unsigned long);
    int fDelete(L_DATA);
    L_DATA Push(L_DATA);
    L_DATA Pop(void);
    L_DATA Get(unsigned long);
    int Enqueue(L_DATA);
    L_DATA Dequeue(void);
    unsigned long Is_In_List(L_DATA);
    void delete_items();
};

template <class L_DATA>
class DL_Indexed_List : virtual public DLList<L_DATA> {
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

template <class L_DATA> class DLList_Iter {
private:
    DLList<L_DATA>  *list;
    DLItem<L_DATA> *current;
    bool end_reached;
public:
    DLList_Iter(void);
    ~DLList_Iter() {
        end_reached = true;
    };
    L_DATA Next(void);
    L_DATA Previous(void);
    L_DATA First(DLList<L_DATA> *l);
    L_DATA Last(DLList<L_DATA> *l);
    bool End(void) {
        return end_reached;
    }
    DLItem<L_DATA> *Get_Current(void) {
        return current;
    }
    L_DATA Get_Current_Item(void) {
        return current->item;
    }
    void Set_Current(DLItem<L_DATA> *c) {
        current = c;
    }
    void Set_Status(bool s) {
        end_reached = s;
    }
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

class NNode {
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
    unsigned long Get_Index(void)  {
        return (index);
    }
    unsigned long Get_ClusterIndex(void) {
        return (cluster_index);
    }
    unsigned long Get_Marker(void) {
        return marker;
    }
    void Set_Marker(unsigned long m) {
        marker = m;
    }
    unsigned long Get_Affiliations(void) {
        return affiliations;
    }
    void Set_Affiliations(unsigned long m) {
        affiliations = m;
    }
    void Set_ClusterIndex(unsigned long ci) {
        cluster_index = ci;
        return;
    }
    void Set_Index(unsigned long i) {
        index = i;
        return;
    }
    unsigned long Get_Degree(void) {
        return (neighbours->Size());
    }
    char *Get_Name(void) {
        return name;
    }
    void Set_Name(char* n) {
        strcpy(name, n);
    }
    double Get_Links_Among_Neigbours(void);
    double Get_Clustering(void);
    double Get_Weight(void) {
        return weight;
    }
    double Get_Affinity(void) {
        return affinity;
    }
    unsigned long *Get_StateHistory(void) {
        return state_history;
    }
    void Add_StateHistory(unsigned int q);
    //  double Get_OldWeight(void) {return old_weight;}
    void Set_Weight(double w) {
        weight = w;
    }
    void Set_Affinity(double w) {
        affinity = w;
    }

    //  void Set_OldWeight(double w) {old_weight=w;}
    long Get_Distance(void) {
        return distance;
    }
    void Set_Distance(long d) {
        distance = d;
    }
    int  Connect_To(NNode*, double);
    DLList<NNode*> *Get_Neighbours(void) {
        return neighbours;
    }
    DLList<NLink*> *Get_Links(void) {
        return n_links;
    }
    int  Disconnect_From(NNode*);
    int  Disconnect_From_All(void);
    bool Is_Linked_To(NNode*);
    RGBcolor Get_Color(void) {
        return color;
    }
    void Set_Color(RGBcolor c);
    NLink *Get_LinkToNeighbour(NNode *neighbour);
};

//#####################################################################################################

class NLink {
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
    unsigned long Get_Start_Index(void)  {
        return (start->Get_Index());
    }
    unsigned long Get_End_Index(void)    {
        return (end->Get_Index());
    }
    NNode *Get_Start(void) {
        return (start);
    }
    NNode *Get_End(void) {
        return (end);
    }
    double Get_Weight(void) {
        return weight;
    }
    void Set_Weight(double w) {
        weight = w;
    }
    double Get_OldWeight(void) {
        return old_weight;
    }
    void Set_OldWeight(double w) {
        old_weight = w;
    }
    unsigned long Get_Marker(void) {
        return marker;
    }
    void Set_Marker(unsigned long m) {
        marker = m;
    }
    unsigned long Get_Index() {
        return index;
    }
    void Set_Index(unsigned long i) {
        index = i;
    }
};

//#####################################################################################################

template <class L_DATA>  class ClusterList : public DLList<L_DATA> {
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
    long Get_Links_OOC(void) {
        return (links_out_of_cluster);
    }
    void Set_Links_OOC(long looc) {
        links_out_of_cluster = looc;
    }
    unsigned long Get_Links_IC(void) {
        return (links_inside_cluster);
    }
    unsigned long Get_Frequency(void) {
        return (frequency);
    }
    void IncreaseFrequency(void) {
        frequency++;
    }
    void Set_Links_IC(unsigned long lic) {
        links_inside_cluster = lic;
    }
    double Get_Energy(void) {
        return (cluster_energy);
    }
    void Set_Energy(double e) {
        cluster_energy = e;
    }
    DLList<L_DATA> *Get_Candidates(void) {
        return candidates;
    }
    bool operator<(ClusterList<L_DATA> &b);
    bool operator==(ClusterList <L_DATA> &b);
    long Get_Marker(void) {
        return marker;
    }
    void Set_Marker(long m) {
        marker = m;
    }
};
//#####################################################################################################
template <class L_DATA>
class DL_Node_List : virtual public DL_Indexed_List<NNode*> {
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



struct cluster_join_move {
    ClusterList<NNode*> *c1;
    ClusterList<NNode*> *c2;
    double joint_energy;
    long joint_looc;
    unsigned long joint_lic;
} ;

struct network {
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

template <class DATA>
HugeArray<DATA>::HugeArray(void) {
    max_bit_left = 1 << 31; //wir setzen das 31. Bit auf 1
    size = 2;
    max_index = 0;
    highest_field_index = 0;
    data = new DATA[2]; //ein extra Platz fuer das Nullelement
    data[0] = 0;
    data[1] = 0;
    for (int i = 0; i < 32; i++) {
        fields[i] = NULL;
    }
    fields[highest_field_index] = data;
}

template <class DATA> HugeArray<DATA>::~HugeArray(void) {
    for (unsigned int i = 0; i <= highest_field_index; i++) {
        data = fields[i];
        delete [] data;
    }
}

template <class DATA>
HUGE_INDEX HugeArray<DATA>::get_huge_index(unsigned long index) {
    HUGE_INDEX h_index;
    unsigned int shift_index = 0;
    unsigned long help_index;
    help_index = index;
    if (index < 2) {
        h_index.field_index = 0;
        h_index.in_field_index = index;
        return h_index;
    }
    // wie oft muessen wir help_index nach links shiften, damit das 31. Bit gesetzt ist??
    while (!(max_bit_left & help_index)) {
        help_index <<= 1;
        shift_index++;
    }
    h_index.field_index = 31 - shift_index;   // das hoechste  besetzte Bit im Index
    help_index = 1 << h_index.field_index;  // in help_index wird das hoechste besetzte Bit von Index gesetzt
    h_index.in_field_index = (index ^ help_index); // index XOR help_index, womit alle bits unter dem hoechsten erhalten bleiben
    return h_index;
}

template <class DATA>
DATA &HugeArray<DATA>::Set(unsigned long int index) {
    HUGE_INDEX h_index;
    unsigned long data_size;
    while (size < index + 1) {
        highest_field_index++;
        data_size = 1 << highest_field_index;
        data = new DATA[data_size];
        for (unsigned long i = 0; i < data_size; i++) {
            data[i] = 0;
        }
        size = size + data_size; //overflow noch abfangen
        //printf("Vergroesserung auf: %u bei index %u\n",size,index);
        fields[highest_field_index] = data;
    }
    h_index = get_huge_index(index);
//printf("index %lu = %lu . %lu\n",index,h_index.field_index,h_index.in_field_index);
    data = fields[h_index.field_index];
    if (max_index < index) {
        max_index = index;
    }
    return (data[h_index.in_field_index]);
}

template <class DATA>
DATA HugeArray<DATA>::Get(unsigned long index) {
    return (Set(index));
}


template <class DATA>
DATA &HugeArray<DATA>::operator[](unsigned long index) {
    return (Set(index));
}


//###############################################################################
template <class L_DATA>
DLItem<L_DATA>::DLItem(L_DATA i, unsigned long ind) : item(i), index(ind), previous(0), next(0) {
}

template <class L_DATA>
DLItem<L_DATA>::DLItem(L_DATA i, unsigned long ind, DLItem<L_DATA> *p, DLItem<L_DATA> *n) : item(i), index(ind), previous(p), next(n) {
}

template <class L_DATA>
DLItem<L_DATA>::~DLItem() {
//delete item;      //eigentlich muessten wir pruefen, ob item ueberhaupt ein Pointer ist...
//previous=NULL;
//next=NULL;
}


//######################################################################################################################
template <class L_DATA>
DLList<L_DATA>::DLList(void) {
    head = tail = NULL;
    number_of_items = 0;
    head = new DLItem<L_DATA>(NULL, 0); //fuer head und Tail gibt es das gleiche Array-Element!! Vorsicht!!
    tail = new DLItem<L_DATA>(NULL, 0);
    if ( !head || !tail ) {
        if (head) {
            delete (head);
        }
        if (tail) {
            delete (tail);
        }
        return;
    }  else {
        head->next = tail;
        tail->previous = head;
    }
}

template <class L_DATA>
DLList<L_DATA>::~DLList() {
    DLItem<L_DATA> *cur = head, *next;
    while (cur) {
        next = cur->next;
        delete (cur);
        cur = next;
    }
    number_of_items = 0;
    //  printf("Liste Zerstoert!\n");
}

template <class L_DATA>
void DLList<L_DATA>::delete_items() {
    DLItem<L_DATA> *cur, *next;
    cur = this->head;
    while (cur) {
        next = cur->next;
        cur->del();
        cur = next;
    }
    this->number_of_items = 0;
}

//privates Insert
template <class L_DATA>
DLItem<L_DATA> *DLList<L_DATA>::pInsert(L_DATA data, DLItem<L_DATA> *pos) {
    DLItem<L_DATA> *i = new DLItem<L_DATA>(data, number_of_items + 1, pos->previous, pos);
    if (i) {
        pos->previous->next = i;
        pos->previous = i;
        number_of_items++;
        return (i);
    } else {
        return (0);
    }
}
//privates delete
template <class L_DATA>
L_DATA DLList<L_DATA>::pDelete(DLItem<L_DATA> *i) {
    L_DATA data = i->item;
    i->previous->next = i->next;
    i->next->previous = i->previous;
//  array[i->index]=0;
    delete (i);
    number_of_items--;
    return (data);
}
//oeffentliches Insert
template <class L_DATA>
int DLList<L_DATA>::Insert(L_DATA data, unsigned long pos) {
    if ((pos < 0) || (pos > (number_of_items))) {
        return (0);
    }
    DLItem<L_DATA> *cur = head;
    while (pos--) {
        cur = cur->next;
    }
    return (pInsert(data, cur) != 0);
}
//oeffentliche Delete
template <class L_DATA>
int DLList<L_DATA>::Delete(unsigned long pos) {
    if ((pos < 0) || (pos > (number_of_items))) {
        return (0);
    }
    DLItem<L_DATA> *cur = head;
    while (pos--) {
        cur = cur->next;
    }
    return (pDelete(cur) != 0);
}

//oeffentliche Delete
template <class L_DATA>
int DLList<L_DATA>::fDelete(L_DATA data) {
    if ((number_of_items == 0) || (!data)) {
        return (0);
    }
    DLItem<L_DATA> *cur;
    cur = head->next;
    while ((cur != tail) && (cur->item != data)) {
        cur = cur->next;
    }
    if (cur != tail) {
        return (pDelete(cur) != 0);
    }
    return (0);
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Push(L_DATA data) {
    DLItem<L_DATA> *tmp;
    tmp = pInsert(data, tail);
    if (tmp) {
        return (tmp->item);
    }
    return (0);
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Pop(void) {
    return (pDelete(tail->previous));
}


template <class L_DATA>
L_DATA DLList<L_DATA>::Get(unsigned long pos) {
    if ((pos < 1) || (pos > (number_of_items + 1))) {
        return (0);
    }
//  return(array[pos]->item);
    DLItem<L_DATA> *cur = head;
    while (pos--) {
        cur = cur->next;
    }
    return (cur->item);
}


template <class L_DATA>
int DLList<L_DATA>::Enqueue(L_DATA data) {
    return (pInsert(data, tail) != 0);
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Dequeue(void) {
    return (pDelete(head->next));
}

//gibt Index des gesuchte Listenelement zurueck, besser waere eigentlich zeiger
template <class L_DATA>
unsigned long DLList<L_DATA>::Is_In_List(L_DATA data) {
    DLItem<L_DATA> *cur = head, *next;
    unsigned long pos = 0;
    while (cur) {
        next = cur->next;
        if (cur->item == data) {
            return (pos) ;
        }
        cur = next;
        pos++;
    }
    return (0);
}

//######################################################################################################################
template <class L_DATA>
DL_Indexed_List<L_DATA>::DL_Indexed_List(void) : DLList<L_DATA>() {
    last_index = 0;
}

template <class L_DATA>
DL_Indexed_List<L_DATA>::~DL_Indexed_List() {
    /* This is already done by the DLList destructor */
    /*   DLItem<L_DATA> *cur, *next; */
    /*   cur=this->head; */
    /*   while (cur) */
    /*     { */
    /*       next=cur->next; */
    /*       delete(cur); */
    /*       cur=next; */
    /*     } */
    /*     this->number_of_items=0; */
    //  printf("Liste Zerstoert!\n");
}

//privates Insert
template <class L_DATA>
DLItem<L_DATA> *DL_Indexed_List<L_DATA>::pInsert(L_DATA data, DLItem<L_DATA> *pos) {
    DLItem<L_DATA> *i = new DLItem<L_DATA>(data, last_index, pos->previous, pos);
    if (i) {
        pos->previous->next = i;
        pos->previous = i;
        this->number_of_items++;
        array[last_index] = i;
        last_index++;
        return (i);
    } else {
        return (0);
    }
}
//privates delete
template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::pDelete(DLItem<L_DATA> *i) {
    L_DATA data = i->item;
    i->previous->next = i->next;
    i->next->previous = i->previous;
    array[i->index] = 0;
    last_index = i->index;
    delete (i);
    this->number_of_items--;
    return (data);
}
template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Push(L_DATA data) {
    DLItem<L_DATA> *tmp;
    tmp = pInsert(data, this->tail);
    if (tmp) {
        return (tmp->item);
    }
    return (0);
}

template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Pop(void) {
    return (pDelete(this->tail->previous));
}

template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Get(unsigned long pos) {
    if (pos > this->number_of_items - 1) {
        return (0);
    }
    return (array[pos]->item);
}

//#######################################################################################

//************************************************************************************************************
template <class L_DATA>
ClusterList<L_DATA>::ClusterList(void) : DLList<L_DATA>() {
    links_out_of_cluster = 0;
    links_inside_cluster = 0;
    frequency = 1;
    cluster_energy = 1e30;
    candidates = new DLList<L_DATA>();
    marker = 0;
}

template <class L_DATA>
ClusterList<L_DATA>::~ClusterList() {
    while (candidates->Size()) {
        candidates->Pop();
    }
    delete candidates;
}


template <class L_DATA>
bool ClusterList<L_DATA>::operator==(ClusterList<L_DATA> &b) {
    bool found = false;
    L_DATA n_cur, n_cur_b;
    DLList_Iter<L_DATA> a_iter, b_iter;

    if (this->Size() != b.Size()) {
        return false;
    }

    n_cur = a_iter.First(this);
    while (!(a_iter.End())) {
        found = false;
        n_cur_b = b_iter.First(&b);
        while (!(b_iter.End()) && !found) {
            if (n_cur == n_cur_b) {
                found = true;
            }
            n_cur_b = b_iter.Next();
        }
        if (!found) {
            return false;
        }
        n_cur = a_iter.Next();
    }
    return (found);
}
//A<B ist Wahr, wenn A echte Teilmenge von B ist
template <class L_DATA>
bool ClusterList<L_DATA>::operator<(ClusterList<L_DATA> &b) {
    bool found = false;
    L_DATA n_cur, n_cur_b;
    DLList_Iter<L_DATA> a_iter, b_iter;

    if (this->Size() >= b.Size()) {
        return false;
    }
    n_cur = a_iter.First(this);
    while (!(a_iter.End())) {
        found = false;
        n_cur_b = b_iter.First(&b);
        while (!(b_iter.End()) && !found) {
            if (n_cur == n_cur_b) {
                found = true;
            }
            n_cur_b = b_iter.Next();
        }
        if (!found) {
            return false;
        }
        n_cur = a_iter.Next();
    }
    return (found);
}

//#####################################################################################
template <class L_DATA>
DLList_Iter<L_DATA>::DLList_Iter() {
    list = NULL;
    current = NULL;
    end_reached = true;
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::Next(void) {
    current = current->next;
    if (current == (list->tail)) {
        end_reached = true;
    }
    return (current->item);
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::Previous(void) {
    current = current->previous;
    if (current == (list->head)) {
        end_reached = true;
    }
    return (current->item);
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::First(DLList<L_DATA> *l) {
    list = l;
    current = list->head->next;
    if (current == (list->tail)) {
        end_reached = true;
    } else {
        end_reached = false;
    }
    return (current->item);
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::Last(DLList<L_DATA> *l) {
    list = l;
    current = list->tail->previous;
    if (current == (list->head)) {
        end_reached = true;    // falls die List leer ist
    } else {
        end_reached = false;
    }
    return (current->item);
}

template <class L_DATA>
bool DLList_Iter<L_DATA>::Swap(DLList_Iter<L_DATA> b) {
    L_DATA h;
    if (list != b.list) {
        return false;    //elemeten muessen aus der gleichen List stammen
    }
    if (end_reached || b.end_reached) {
        return false;
    }
    h = current->item; current->item = b.current->item; b.current->item = h;
    return true;
}

#endif

