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
    HUGE_INDEX get_huge_index(unsigned long) const;
    DATA &Set(unsigned long);
    DATA Get(unsigned long);
    HugeArray();
    ~HugeArray();
    DATA &operator[](unsigned long);
    unsigned long Size() const {
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
    DLList();
    ~DLList();
    unsigned long Size() const {
        return number_of_items;
    }
    int Insert(L_DATA, unsigned long);
    int Delete(unsigned long);
    int fDelete(L_DATA);
    L_DATA Push(L_DATA);
    L_DATA Pop();
    L_DATA Get(unsigned long);
    int Enqueue(L_DATA);
    L_DATA Dequeue();
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
    DL_Indexed_List();
    L_DATA Push(L_DATA);
    L_DATA Pop();
    L_DATA Get(unsigned long);
};

//#####################################################################################################

template <class L_DATA> class DLList_Iter {
private:
    DLList<L_DATA>  *list;
    DLItem<L_DATA> *current;
    bool end_reached;
public:
    DLList_Iter();
    ~DLList_Iter() {
        end_reached = true;
    };
    L_DATA Next();
    L_DATA Previous();
    L_DATA First(DLList<L_DATA> *l);
    L_DATA Last(DLList<L_DATA> *l);
    bool End() const {
        return end_reached;
    }
    DLItem<L_DATA> *Get_Current() {
        return current;
    }
    L_DATA Get_Current_Item() {
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

class NLink;

class NNode {
    friend class NLink;
private :
    unsigned long index;
    unsigned long cluster_index;
    unsigned long marker = 0;
    double clustering = 0.0;
    double weight = 0.0;
    double affinity = 0.0;

    DLList<NNode*> *neighbours;    //list with pointers to neighbours
    DLList<NLink*> *n_links;
    DLList<NLink*> *global_link_list;
    char name[255];
public :
    NNode(unsigned long, unsigned long, DLList<NLink *> *, const char *);
    ~NNode();
    unsigned long Get_Index() const {
        return index;
    }
    unsigned long Get_ClusterIndex() const {
        return cluster_index;
    }
    unsigned long Get_Marker() const {
        return marker;
    }
    void Set_Marker(unsigned long m) {
        marker = m;
    }
    void Set_ClusterIndex(unsigned long ci) {
        cluster_index = ci;
    }
    void Set_Index(unsigned long i) {
        index = i;
    }
    unsigned long Get_Degree() const {
        return (neighbours->Size());
    }
    char *Get_Name() {
        return name;
    }
    void Set_Name(const char *n) {
        strcpy(name, n);
    }
    double Get_Links_Among_Neigbours();
    double Get_Clustering();
    double Get_Weight() const {
        return weight;
    }
    double Get_Affinity() const {
        return affinity;
    }

    void Set_Weight(double w) {
        weight = w;
    }
    void Set_Affinity(double w) {
        affinity = w;
    }
    int  Connect_To(NNode*, double);
    DLList<NNode*> *Get_Neighbours() {
        return neighbours;
    }
    DLList<NLink*> *Get_Links() {
        return n_links;
    }
    int  Disconnect_From(NNode*);
    int  Disconnect_From_All();
    bool Is_Linked_To(NNode*);
    NLink *Get_LinkToNeighbour(NNode *neighbour);
};

//#####################################################################################################

class NLink {
    friend class NNode;
private :
    NNode *start;
    NNode *end;
    double weight;
    unsigned long index;
    unsigned long marker;
public :
    NLink( NNode*, NNode*, double);
    ~NLink();
    unsigned long Get_Start_Index() const {
        return start->Get_Index();
    }
    unsigned long Get_End_Index() const {
        return end->Get_Index();
    }
    NNode *Get_Start() {
        return start;
    }
    NNode *Get_End() {
        return end;
    }
    double Get_Weight() const {
        return weight;
    }
    void Set_Weight(double w) {
        weight = w;
    }
};

//#####################################################################################################

template <class L_DATA>  class ClusterList : public DLList<L_DATA> {
    friend class DLList_Iter<L_DATA>;
private:
    long links_out_of_cluster = 0;
    unsigned long links_inside_cluster = 0;
    unsigned long frequency = 1;
    double cluster_energy = 1e30;
    DLList<L_DATA> *candidates;
    long marker = 0;
public:
    ClusterList();
    ~ClusterList();
    long Get_Links_OOC() {
        return links_out_of_cluster;
    }
    void Set_Links_OOC(long looc) {
        links_out_of_cluster = looc;
    }
    unsigned long Get_Links_IC() {
        return links_inside_cluster;
    }
    unsigned long Get_Frequency() {
        return frequency;
    }
    void IncreaseFrequency() {
        frequency++;
    }
    void Set_Links_IC(unsigned long lic) {
        links_inside_cluster = lic;
    }
    double Get_Energy() {
        return cluster_energy;
    }
    void Set_Energy(double e) {
        cluster_energy = e;
    }
    DLList<L_DATA> *Get_Candidates() {
        return candidates;
    }
    bool operator<(ClusterList<L_DATA> &b);
    bool operator==(ClusterList <L_DATA> &b);

};
//#####################################################################################################
template <class L_DATA>
class DL_Node_List : virtual public DL_Indexed_List<NNode*> {
    friend class DLList_Iter<L_DATA>;
private:
    DLItem<L_DATA> *pInsert(NNode*, DLItem<NNode*>*);
    NNode* pDelete(DLItem<NNode*>*);
    HugeArray<DLItem<NNode*>*> array;
public:
    DL_Node_List();
    ~DL_Node_List();
    NNode* Push(NNode*);
    NNode* Pop();
    NNode* Get(unsigned long);
    int Delete(unsigned long);

};
//#####################################################################################################


struct network {
    DL_Indexed_List<NNode*> *node_list;
    DL_Indexed_List<NLink*> *link_list;
    DL_Indexed_List<ClusterList<NNode*>*> *cluster_list;
    unsigned long diameter;
    double sum_weights;

    network() {
        node_list   = new DL_Indexed_List<NNode*>();
        link_list   = new DL_Indexed_List<NLink*>();
        cluster_list = new DL_Indexed_List<ClusterList<NNode*>*>();
    }

    ~network() {
        ClusterList<NNode*> *cl_cur;

        while (link_list->Size()) {
            delete link_list->Pop();
        }
        while (node_list->Size()) {
            delete node_list->Pop();
        }
        while (cluster_list->Size()) {
            cl_cur = cluster_list->Pop();
            while (cl_cur->Size()) {
                cl_cur->Pop();
            }
            delete cl_cur;
        }
        delete link_list;
        delete node_list;
        delete cluster_list;
    }
};

template <class DATA>
HugeArray<DATA>::HugeArray() {
    max_bit_left = 1UL << 31; //wir setzen das 31. Bit auf 1
    size = 2;
    max_index = 0;
    highest_field_index = 0;
    data = new DATA[2]; //ein extra Platz fuer das Nullelement
    data[0] = 0;
    data[1] = 0;
    for (int i = 0; i < 32; i++) {
        fields[i] = nullptr;
    }
    fields[highest_field_index] = data;
}

template <class DATA> HugeArray<DATA>::~HugeArray() {
    for (unsigned int i = 0; i <= highest_field_index; i++) {
        data = fields[i];
        delete [] data;
    }
}

template <class DATA>
HUGE_INDEX HugeArray<DATA>::get_huge_index(unsigned long index) const {
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
    help_index = 1UL << h_index.field_index;  // in help_index wird das hoechste besetzte Bit von Index gesetzt
    h_index.in_field_index = (index ^ help_index); // index XOR help_index, womit alle bits unter dem hoechsten erhalten bleiben
    return h_index;
}

template <class DATA>
DATA &HugeArray<DATA>::Set(unsigned long int index) {
    HUGE_INDEX h_index;
    unsigned long data_size;
    while (size < index + 1) {
        highest_field_index++;
        data_size = 1UL << highest_field_index;
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
DLItem<L_DATA>::DLItem(L_DATA i, unsigned long ind) : item(i), index(ind), previous(nullptr), next(nullptr) {
}

template <class L_DATA>
DLItem<L_DATA>::DLItem(L_DATA i, unsigned long ind, DLItem<L_DATA> *p, DLItem<L_DATA> *n) : item(i), index(ind), previous(p), next(n) {
}

//######################################################################################################################
template <class L_DATA>
DLList<L_DATA>::DLList() {
    head = tail = nullptr;
    number_of_items = 0;
    head = new DLItem<L_DATA>(NULL, 0); //fuer head und Tail gibt es das gleiche Array-Element!! Vorsicht!!
    tail = new DLItem<L_DATA>(NULL, 0);
    if ( !head || !tail ) {
        if (head) {
            delete head;
        }
        if (tail) {
            delete tail;
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
        delete cur;
        cur = next;
    }
    number_of_items = 0;
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
    delete i;
    number_of_items--;
    return data;
}
//oeffentliches Insert
template <class L_DATA>
int DLList<L_DATA>::Insert(L_DATA data, unsigned long pos) {
    if ((pos < 0) || (pos > number_of_items)) {
        return 0;
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
L_DATA DLList<L_DATA>::Pop() {
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
L_DATA DLList<L_DATA>::Dequeue() {
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
            return pos ;
        }
        cur = next;
        pos++;
    }
    return (0);
}

//######################################################################################################################
template <class L_DATA>
DL_Indexed_List<L_DATA>::DL_Indexed_List() : DLList<L_DATA>() {
    last_index = 0;
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
        return i;
    } else {
        return 0;
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
    delete i;
    this->number_of_items--;
    return data;
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
L_DATA DL_Indexed_List<L_DATA>::Pop() {
    return pDelete(this->tail->previous);
}

template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Get(unsigned long pos) {
    if (pos > this->number_of_items - 1) {
        return 0;
    }
    return (array[pos]->item);
}

//#######################################################################################

//************************************************************************************************************
template <class L_DATA>
ClusterList<L_DATA>::ClusterList() : DLList<L_DATA>() {
    candidates = new DLList<L_DATA>();
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
L_DATA DLList_Iter<L_DATA>::Next() {
    current = current->next;
    if (current == (list->tail)) {
        end_reached = true;
    }
    return (current->item);
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::Previous() {
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
