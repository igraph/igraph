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

#include "igraph_types.h"

#include <cassert>
#include <cstring>

// In igraph, we set node names to be a string representation of the one-based
// vertex ID. This takes at most 20 characters. Add one for a potential sign
// (should not happen) and one more for the null terminator.
#define SPINGLASS_MAX_NAME_LEN 22

//###########################################################################################

struct HUGE_INDEX {
    unsigned int field_index;
    igraph_integer_t in_field_index;
};

template <class DATA>
class HugeArray {
    igraph_integer_t size = 2;
    unsigned int highest_field_index = 0;
    const igraph_integer_t max_bit_left = 1UL << 31; //wir setzen das 31. Bit auf 1
    igraph_integer_t max_index = 0;
    DATA *data;
    DATA *fields[32];
public:
    HugeArray();
    HugeArray(const HugeArray &) = delete;
    HugeArray & operator = (const HugeArray &) = delete;
    ~HugeArray();
    HUGE_INDEX get_huge_index(igraph_integer_t) const;
    DATA &Set(igraph_integer_t index);
    DATA Get(igraph_integer_t index) { return Set(index); }
    DATA &operator[](igraph_integer_t index) { return Set(index); }
    igraph_integer_t Size() const { return max_index; }
} ;

//###############################################################################################
template <class L_DATA > class DLList;
template <class L_DATA > class DL_Indexed_List;
template <class L_DATA > using ClusterList= DLList<L_DATA>;
template <class L_DATA > class DLList_Iter;

template <class L_DATA>
class DLItem {
    friend class DLList<L_DATA> ;
    friend class DL_Indexed_List<L_DATA>;
    friend class DLList_Iter<L_DATA>;

    L_DATA  item;
    igraph_integer_t index;
    DLItem *previous;
    DLItem *next;
    DLItem(L_DATA i, igraph_integer_t ind);
    DLItem(L_DATA i, igraph_integer_t ind, DLItem<L_DATA> *p, DLItem<L_DATA> *n);
public:
    void del() {
        delete item;
    }
};

template <class L_DATA >
class DLList {
    friend class DLList_Iter<L_DATA>;
protected:
    DLItem<L_DATA> *head;
    DLItem<L_DATA> *tail;
    igraph_integer_t number_of_items = 0;
    virtual DLItem<L_DATA> *pInsert(L_DATA, DLItem<L_DATA>*);
    virtual L_DATA pDelete(DLItem<L_DATA>*);
public:
    DLList();
    DLList(const DLList &) = delete;
    DLList & operator = (const DLList &) = delete;
    virtual ~DLList();
    igraph_integer_t Size() const {
        return number_of_items;
    }
    int fDelete(L_DATA);
    virtual L_DATA Push(L_DATA);
    virtual L_DATA Pop();
    virtual L_DATA Get(igraph_integer_t);
    igraph_integer_t Is_In_List(L_DATA);
    void delete_items();
};

template <class L_DATA>
class DL_Indexed_List : public DLList<L_DATA> {
    DLItem<L_DATA> *pInsert(L_DATA, DLItem<L_DATA>*) final;
    L_DATA pDelete(DLItem<L_DATA>*) final;
    HugeArray<DLItem<L_DATA>*> array;
    igraph_integer_t last_index = 0;

public:
    DL_Indexed_List() = default;
    L_DATA Push(L_DATA) final;
    L_DATA Pop() final;
    L_DATA Get(igraph_integer_t) final;
};

//#####################################################################################################

template <class L_DATA> class DLList_Iter {
    const DLList<L_DATA> *list = nullptr;
    const DLItem<L_DATA> *current = nullptr;
    bool end_reached = true;

public:
    L_DATA Next();
    L_DATA Previous();
    L_DATA First(const DLList<L_DATA> *l);
    L_DATA Last(const DLList<L_DATA> *l);
    bool End() const {
        return end_reached;
    }
    bool Swap(DLList_Iter<L_DATA>);  //swapt die beiden Elemente, wenn sie in der gleichen Liste stehen!!
};

//#####################################################################################################

class NLink;

class NNode {
    igraph_integer_t index;
    igraph_integer_t cluster_index;
    igraph_integer_t marker = 0;
    double weight = 0.0;

    DLList<NNode*> neighbours;    //list with pointers to neighbours
    DLList<NLink*> n_links;
    DLList<NLink*> *global_link_list;
    char name[SPINGLASS_MAX_NAME_LEN];
public :
    NNode(igraph_integer_t ind, igraph_integer_t c_ind, DLList<NLink *> *ll, const char *n) :
            index(ind), cluster_index(c_ind), global_link_list(ll)
    {
        strcpy(name, n);
    }
    NNode(const NNode &) = delete;
    NNode &operator=(const NNode &) = delete;
    ~NNode() { Disconnect_From_All(); }

    igraph_integer_t Get_Index() const {
        return index;
    }
    igraph_integer_t Get_ClusterIndex() const {
        return cluster_index;
    }
    igraph_integer_t Get_Marker() const {
        return marker;
    }
    void Set_Marker(igraph_integer_t m) {
        marker = m;
    }
    void Set_ClusterIndex(igraph_integer_t ci) {
        cluster_index = ci;
    }
    igraph_integer_t Get_Degree() const {
        return (neighbours.Size());
    }
    const char *Get_Name() {
        return name;
    }
    void Set_Name(const char *n) {
        strcpy(name, n);
    }
    double Get_Weight() const {
        return weight;
    }
    void Set_Weight(double w) {
        weight = w;
    }
    int Connect_To(NNode*, double);
    const DLList<NNode*> *Get_Neighbours() const {
        return &neighbours;
    }
    const DLList<NLink*> *Get_Links() const {
        return &n_links;
    }
    igraph_integer_t Disconnect_From(NNode*);
    igraph_integer_t Disconnect_From_All();
    NLink *Get_LinkToNeighbour(const NNode *neighbour);
};

//#####################################################################################################

class NLink {

    NNode *start;
    NNode *end;
    double weight;

public :
    NLink(NNode *s, NNode *e, double w) : start(s), end(e), weight(w) { }
    NLink(const NLink &) = delete;
    NLink & operator = (const NLink &) = delete;
    ~NLink() { start->Disconnect_From(end); }

    NNode *Get_Start() { return start; }
    NNode *Get_End() { return end; }
    const NNode *Get_Start() const { return start; }
    const NNode *Get_End() const { return end; }

    double Get_Weight() const { return weight; }
};

//#####################################################################################################

struct network {
    DL_Indexed_List<NNode*> node_list;
    DL_Indexed_List<NLink*> link_list;
    DL_Indexed_List<ClusterList<NNode*>*> cluster_list;
    double sum_weights;

    network() = default;
    network (const network &) = delete;
    network & operator = (const network &) = delete;

    ~network() {
        ClusterList<NNode*> *cl_cur;

        while (link_list.Size()) {
            delete link_list.Pop();
        }
        while (node_list.Size()) {
            delete node_list.Pop();
        }
        while (cluster_list.Size()) {
            cl_cur = cluster_list.Pop();
            while (cl_cur->Size()) {
                cl_cur->Pop();
            }
            delete cl_cur;
        }
    }
};

template <class DATA>
HugeArray<DATA>::HugeArray() {
    data = new DATA[2]; //ein extra Platz fuer das Nullelement
    data[0] = 0;
    data[1] = 0;
    for (auto & field : fields) {
        field = nullptr;
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
HUGE_INDEX HugeArray<DATA>::get_huge_index(igraph_integer_t index) const {
    HUGE_INDEX h_index;
    unsigned int shift_index = 0;
    igraph_integer_t help_index;
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
    help_index = igraph_integer_t(1) << h_index.field_index;  // in help_index wird das hoechste besetzte Bit von Index gesetzt
    h_index.in_field_index = (index ^ help_index); // index XOR help_index, womit alle bits unter dem hoechsten erhalten bleiben
    return h_index;
}

template <class DATA>
DATA &HugeArray<DATA>::Set(igraph_integer_t index) {
    igraph_integer_t data_size;
    while (size < index + 1) {
        highest_field_index++;
        data_size = 1UL << highest_field_index;
        data = new DATA[data_size];
        for (igraph_integer_t i = 0; i < data_size; i++) {
            data[i] = 0;
        }
        size = size + data_size; //overflow noch abfangen
        fields[highest_field_index] = data;
    }
    HUGE_INDEX h_index = get_huge_index(index);
    data = fields[h_index.field_index];
    if (max_index < index) {
        max_index = index;
    }
    return data[h_index.in_field_index];
}


//###############################################################################
template <class L_DATA>
DLItem<L_DATA>::DLItem(L_DATA i, igraph_integer_t ind) :
    item(i), index(ind), previous(nullptr), next(nullptr) { }

template <class L_DATA>
DLItem<L_DATA>::DLItem(L_DATA i, igraph_integer_t ind, DLItem<L_DATA> *p, DLItem<L_DATA> *n) :
    item(i), index(ind), previous(p), next(n) { }

//######################################################################################################################
template <class L_DATA>
DLList<L_DATA>::DLList() {
    head = new DLItem<L_DATA>(NULL, 0); //fuer head und Tail gibt es das gleiche Array-Element!! Vorsicht!!
    tail = new DLItem<L_DATA>(NULL, 0);

    head->next = tail;
    tail->previous = head;
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
    auto *i = new DLItem<L_DATA>(data, number_of_items + 1, pos->previous, pos);
    pos->previous->next = i;
    pos->previous = i;
    number_of_items++;
    return i;
}
//privates delete
template <class L_DATA>
L_DATA DLList<L_DATA>::pDelete(DLItem<L_DATA> *i) {
    assert(number_of_items > 0);
    L_DATA data = i->item;
    i->previous->next = i->next;
    i->next->previous = i->previous;
    delete i;
    number_of_items--;
    return data;
}

//oeffentliche Delete
template <class L_DATA>
int DLList<L_DATA>::fDelete(L_DATA data) {
    if ((number_of_items == 0) || (!data)) {
        return 0;
    }
    DLItem<L_DATA> *cur;
    cur = head->next;
    while ((cur != tail) && (cur->item != data)) {
        cur = cur->next;
    }
    if (cur != tail) {
        return (pDelete(cur) != 0);
    }
    return 0;
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Push(L_DATA data) {
    DLItem<L_DATA> *tmp = pInsert(data, tail);
    return tmp->item;
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Pop() {
    return pDelete(tail->previous);
}


template <class L_DATA>
L_DATA DLList<L_DATA>::Get(igraph_integer_t pos) {
    if ((pos < 1) || (pos > (number_of_items + 1))) {
        return 0;
    }
    DLItem<L_DATA> *cur = head;
    while (pos--) {
        cur = cur->next;
    }
    return (cur->item);
}

//gibt Index des gesuchte Listenelement zurueck, besser waere eigentlich zeiger
template <class L_DATA>
igraph_integer_t DLList<L_DATA>::Is_In_List(L_DATA data) {
    DLItem<L_DATA> *cur = head, *next;
    igraph_integer_t pos = 0;
    while (cur) {
        next = cur->next;
        if (cur->item == data) {
            return pos ;
        }
        cur = next;
        pos++;
    }
    return 0;
}

//######################################################################################################################

//privates Insert
template <class L_DATA>
DLItem<L_DATA> *DL_Indexed_List<L_DATA>::pInsert(L_DATA data, DLItem<L_DATA> *pos) {
    auto *i = new DLItem<L_DATA>(data, last_index, pos->previous, pos);
    pos->previous->next = i;
    pos->previous = i;
    this->number_of_items++;
    array[last_index] = i;
    last_index++;
    return i;
}
//privates delete
template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::pDelete(DLItem<L_DATA> *i) {
    assert(this->number_of_items > 0);
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
    return tmp->item;
}

template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Pop() {
    return pDelete(this->tail->previous);
}

template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Get(igraph_integer_t pos) {
    if (pos > this->number_of_items - 1) {
        return 0;
    }
    return array[pos]->item;
}

//#####################################################################################

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
L_DATA DLList_Iter<L_DATA>::First(const DLList<L_DATA> *l) {
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
L_DATA DLList_Iter<L_DATA>::Last(const DLList<L_DATA> *l) {
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
