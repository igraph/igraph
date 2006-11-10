/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
                          NetDataTypes.cpp  -  description
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "NetDataTypes.h"

//#################################################################################
template <class DATA>
HugeArray<DATA>::HugeArray(void)
{
 max_bit_left=1<<31;   //wir setzen das 31. Bit auf 1
 size=2;
 max_index=0;
 highest_field_index=0;
 data=new DATA[2]; //ein extra Platz fuer das Nullelement
 data[0]=0;
 data[1]=0; 
 for (int i=0; i<32; i++) fields[i]=NULL;
 fields[highest_field_index]=data;
}

template <class DATA> HugeArray<DATA>::~HugeArray(void)
{
 for (unsigned int i=0; i<=highest_field_index; i++)
 {
   data=fields[i];
   delete(data);
 }
}

template <class DATA>
HUGE_INDEX HugeArray<DATA>::get_huge_index(unsigned long index)
{
  HUGE_INDEX h_index;
  unsigned int shift_index=0;
  unsigned long help_index;
  help_index=index;
  if (index<2) {
    h_index.field_index=0;
    h_index.in_field_index=index;
    return h_index;
  }
  // wie oft muessen wir help_index nach links shiften, damit das 31. Bit gesetzt ist??
  while (!(max_bit_left & help_index))
  {
    help_index <<= 1;
    shift_index++;
  }
  h_index.field_index=31-shift_index;       // das hoechste  besetzte Bit im Index
  help_index=1 << h_index.field_index;    // in help_index wird das hoechste besetzte Bit von Index gesetzt
  h_index.in_field_index=(index ^ help_index);   // index XOR help_index, womit alle bits unter dem hoechsten erhalten bleiben
  return h_index;
}

template <class DATA>
DATA &HugeArray<DATA>::Set(unsigned long int index)
{
 HUGE_INDEX h_index;
 unsigned long data_size;
 while (size<index+1)
 {
    highest_field_index++;
    data_size=1<<highest_field_index;
    data=new DATA[data_size];
    for (unsigned long i=0; i<data_size; i++)
    {
      data[i]=0;
    }
    size=size+data_size;   //overflow noch abfangen
    //printf("Vergroesserung auf: %u bei index %u\n",size,index);
    fields[highest_field_index]=data;
 }
 h_index=get_huge_index(index);
 //printf("index %lu = %lu . %lu\n",index,h_index.field_index,h_index.in_field_index);
 data=fields[h_index.field_index];
 if (max_index<index) {max_index=index;}
 return(data[h_index.in_field_index]);
}

template <class DATA>
DATA HugeArray<DATA>::Get(unsigned long index)
{
  return(Set(index));
}


template <class DATA>
DATA &HugeArray<DATA>::operator[](unsigned long index)
{
  return(Set(index));
}

                                                                                                                                                                                                                                                                                                                                                                  
//###############################################################################
template <class L_DATA>
DLItem<L_DATA>::DLItem(L_DATA i, unsigned long ind) : item(i), index(ind), previous(0), next(0)
{
}

template <class L_DATA>
DLItem<L_DATA>::DLItem(L_DATA i, unsigned long ind, DLItem<L_DATA> *p, DLItem<L_DATA> *n) : item(i), index(ind), previous(p), next(n)
{
}

template <class L_DATA>
DLItem<L_DATA>::~DLItem()
{
//  delete item;      //eigentlich muessten wir pruefen, ob item ueberhaupt ein Pointer ist...
 //previous=NULL;
 //next=NULL;
}


//######################################################################################################################
template <class L_DATA>
DLList<L_DATA>::DLList(void)
{
   head=tail=NULL;
   number_of_items=0;
   head=new DLItem<L_DATA>(NULL,0);  //fuer head und Tail gibt es das gleiche Array-Element!! Vorsicht!!
   tail=new DLItem<L_DATA>(NULL,0);
   if ( !head || !tail )
     {
      if (head) delete(head);
      if (tail) delete(tail);
      return;
     }  else {
     head->next=tail;
     tail->previous=head;
   }
}

template <class L_DATA>
DLList<L_DATA>::~DLList()
{
  DLItem<L_DATA> *cur=head, *next;
  while (cur)
    {
      next=cur->next;
      delete(cur);
      cur=next;
    }
    number_of_items=0;
  //  printf("Liste Zerstoert!\n");
}
//privates Insert
template <class L_DATA>
DLItem<L_DATA> *DLList<L_DATA>::pInsert(L_DATA data, DLItem<L_DATA> *pos)
{
   DLItem<L_DATA> *i=new DLItem<L_DATA>(data, number_of_items+1, pos->previous, pos);
   if (i)
      {
      pos->previous->next=i;
      pos->previous=i;
      number_of_items++;
      return(i);
      }
   else return(0);
}
//privates delete
template <class L_DATA>
L_DATA DLList<L_DATA>::pDelete(DLItem<L_DATA> *i)
{
  L_DATA data=i->item;
  i->previous->next=i->next;
  i->next->previous=i->previous;
//  array[i->index]=0;
  delete(i);
  number_of_items--;
  return(data);
}
//oeffentliches Insert
template <class L_DATA>
int DLList<L_DATA>::Insert(L_DATA data, unsigned long pos)
{
  if ((pos<0)||(pos>(number_of_items))) return(0);
  DLItem<L_DATA> *cur=head;
  while(pos--) cur=cur->next;
  return(pInsert(data,cur)!=0);
}
//oeffentliche Delete
template <class L_DATA>
int DLList<L_DATA>::Delete(unsigned long pos)
{
  if ((pos<0)||(pos>(number_of_items))) return(0);
  DLItem<L_DATA> *cur=head;
  while(pos--) cur=cur->next;
  return(pDelete(cur)!=0);
}

//oeffentliche Delete
template <class L_DATA>
int DLList<L_DATA>::fDelete(L_DATA data)
{
  if ((number_of_items==0) || (!data)) return(0);
  DLItem<L_DATA> *cur;
  cur=head->next;
  while ((cur!=tail) && (cur->item!=data)) cur=cur->next;
  if (cur!=tail) return(pDelete(cur)!=0);
  return(0);
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Push(L_DATA data)
{
  DLItem<L_DATA> *tmp;
  tmp=pInsert(data,tail);
  if (tmp) return (tmp->item);
  return(0);
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Pop(void)
{
  return(pDelete(tail->previous));
}


template <class L_DATA>
L_DATA DLList<L_DATA>::Get(unsigned long pos)
{
  if ((pos<1)||(pos>(number_of_items+1))) return(0);
//  return(array[pos]->item);
  DLItem<L_DATA> *cur=head;
  while(pos--) cur=cur->next;
  return(cur->item);
}


template <class L_DATA>
int DLList<L_DATA>::Enqueue(L_DATA data)
{
  return(pInsert(data,tail)!=0);
}

template <class L_DATA>
L_DATA DLList<L_DATA>::Dequeue(void)
{
  return(pDelete(head->next));
}

//gibt Index des gesuchte Listenelement zurueck, besser waere eigentlich zeiger
template <class L_DATA>
unsigned long DLList<L_DATA>::Is_In_List(L_DATA data)
{
  DLItem<L_DATA> *cur=head, *next;
  unsigned long pos=0;
  while (cur)
    {
      next=cur->next;
      if (cur->item==data) return(pos) ;
      cur=next;
      pos++;
    }
  return(0);
}

//######################################################################################################################
template <class L_DATA>
DL_Indexed_List<L_DATA>::DL_Indexed_List(void) : DLList<L_DATA>()
{
  last_index=0;
}

template <class L_DATA>
DL_Indexed_List<L_DATA>::~DL_Indexed_List()
{
  DLItem<L_DATA> *cur, *next;
  cur=this->head;
  while (cur)
    {
      next=cur->next;
      delete(cur);
      cur=next;
    }
    this->number_of_items=0;
  //  printf("Liste Zerstoert!\n");
}

//privates Insert
template <class L_DATA>
DLItem<L_DATA> *DL_Indexed_List<L_DATA>::pInsert(L_DATA data, DLItem<L_DATA> *pos)
{
    DLItem<L_DATA> *i=new DLItem<L_DATA>(data, last_index, pos->previous, pos);
   if (i)
      {
      pos->previous->next=i;
      pos->previous=i;
      this->number_of_items++;
      array[last_index]=i;
      last_index++;
      return(i);
      }
   else return(0);
}
//privates delete
template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::pDelete(DLItem<L_DATA> *i)
{
  L_DATA data=i->item;
  i->previous->next=i->next;
  i->next->previous=i->previous;
  array[i->index]=0;
  last_index=i->index;
  delete(i);
  this->number_of_items--;
  return(data);
}
template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Push(L_DATA data)
{
  DLItem<L_DATA> *tmp;
  tmp=pInsert(data,this->tail);
  if (tmp) return (tmp->item);
  return(0);
}

template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Pop(void)
{
  return(pDelete(this->tail->previous));
}

template <class L_DATA>
L_DATA DL_Indexed_List<L_DATA>::Get(unsigned long pos)
{
  if ((pos<0)||(pos>(this->number_of_items-1))) return(0);
  return(array[pos]->item);
}

//#######################################################################################

//************************************************************************************************************
template <class L_DATA>
ClusterList<L_DATA>::ClusterList(void) : DLList<L_DATA>()
{
  links_out_of_cluster=0;
  links_inside_cluster=0;
  frequency=1;
  cluster_energy=1e30;
  candidates=new DLList<L_DATA>();
  marker=0;
}

template <class L_DATA>
ClusterList<L_DATA>::~ClusterList() 
{
  while (candidates->Size()) {
    candidates->Pop();
  }
  delete candidates;
}


template <class L_DATA>
bool ClusterList<L_DATA>::operator==(ClusterList<L_DATA> &b)
{
bool found=false;
L_DATA n_cur, n_cur_b;
DLList_Iter<L_DATA> a_iter,b_iter;

if (this->Size()!=b.Size()) return false;

n_cur=a_iter.First(this);
while (!(a_iter.End()))
{
  found=false;
  n_cur_b=b_iter.First(&b);
  while (!(b_iter.End()) && !found)
  {
    if (n_cur==n_cur_b) found=true;
    n_cur_b=b_iter.Next();
  }
  if (!found) return false;
  n_cur=a_iter.Next();
}
return(found);
}
//A<B ist Wahr, wenn A echte Teilmenge von B ist
template <class L_DATA>
bool ClusterList<L_DATA>::operator<(ClusterList<L_DATA> &b)
{
bool found=false;
L_DATA n_cur, n_cur_b;
DLList_Iter<L_DATA> a_iter, b_iter;

if (this->Size()>=b.Size()) return false;
n_cur=a_iter.First(this);
while (!(a_iter.End()))
{
  found=false;
  n_cur_b=b_iter.First(&b);
  while (!(b_iter.End()) && !found)
  {
    if (n_cur==n_cur_b) found=true;
    n_cur_b=b_iter.Next();
  }
  if (!found) return false;
  n_cur=a_iter.Next();
}
return(found);
}
//###############################################################################
//Constructor
NNode::NNode(unsigned long ind, unsigned long c_ind, DLList<NLink*> *ll, char* n, int states)
{
   index=ind;
   cluster_index=c_ind;
   neighbours = new DLList<NNode*>();
   n_links = new DLList<NLink*>();
   global_link_list=ll;
   strcpy(name,n);
   color.red=0;
   color.green=0;
   color.blue=0;
   strcpy(color.pajek_c,"Green");
   clustering=0.0;
   marker=0;
   affiliations=0;
   weight=0.0;
   affinity=0.0;
   distance=0;
   max_states=states;
   state_history=new unsigned long[states+1];
}

//Destructor
NNode::~NNode()
{
  Disconnect_From_All();
  delete neighbours;
  delete n_links;
  delete state_history;
  neighbours=NULL;
  n_links=NULL;
  state_history=NULL;
}

void NNode::Add_StateHistory(unsigned int state)
{
  if (max_states>=state) {
    state_history[state]++;
  }
}

void NNode::Set_Color(RGBcolor c)
{
  color.red=c.red; color.blue=c.blue; color.green=c.green;
  strcpy(color.pajek_c,c.pajek_c);
}

int NNode::Connect_To(NNode* neighbour, double weight)
{
  NLink *link;
  //sollen doppelte Links erlaubt sein??  NEIN
  if (!neighbour) return 0;
  if (!(neighbours->Is_In_List(neighbour)) && (neighbour!=this))
    {
      neighbours->Push(neighbour);        // nachbar hier eintragen
      neighbour->neighbours->Push(this); // diesen knoten beim nachbarn eintragen

      link=new NLink(this,neighbour, weight);        //link erzeugen
      global_link_list->Push(link);                        // in globaler liste eintragen
      n_links->Push(link);                                   // bei diesem Knoten eintragen
      neighbour->n_links->Push(link);                  // beim nachbarn eintragen

      return(1);
    }
   return(0);
}

NLink *NNode::Get_LinkToNeighbour(NNode* neighbour)
{
  DLList_Iter<NLink*> iter;
  NLink *l_cur, *link;
  bool found=false;
  // finde einen bestimmten Link aus der Liste der links eines Knotens
  l_cur=iter.First(n_links);
  while (!iter.End() && !found)
  {
    if (((l_cur->Get_Start()==this) && (l_cur->Get_End()==neighbour)) || ((l_cur->Get_End()==this) && (l_cur->Get_Start()==neighbour)))
    {
      found=true;
      link=l_cur;
    }
    l_cur=iter.Next();
  }
  if (found) return link; else return NULL;
}

int NNode::Disconnect_From(NNode* neighbour)
{
  //sollen doppelte Links erlaubt sein??  s.o.
  if (!neighbours) return 0;
  neighbours->fDelete(neighbour);
  n_links->fDelete(Get_LinkToNeighbour(neighbour));
  neighbour->n_links->fDelete(neighbour->Get_LinkToNeighbour(this));
  neighbour->neighbours->fDelete(this);
  return 1;
}

int NNode::Disconnect_From_All()
{
 int number_of_neighbours=0;
 while (neighbours->Size()) {
      Disconnect_From(neighbours->Pop());
      number_of_neighbours++;
 }
 return(number_of_neighbours) ;
}

/*
int NNode::Disconnect_From_All_Grandchildren()
{
 int n_l=links->Size();
 unsigned long pos=0;
 while ((n_l--)>1) {  //alle bis auf das erste loeschen
      pos=(links->Get(n_l+1))->links->Is_In_List(this);
     // printf("%d %d\n",n_l,pos);
      (links->Get(n_l+1))->links->Delete(pos);
  }
 return(pos) ;
}
*/

double NNode::Get_Links_Among_Neigbours(void)
{
//  long neighbours1, neighbours2;
  double lam=0;
  DLList_Iter<NNode*> iter1, iter2;
//  neighbours1=neighbours->Size();        //so viele Nachbarn hat die Betrachtete Node
  NNode *step1,*step2;
  step1=iter1.First(neighbours);
  while (!iter1.End()) //  for (int n1=1;n1<=neighbours1; n1++)
  {
    //step1=neighbours->Get(n1);
    //neighbours2=step1->neighbours->Size();  //so viele Nachbarn hat der n1-ste Nachbar
    step2=iter2.First(step1->Get_Neighbours());
    while (!iter2.End()) //for (int n2=1;n2<=neighbours2; n2++)
    {
        //step2=step1->neighbours->Get(n2);
        if (step2->Get_Neighbours()->Is_In_List(this)) {lam++;}
        step2=iter2.Next();
    }
  step1=iter1.Next();
  }
  return(lam/2.0);
}


double NNode::Get_Clustering()
{
  double c;
  unsigned long k;
  k=neighbours->Size();
  if (k<=1) return(0);
  c=2.0*Get_Links_Among_Neigbours()/double(k*k-k);
  return(c);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Constructor
NLink::NLink(NNode *s, NNode *e, double w)
{
   start=s;
   end=e;
   weight=w;
   old_weight=0;
   marker=0;
}

//Destructor
NLink::~NLink()
{
  if (start && end) start->Disconnect_From(end);
}
//#####################################################################################
template <class L_DATA>
DLList_Iter<L_DATA>::DLList_Iter()
{
  list=NULL;
  current=NULL;
  end_reached=true;
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::Next(void)
{
  current=current->next;
  if (current==(list->tail)) end_reached=true;
  return(current->item);
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::Previous(void)
{
  current=current->previous;
  if (current==(list->head)) end_reached=true;
  return(current->item);
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::First(DLList<L_DATA> *l)
{
    list=l;
    current=list->head->next;
    if (current==(list->tail)) end_reached=true;
    else end_reached=false;
    return(current->item);
}

template <class L_DATA>
L_DATA DLList_Iter<L_DATA>::Last(DLList<L_DATA> *l)
{
    list=l;
    current=list->tail->previous;
    if (current==(list->head)) end_reached=true;         // falls die List leer ist
    else end_reached=false;
    return(current->item);
}

template <class L_DATA>
bool DLList_Iter<L_DATA>::Swap(DLList_Iter<L_DATA> b)
{
  L_DATA h;
  if (list!=b.list) return false; //elemeten muessen aus der gleichen List stammen
  if (end_reached || b.end_reached) return false;
  h=current->item; current->item=b.current->item; b.current->item=h;
  return true;
} 
//########################################################################################

 void test2(void)
 {
   
   char *test;
   int *test31, inum;
   unsigned long *test2;
   long *test12;
   double *test3, dnum;
   NNode *test4;
   NLink *link;
   char *test_c;
   DLList_Iter<NNode*> *n_iter, *iter;
   DLList_Iter<ClusterList<NNode*>*> *cl_iter;
   DLList_Iter<DLList<NNode*>*> *dl_iter;
    DLList_Iter<NLink*> *l_iter;
   DLList_Iter<char*> *ch_iter;
   DLList_Iter<double*> *dbl_iter;
   DLList_Iter<int*> *i_iter;   
   DLList_Iter<unsigned long*> *int_iter;
   DLList_Iter<long*> *long_iter;   
   DLList_Iter<cluster_join_move*> *m_iter;
   
   HugeArray<int> array;
   HugeArray<char*> c_array;
   HugeArray<double> d_array;
   HugeArray<double> *ptr_d_array;
   HugeArray< HugeArray<double>* > h_array;
   HugeArray<NNode*> n_array;
   HugeArray<NLink*> l_array;    

   HUGE_INDEX h;
   
   h=array.get_huge_index(10);
   h=c_array.get_huge_index(10);
   h=n_array.get_huge_index(10);   
   h=d_array.get_huge_index(10);
   h=h_array.get_huge_index(10);
   h=l_array.get_huge_index(10);

   n_array[1]=test4;
   n_array.Set(1)=test4;
   test4=n_array.Get(1);
   array[1]=1;
   array.Set(1)=1;
   inum=array.Get(1);

   l_array[1]=link;
   l_array.Set(1)=link;
   link=l_array.Get(1);
   
   d_array[1]=1;
   d_array.Set(1)=1;
   dnum=d_array.Get(1);
   ptr_d_array=new HugeArray<double>;
   ptr_d_array->Set(1)=12.4;
   h_array[1]=ptr_d_array;
   h_array[1]->Set(1)=1.0;
   dnum=h_array[1]->Get(1);
   
   c_array[1]=test_c;
   c_array.Set(1)=test_c;
   test_c=c_array[1];
   test_c=c_array.Get(1);

   DLList<cluster_join_move*> *moveset;
   moveset=new DLList<cluster_join_move*>();
   m_iter=new DLList_Iter<cluster_join_move*>();
   cluster_join_move *move;
   moveset->Is_In_List(move);
   moveset->Push(move);
   move=moveset->Pop();
 //  move=moveset->Get(1);
   moveset->Delete(1);
   moveset->fDelete(move);
   move=m_iter->First(moveset);
   move=m_iter->Next();
   delete move;
   delete moveset;

      
   DLList<char*> *t;
   t=new DLList<char*>();
   t->Push(test);
//   t->Get(1);
   t->Pop();
   t->Delete(1);
   ch_iter=new DLList_Iter<char*>();
   test=ch_iter->First(t);
   test=ch_iter->Next();
   delete(t);
   
   DLList<unsigned long*> *t2;
   t2=new DLList<unsigned long*>();
   t2->Push(test2);
//   t2->Get(1);
   t2->Pop();
   int_iter=new DLList_Iter<unsigned long*>();
   test2=int_iter->First(t2);
   test2=int_iter->Next();
   delete t2;

   
   DLList<double*> *t3;
   t3=new DLList<double*>();
   t3->Push(test3);
//   t3->Get(1);
   t3->Pop();
   dbl_iter=new DLList_Iter<double*>();
   test3=dbl_iter->First(t3);
   test3=dbl_iter->Next();
   delete t3;

   DL_Indexed_List<double*> *t32;
   t32=new DL_Indexed_List<double*>();
   t32->Push(test3);
   t32->Get(1);
   t32->Pop();
   dbl_iter=new DLList_Iter<double*>();
   test3=dbl_iter->First(t32);
   test3=dbl_iter->Next();
   delete t32;


   DL_Indexed_List<int*> *t33;
   t33=new DL_Indexed_List<int*>();
   t33->Push(test31);
   t33->Get(1);
   t33->Pop();
   i_iter=new DLList_Iter<int*>();
   test31=i_iter->First(t33);
   test31=i_iter->Next();
   delete t33;

   
   
   DLList<long*> *t12;
   t12=new DLList<long*>();
   t12->Push(test12);
   t12->Get(1);
   t12->Pop();
   long_iter=new DLList_Iter<long*>();
   test12=long_iter->First(t12);
   test12=long_iter->Next();
   delete t12;

   DL_Indexed_List<unsigned long*> *t15;
   t15=new DL_Indexed_List<unsigned long*>();
   unsigned long *test13;
   t15->Push(test13);
   t15->Get(1);
   t15->Pop();
   DLList_Iter<unsigned long*> *ulong_iter;
   ulong_iter=new DLList_Iter<unsigned long*>();
   test13=ulong_iter->First(t15);
   test13=ulong_iter->Next();
   delete t15;

   DL_Indexed_List<unsigned int*> *t16;
   t16=new DL_Indexed_List<unsigned int*>();
   unsigned int *test17;
   t16->Push(test17);
   t16->Get(1);
   t16->Pop();
   DLList_Iter<unsigned int*> *uint_iter;
   uint_iter=new DLList_Iter<unsigned int*>();
   test17=uint_iter->First(t16);
   test17=uint_iter->Next();
   delete t16;
   
      
   DLList<NNode*> *t4;
   t4=new DLList<NNode*>();
   iter=new DLList_Iter<NNode*>();
   test4=iter->First(t4);
   test4=iter->Last(t4);   
   test4=iter->Next();
   test4=iter->Previous();
   iter->Swap(*iter); 
   t4->Push(test4);
   t4->Push(test4);
   t4->Get(1);
   t4->Pop();
   t4->Delete(1);
   t4->Enqueue(test4);
   test4=t4->Dequeue();
   t4->Insert(test4,1);
   t4->fDelete(test4);
   delete(t4);

   DL_Indexed_List<NNode*> *t42;
   t42=new DL_Indexed_List<NNode*>();
   t42->Push(test4);
   t42->Get(1);
   t42->Pop();
   delete(t42);

   
   DLList<NLink*> *t5;
   t5=new DLList<NLink*>();
   link=new NLink(test4,test4,1);
   t5->Push(link);
   link=t5->Pop();
   link=t5->Get(1);
   l_iter=new DLList_Iter<NLink*>();
   link=l_iter->First(t5);
   link=l_iter->Last(t5);   
   link=l_iter->Next();
   link=l_iter->Previous();
   l_iter->Swap(*l_iter);
   delete t5;

   DL_Indexed_List<NLink*> *t52;
   t52=new DL_Indexed_List<NLink*>();
   t52->Push(link);
   link=t52->Pop();
   link=t52->Get(1);
   l_iter=new DLList_Iter<NLink*>();
   link=l_iter->First(t52);
   link=l_iter->Last(t52);
   link=l_iter->Next();
   link=l_iter->Previous();
   l_iter->Swap(*l_iter);
   delete t52;
   

   DLList<DLList<NNode*>*> *t6;
   t6=new DLList<DLList<NNode*>*>();
   dl_iter = new DLList_Iter<DLList<NNode*>*>();
   t4=dl_iter->First(t6);
   t4=dl_iter->Next();
   t6->Push(t4);
   t6->Get(1);
   t6->Pop();
   delete t6;
         
   DLList<ClusterList<NNode*>*> *t7;
   t7=new DLList<ClusterList<NNode*>*>();
   cl_iter = new DLList_Iter<ClusterList<NNode*>*>();
   ClusterList<NNode*> *t71;
   n_iter = new DLList_Iter<NNode*>();
   test4=n_iter->First(t71);
   test4=n_iter->Next();
   t7->Insert(t71,1);
   delete t7;

   DL_Indexed_List<ClusterList<NNode*>*> *t72;
   t72=new DL_Indexed_List<ClusterList<NNode*>*>();
   cl_iter = new DLList_Iter<ClusterList<NNode*>*>();
   t72->Push(t71);
   t71=t72->Pop();
   t71=t72->Get(1);
   delete t72;
   

   t71=new ClusterList<NNode*>();
   if (*t71<*t71) printf("Baeh");
   if (*t71==*t71) printf("Buuh");
   t71=cl_iter->First(t7);
   t71=cl_iter->Next();
   t71->Insert(test4,1);
   delete t71;
      
   t71->Push(test4);
   t71->fDelete(test4);
   t7->Push(t71);
   t7->Is_In_List(t71);
   t7->Get(1);
   t7->Pop();
   t7->Delete(0);
   t7->fDelete(t71);

 
    DLList_Iter<ClusterList<NNode*>*> *c_iter;
     c_iter=new DLList_Iter<ClusterList<NNode*>*>();
     
  
   printf("Test2\n");
 }


