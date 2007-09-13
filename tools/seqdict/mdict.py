################################################################################
#  Sequential Dictionary Class                                                 #
#                                                                              #
#  by Wolfgang Grafen                                                          #
#                                                                              #
#  Version 0.2    11. February 2004
#                                                                              #
#  email to: WolfgangGrafen@gmx.de                                             #
#                                                                              #
################################################################################

from ndict import seqdict   #Sequential Single Value Dictionary
from UserList import UserList

class MyUserList(UserList):
    from UserList import UserList
    def __init__(self,parent,liste=None):
      UserList.__init__(self,liste)
      self.parent = parent  #remember parent for call-back
    def __delitem__(self, i):
      del self.data[i]
      if self.data==[]:     #call-back, deletes item of parent
        index = self.parent.values().index([])
        del self.parent[index:index+1]
        
class mseqdict(seqdict):    #Sequential Multiple Value Dictionary
  def __init__(self,List=[],Dict={}):
    self.list = []
    self.dict = {}
    if not List:
      pass     
    elif type(List)==type({}):
      for key,value in List.items():
        self.__setitem__(key,value)
    elif List and not Dict: #dict.items()
      for key,value in List:
        if isinstance(value,MyUserList):
          for v in value:
            self.__setitem__(key,v)
        else:
          self.__setitem__(key,value)     
    elif type(List)==type(Dict)==type([]):
      for key,value in map(None,List,Dict):
        self.__setitem__(key,value)
    else:
      if isinstance(Dict.values()[0],MyUserList):
        self.dict = Dict
        self.list = List
      else:
        for key in List:
          value = Dict[key]
          if type(value)==type([]):
            for v in value:
              self.__setitem__(key,v)
          else:
            self.__setitem__(key,value)     

    self_list = self.list
    self_dict = self.dict
    for k in self_list:
        assert self_dict.has_key(k),"key %r not in self.dict" % k

    for k in self_dict.keys():
        if k not in self_list:
            self_list.append(k)

  def __setitem__(self,key,value):
    if not self.dict.has_key(key):
      self.list.append(key)
      if isinstance(value,MyUserList):
        self.dict[key] = value
      else:
        self.dict[key]=MyUserList(self,[value])
    else:
      values = self.dict[key]
      if isinstance(value,MyUserList):
        for v in value:
          if not v in values:
            values.extend(MyUserList(self,[v]))
      else:
        #if not value in values:
        for v in values:
            if v is value:
                break
            values.extend(MyUserList(self,[value]))
  def __delitem__(self, key):
    del self.dict[key]
    self.list.remove(key)      

  def append(self,key,value):
    self.__setitem__(key,value)                 
  def __setslice__(self,start,stop,newdict):
    start = max(start,0); stop = max(stop,0)
    delindexes = []
    for key in newdict.keys():
      if self.dict.has_key(key):
        index = self.list.index(key)
        delindexes.append(index)
        if index < start:
          start = start - 1
          stop  = stop  - 1
        elif index >= stop:
          pass
        else:
          stop  = stop - 1
      else:
        self.dict[key]=UserList(self)
    delindexes.sort()
    delindexes.reverse()
    for index in delindexes:
      key = self.list[index]
      #del   self.dict[key]
      del  self.list[index]
    self.list[start:stop] = newdict.list[:]
    self.dict.update(newdict.dict)  
  def copy(self):
    values = map(lambda x:x[:],self.values())
    return self.__class__(self.list,values)
  def count(self,value):
    vallist = self.dict.values()
    return map(lambda x,y=value:x.count(y),vallist).count(1)
  def filter(self,function,filtervalues=0):
    if   filtervalues == 1: #consider key and all keyvalues at once
      dict = self.__class__()
      for key,values in self.items():
        if function(key,values):
          dict[key]=values
      return dict
    elif filtervalues == 2: #consider key and every keyvalue for itself
      dict = self.__class__()
      for key,values in self.items():
        for value in values:
          if function(key,value):
            dict[key]=value
      return dict
    else:                   #consider key only
      liste=filter(function,self.list)
      dict = {}
      for i in liste:
        dict[i]=self.dict[i]
      return self.__class__(liste,dict)
  def map(self,function,mapvalues=2):
    if   mapvalues == 1:    #consider key and all keyvalues at once
      dict = self.__class__()
      for key,values in self.items():
        k,v = function(key,values)
        dict[k]=v
      return dict
    else: #if mapvalues!=1: #consider key and every keyvalue for itself
      dict = self.__class__()
      for key,values in self.items():
        for value in values:
          k,v = function(key,value)
          dict[k]=v
      return dict
  def pop(self,key='...None',value='...None'):
    if value=='...None':
      if key=='...None':
        pos = -1
        key = self.list[pos]
      else:
        pos = self.list.index(key)
      tmp = self.dict[key]
      del self.dict[key]
      return {self.list.pop(pos):tmp}
    else:
      val = self.dict[key]
      index = val.index(value)
      tmp = val[index]
      del val[index]
      return {key:tmp}
  def remove(self,key,value='...None'):
    if value=='...None':
      del self[key]
    else:
      index = self[key].index(value) 
      del self[key][index]
  def sort(self,func1=None,func2=None):
    if not func1:
      self.list.sort()
    else:
      apply(self.list.sort,[func1])
    if func2:
      for value in self.values():
        apply(value.sort,[func2])
      
  def swap(self):
    tmp = self.__class__()
    for key,values in self.items():
      for value in values:
        tmp[value]=key
    self.list,self.dict = tmp.list,tmp.dict
    del tmp
  
  def __repr__(self):return 'mseqdict(\n%s,\n%s)'%(self.list,self.dict)
