<?xml version="1.0" encoding="UTF-8" ?>
<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema">

  <xsd:annotation>
    <xsd:documentation xml:lang="en">
      Nexus network repository, dataset list
      http://nexus.igraph.org 
    </xsd:documentation>
  </xsd:annotation>

  <xsd:element name="datasets" type="datasetsType" />

  <xsd:complexType name="datasetsType">
    <xsd:sequence>
      <xsd:element name="size" type="xsd:nonNegativeInteger" />
      <xsd:element name="totalsize" type="xsd:nonNegativeInteger" />
      <xsd:element name="offset" type="xsd:nonNegativeInteger" />
      <xsd:element name="limit" type="xsd:nonNegativeInteger" />
      <xsd:element name="dataset" type="datasetSummaryType" minOccurs="0" 
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>

  <xsd:complexType name="datasetSummaryType">
    <xsd:sequence>
      <xsd:element name="id" type="xsd:positiveInteger" />
      <xsd:element name="name" type="xsd:string" />
      <xsd:element name="description" type="xsd:string" />
      <xsd:element name="vertices" type="xsd:nonNegativeInteger" />
      <xsd:element name="edges" type="xsd:nonNegativeInteger" />
      <xsd:element name="tags" type="tagsType" />
      <xsd:element name="date" type="xsd:date" />      
    </xsd:sequence>
  </xsd:complexType>

  <xsd:complexType name="tagsType">
    <xsd:sequence>
      <xsd:element name="tag" type="xsd:string" minOccurs="0" 
		   maxOccurs="unbounded" />
    </xsd:sequence>
  </xsd:complexType>

</xsd:schema>
